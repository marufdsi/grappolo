// ***********************************************************************
//
//            Grappolo: A C++ library for graph clustering
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2014) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include "defs.h"
#include "basic_comm.h"
#include "basic_util.h"
#include <sstream>
#include <fstream>
#include <sys/stat.h>
using namespace std;

std::vector<std::string> split(std::string str, char delim) {
    std::stringstream ss(str);
    std::string token;
    std::vector<std::string>tokens;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

//WARNING: This will overwrite the original graph data structure to
//         minimize memory footprint
// Return: C_orig will hold the cluster ids for vertices in the original graph
//         Assume C_orig is initialized appropriately
//WARNING: Graph G will be destroyed at the end of this routine
void runMultiPhaseBasic(graph *G, long *C_orig, int basicOpt, long minGraphSize,
                        double threshold, double C_threshold, int numThreads, int threadsOpt, char *graphName) {
    double totTimeClustering = 0, totTimeBuildingPhase = 0, totTimeColoring = 0, tmpTime = 0;
    int tmpItr = 0, totItr = 0;
    long NV = G->numVertices;


    /* Step 1: Find communities */
    f_weight prevMod = -1;
    f_weight currMod = -1;
    long phase = 1;

    graph *Gnew; //To build new hierarchical graphs
    long numClusters;
    long *C = (long *) malloc(NV * sizeof(long));
    assert(C != 0);
#pragma omp parallel for
    for (long i = 0; i < NV; i++) {
        C[i] = -1;
    }

    while (1) {
        printf("===============================\n");
        printf("Phase %ld\n", phase);
        printf("===============================\n");
        prevMod = currMod;

        bool change = false;
        if (basicOpt == 1) {
            currMod = parallelLouvianMethodNoMap(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else if (threadsOpt == 1) {
            currMod = vectorizedLouvianMethod(G, C, numThreads, currMod, (f_weight)threshold, &tmpTime, &tmpItr, &change);
            //currMod = parallelLouvianMethodApprox(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else {
            currMod = parallelLouvianMethodScale(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }

        totTimeClustering += tmpTime;
        totItr += tmpItr;

        //Renumber the clusters contiguiously
        numClusters = renumberClustersContiguously(C, G->numVertices);
        printf("Number of unique clusters: %ld\n", numClusters);

        //printf("About to update C_orig\n");
        //Keep track of clusters in C_orig
        if (phase == 1) {
#pragma omp parallel for
            for (long i = 0; i < NV; i++) {
                C_orig[i] = C[i]; //After the first phase
            }
        } else {
#pragma omp parallel for
            for (long i = 0; i < NV; i++) {
                assert(C_orig[i] < G->numVertices);
                if (C_orig[i] >= 0)
                    C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
            }
        }
        printf("Done updating C_orig\n");

        //Break if too many phases or iterations
        if ((phase > 200) || (totItr > 100000)) {
            break;
        }

        //Check for modularity gain and build the graph for next phase
        //In case coloring is used, make sure the non-coloring routine is run at least once
        if (change /*(currMod - prevMod) > threshold*/) {
            Gnew = (graph *) malloc(sizeof(graph));
            assert(Gnew != 0);
            tmpTime = buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
            totTimeBuildingPhase += tmpTime;
            //Free up the previous graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew; //Swap the pointers
            G->edgeListPtrs = Gnew->edgeListPtrs;
            G->edgeList = Gnew->edgeList;

            //Free up the previous cluster & create new one of a different size
            free(C);
            C = (long *) malloc(numClusters * sizeof(long));
            assert(C != 0);

#pragma omp parallel for
            for (long i = 0; i < numClusters; i++) {
                C[i] = -1;
            }
            phase++; //Increment phase number
        } else {
            break; //Modularity gain is not enough. Exit.
        }

    } //End of while(1)
    std::vector<std::string> parts = split(std::string(graphName), '/');
    std::ofstream resultCSV;
    std::string folderName = "Results/";
    std::string fileName = "Grappolo_Lovain_Result_SFP.csv";
    if (mkdir(folderName.c_str(), 0777) == -1)
        std::cout << "Directory " << folderName << " is already exist" << std::endl;
    else
        std::cout << "Directory " << folderName << " created" << std::endl;
    std::ifstream infile(folderName + fileName);
    resultCSV.open(folderName + fileName, std::ios_base::out | std::ios_base::app | std::ios_base::ate);

    if (!infile.good()) {
        resultCSV
                << "GraphName,Threads,Phases,TotalIterations,Clusters,Modularity,ClusteringTIme,CoarseningTime,TotalTime,Threshold,DataType"
                << std::endl;
    }
    infile.close();
    resultCSV << split(parts[parts.size()-1], '.')[0] << "," << numThreads << "," << phase << "," << totItr << "," << numClusters << "," << prevMod
              << "," << totTimeClustering << "," << totTimeBuildingPhase << ","
              << totTimeClustering + totTimeBuildingPhase + totTimeColoring << "," << threshold << "," << sizeof(f_weight) << std::endl;
    resultCSV.close();
    printf("********************************************\n");
    printf("*********    Compact Summary   *************\n");
    printf("********************************************\n");
    printf("Number of threads              : %ld\n", numThreads);
    printf("Total number of phases         : %ld\n", phase);
    printf("Total number of iterations     : %ld\n", totItr);
    printf("Final number of clusters       : %ld\n", numClusters);
    printf("Final modularity               : %lf\n", prevMod);
    printf("Total time for clustering      : %lf\n", totTimeClustering);
    printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
    printf("********************************************\n");
    printf("TOTAL TIME                     : %lf\n", (totTimeClustering + totTimeBuildingPhase + totTimeColoring));
    printf("********************************************\n");

    //Clean up:
    free(C);
    if (G != 0) {
        free(G->edgeListPtrs);
        free(G->edgeList);
        free(G);
    }
}//End of runMultiPhaseLouvainAlgorithm()

//WARNING: This will overwrite the original graph data structure to
//         minimize memory footprint
// Return: C_orig will hold the cluster ids for vertices in the original graph
//         Assume C_orig is initialized appropriately
//WARNING: Graph G will be destroyed at the end of this routine
void runMultiPhaseBasic_sfp(graph *G, long *C_orig, int basicOpt, long minGraphSize,
                        f_weight threshold, f_weight C_threshold, int numThreads, int threadsOpt, char *graphName) {
    printf("single floating point version called, threshold: %f \n", threshold);
    double totTimeClustering = 0, totTimeBuildingPhase = 0, totTimeColoring = 0, tmpTime = 0;
    int tmpItr = 0, totItr = 0;
    long NV = G->numVertices;


    /* Step 1: Find communities */
    f_weight prevMod = -1;
    f_weight currMod = -1;
    long phase = 1;

    graph *Gnew; //To build new hierarchical graphs
    long numClusters;
    long *C = (long *) malloc(NV * sizeof(long));
    assert(C != 0);
#pragma omp parallel for
    for (long i = 0; i < NV; i++) {
        C[i] = -1;
    }

    while (1) {
        printf("===============================\n");
        printf("Phase %ld\n", phase);
        printf("===============================\n");
        prevMod = currMod;

        bool change = false;
        if (basicOpt == 1) {
            currMod = parallelLouvianMethodNoMap(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else if (threadsOpt == 1) {
            currMod = vectorizedLouvianMethod(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr, &change);
            //currMod = parallelLouvianMethodApprox(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else {
            currMod = parallelLouvianMethodScale(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }

        totTimeClustering += tmpTime;
        totItr += tmpItr;

        //Renumber the clusters contiguiously
        numClusters = renumberClustersContiguously(C, G->numVertices);
        printf("Number of unique clusters: %ld\n", numClusters);

        //printf("About to update C_orig\n");
        //Keep track of clusters in C_orig
        if (phase == 1) {
#pragma omp parallel for
            for (long i = 0; i < NV; i++) {
                C_orig[i] = C[i]; //After the first phase
            }
        } else {
#pragma omp parallel for
            for (long i = 0; i < NV; i++) {
                assert(C_orig[i] < G->numVertices);
                if (C_orig[i] >= 0)
                    C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
            }
        }
        printf("Done updating C_orig\n");

        //Break if too many phases or iterations
        if ((phase > 200) || (totItr > 100000)) {
            break;
        }

        //Check for modularity gain and build the graph for next phase
        //In case coloring is used, make sure the non-coloring routine is run at least once
        if (change /*(currMod - prevMod) > threshold*/) {
            Gnew = (graph *) malloc(sizeof(graph));
            assert(Gnew != 0);
            tmpTime = buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
            totTimeBuildingPhase += tmpTime;
            //Free up the previous graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew; //Swap the pointers
            G->edgeListPtrs = Gnew->edgeListPtrs;
            G->edgeList = Gnew->edgeList;

            //Free up the previous cluster & create new one of a different size
            free(C);
            C = (long *) malloc(numClusters * sizeof(long));
            assert(C != 0);

#pragma omp parallel for
            for (long i = 0; i < numClusters; i++) {
                C[i] = -1;
            }
            phase++; //Increment phase number
        } else {
            break; //Modularity gain is not enough. Exit.
        }

    } //End of while(1)
    std::vector<std::string> parts = split(std::string(graphName), '/');
    std::ofstream resultCSV;
    std::string folderName = "Results/";
    std::string fileName = "Grappolo_Lovain_Result_SFP.csv";
    if (mkdir(folderName.c_str(), 0777) == -1)
        std::cout << "Directory " << folderName << " is already exist" << std::endl;
    else
        std::cout << "Directory " << folderName << " created" << std::endl;
    std::ifstream infile(folderName + fileName);
    resultCSV.open(folderName + fileName, std::ios_base::out | std::ios_base::app | std::ios_base::ate);

    if (!infile.good()) {
        resultCSV
                << "GraphName,Threads,Phases,TotalIterations,Clusters,Modularity,ClusteringTIme,CoarseningTime,TotalTime,Threshold,DataType"
                << std::endl;
    }
    infile.close();
    resultCSV << split(parts[parts.size()-1], '.')[0] << "," << numThreads << "," << phase << "," << totItr << "," << numClusters << "," << prevMod
              << "," << totTimeClustering << "," << totTimeBuildingPhase << ","
              << totTimeClustering + totTimeBuildingPhase + totTimeColoring << "," << threshold << "," << sizeof(f_weight) << std::endl;
    resultCSV.close();
    printf("********************************************\n");
    printf("*********    Compact Summary   *************\n");
    printf("********************************************\n");
    printf("Number of threads              : %ld\n", numThreads);
    printf("Total number of phases         : %ld\n", phase);
    printf("Total number of iterations     : %ld\n", totItr);
    printf("Final number of clusters       : %ld\n", numClusters);
    printf("Final modularity               : %lf\n", prevMod);
    printf("Total time for clustering      : %lf\n", totTimeClustering);
    printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
    printf("********************************************\n");
    printf("TOTAL TIME                     : %lf\n", (totTimeClustering + totTimeBuildingPhase + totTimeColoring));
    printf("********************************************\n");

    //Clean up:
    free(C);
    if (G != 0) {
        free(G->edgeListPtrs);
        free(G->edgeList);
        free(G);
    }
}//End of runMultiPhaseLouvainAlgorithm()

// run one phase of Louvain and return modularity
void runMultiPhaseBasicOnce(graph *G, long *C_orig, int basicOpt, long minGraphSize,
                            double threshold, double C_threshold, int numThreads, int threadsOpt) {
    double totTimeClustering = 0, totTimeBuildingPhase = 0, totTimeColoring = 0, tmpTime = 0;
    int tmpItr = 0, totItr = 0;
    long NV = G->numVertices;

    /* Step 1: Find communities */
    double prevMod = -1;
    double currMod = -1;

    graph *Gnew; //To build new hierarchical graphs
    long numClusters;
    long *C = (long *) malloc(NV * sizeof(long));
    assert(C != 0);
#pragma omp parallel for
    for (long i = 0; i < NV; i++) {
        C[i] = -1;
    }

    // Run just one phase
    {
        bool change = false;
        prevMod = currMod;

        if (basicOpt == 1) {
            currMod = parallelLouvianMethodNoMap(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else if (threadsOpt == 1) {
            currMod = parallelLouvianMethod(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr, &change);
            //currMod = parallelLouvianMethodApprox(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        } else {
            currMod = parallelLouvianMethodScale(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }

        totTimeClustering += tmpTime;
        totItr += tmpItr;

        //Renumber the clusters contiguiously
        numClusters = renumberClustersContiguously(C, G->numVertices);
        printf("Number of unique clusters: %ld\n", numClusters);

        //Keep track of clusters in C_orig
#pragma omp parallel for
        for (long i = 0; i < NV; i++) {
            C_orig[i] = C[i]; //After the first phase
        }
        printf("Done updating C_orig\n");

        //Check for modularity gain and build the graph for next phase
        //In case coloring is used, make sure the non-coloring routine is run at least once
        if (change /*(currMod - prevMod) > threshold*/) {
            Gnew = (graph *) malloc(sizeof(graph));
            assert(Gnew != 0);
            tmpTime = buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
            totTimeBuildingPhase += tmpTime;
            //Free up the previous graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew; //Swap the pointers
            G->edgeListPtrs = Gnew->edgeListPtrs;
            G->edgeList = Gnew->edgeList;

            //Free up the previous cluster & create new one of a different size
            free(C);
            C = (long *) malloc(numClusters * sizeof(long));
            assert(C != 0);

#pragma omp parallel for
            for (long i = 0; i < numClusters; i++) {
                C[i] = -1;
            }
        }

    } //End of while(1)

    printf("********************************************\n");
    printf("***********    After Phase 1   *************\n");
    printf("********************************************\n");
    printf("Number of threads              : %ld\n", numThreads);
    printf("Total number of iterations     : %ld\n", totItr);
    printf("Final number of clusters       : %ld\n", numClusters);
    printf("Final modularity               : %lf\n", currMod);
    printf("Total time for clustering      : %lf\n", totTimeClustering);
    printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
    printf("********************************************\n");
    printf("TOTAL TIME                     : %lf\n", (totTimeClustering + totTimeBuildingPhase + totTimeColoring));
    printf("********************************************\n");

    //Clean up:
    free(C);
    if (G != 0) {
        free(G->edgeListPtrs);
        free(G->edgeList);
        free(G);
    }

}//End of runMultiPhaseLouvainAlgorithm()
