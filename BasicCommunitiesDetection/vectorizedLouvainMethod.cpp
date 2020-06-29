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
#include "utilityClusteringFunctions.h"
#include "basic_comm.h"
#include "../DefineStructure/defs.h"
#include <unistd.h>
using namespace std;
f_weight parallelLouvianMethod_SFP(graph *G, comm_type *C, int nThreads, f_weight Lower,
                               f_weight thresh, double *totTime, int *numItr, bool *change) {

    cout << "modified parallel version called" << endl;
    //    size_t alignment = sysconf(_SC_PAGESIZE);
    size_t alignment = 512;
#ifdef PRINT_DETAILED_STATS_
    printf("Within parallelLouvianMethod()\n");
#endif
    if (nThreads < 1)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(nThreads);
    int nT;
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
#ifdef PRINT_DETAILED_STATS_
    printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
#endif


    double time1, time2, time3, time4; //For timing purposes
    double total = 0, totItr = 0;

    comm_type    NV        = G->numVertices;
    comm_type    NS        = G->sVertices;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;

    comm_type nnz = 0;
    for (comm_type i=0; i<NV; i++) {
        comm_type adj1 = vtxPtr[i];        //Begin
        comm_type adj2 = vtxPtr[i + 1];    //End
        for (comm_type j = adj1; j < adj2; j++) {
            nnz++;
        }
    }
    cout << "NNZ: " << nnz << " NE: " << NE << " twice NE: " << (2*NE) << endl;
    /* Variables for computing modularity */
    comm_type totalEdgeWeightTwice;
    f_weight constantForSecondTerm;
    f_weight prevMod=-1;
    f_weight currMod=-1;
    //double thresMod = 0.000001;
    f_weight thresMod = thresh; //Input parameter
    int numItrs = 0;

    /********************** Initialization **************************/
    time1 = omp_get_wtime();
    //Store the degree of all vertices
    f_weight * vDegree; // = (f_weight *) malloc (NV * sizeof(f_weight));
    posix_memalign((void **) &vDegree, alignment, NV * sizeof(f_weight));
    assert(vDegree != 0);

    //Community info. (ai and size)
    /*Comm *cInfo; // = (Comm *) malloc (NV * sizeof(Comm));
    posix_memalign((void **) &cInfo, alignment, NV * sizeof(Comm));
    assert(cInfo != 0);*/
    /// replace cInfo by the following cInfo_size and cInfo_degree
    comm_type* cInfo_size;
    posix_memalign((void **) &cInfo_size, alignment, NV * sizeof(comm_type));
    assert(cInfo_size != 0);
    f_weight* cInfo_degree;
    posix_memalign((void **) &cInfo_degree, alignment, NV * sizeof(f_weight));
    assert(cInfo_degree != 0);

    //use for updating Community
    /*Comm *cUpdate; // = (Comm*)malloc(NV*sizeof(Comm));
    posix_memalign((void **) &cUpdate, alignment, NV * sizeof(Comm));
    assert(cUpdate != 0);*/
    /// replace cUpdate by the following cUpdate_size and cUpdate_degree
    comm_type* cUpdate_size;
    posix_memalign((void **) &cUpdate_size, alignment, NV * sizeof(comm_type));
    assert(cUpdate_size != 0);
    f_weight* cUpdate_degree;
    posix_memalign((void **) &cUpdate_degree, alignment, NV * sizeof(f_weight));
    assert(cUpdate_degree != 0);

    //use for Modularity calculation (eii)
    f_weight* clusterWeightInternal; // = (f_weight*) malloc (NV*sizeof(f_weight));
    posix_memalign((void **) &clusterWeightInternal, alignment, NV * sizeof(f_weight));
    assert(clusterWeightInternal != 0);

    sumVertexDegreeVec_sfp(vtxInd, vtxPtr, vDegree, NV , cInfo_size, cInfo_degree);	// Sum up the vertex degree

    /*** Compute the total edge weight (2m) and 1/2m ***/
    constantForSecondTerm = calConstantForSecondTerm_sfp(vDegree, NV); // 1 over sum of the degree

    //Community assignments:
    //Store previous iteration's community assignment
    comm_type* pastCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &pastCommAss, alignment, NV * sizeof(comm_type));
    assert(pastCommAss != 0);
    //Store current community assignment
    comm_type* currCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &currCommAss, alignment, NV * sizeof(comm_type));
    assert(currCommAss != 0);
    //Store the target of community assignment
    comm_type* targetCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &targetCommAss, alignment, NV * sizeof(comm_type));
    assert(targetCommAss != 0);

    comm_type* head;
    posix_memalign((void **) &head, alignment, (nnz * sizeof(comm_type)));
    assert(head != 0);
    comm_type* tail;
    posix_memalign((void **) &tail, alignment, (nnz * sizeof(comm_type)));
    assert(tail != 0);
    f_weight* weights;
    posix_memalign((void **) &weights, alignment, (nnz * sizeof(f_weight)));
    assert(weights != 0);
    for (int i = 0; i < nnz; ++i) {
        head[i] = vtxInd[i].head;
        tail[i] = vtxInd[i].tail;
        weights[i] = vtxInd[i].weight;
    }

    //Vectors used in place of maps: Total size = |V|+2*|E| -- The |V| part takes care of self loop
    /*mapElement* clusterLocalMap; // = (mapElement *) malloc ((NV + 2*NE) * sizeof(mapElement));
    posix_memalign((void **) &clusterLocalMap, alignment, ((NV + 2*NE) * sizeof(mapElement)));
    assert(clusterLocalMap != 0);*/
    /// Replace clusterLocalMap by the following cid and Counter
    comm_type* cid;
    posix_memalign((void **) &cid, alignment, ((NV + 2*NE) * sizeof(comm_type)));
    assert(cid != 0);
    f_weight* Counter;
    posix_memalign((void **) &Counter, alignment, ((NV + 2*NE) * sizeof(f_weight)));
    assert(Counter != 0);
    //double* Counter             = (double *)     malloc ((NV + 2*NE) * sizeof(double));     assert(Counter != 0);

    //Initialize each vertex to its own cluster
    //initCommAssOpt(pastCommAss, currCommAss, NV, clusterLocalMapX, vtxPtr, vtxInd, cInfo, constantForSecondTerm, vDegree);

    //Initialize each vertex to its own cluster
//    initCommAss_SFP(pastCommAss, currCommAss, NV);
    initCommAssOptVec_SFP(pastCommAss, currCommAss, NV, cid, Counter, vtxPtr, head, tail, weights, cInfo_size, cInfo_degree, constantForSecondTerm, vDegree);
    time2 = omp_get_wtime();
    printf("Time to initialize: %3.3lf\n", time2-time1);

#ifdef PRINT_DETAILED_STATS_
    printf("========================================================================================================\n");
    printf("Itr      E_xx            A_x2           Curr-Mod         Time-1(s)       Time-2(s)        T/Itr(s)\n");
    printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
    printf("=====================================================\n");
    printf("Itr      Curr-Mod         T/Itr(s)      T-Cumulative\n");
    printf("=====================================================\n");
#endif
    //Start maximizing modularity
    while(true) {
        numItrs++;
        time1 = omp_get_wtime();
        /* Re-initialize datastructures */
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            clusterWeightInternal[i] = 0;
            cUpdate_degree[i] =0;
            cUpdate_size[i] =0;
        }

        bool moved = false;
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            comm_type adj1 = vtxPtr[i];
            comm_type adj2 = vtxPtr[i+1];
            f_weight selfLoop = 0;
            //Build a datastructure to hold the cluster structure of its neighbors
//            map<comm_type, comm_type> clusterLocalMap; //Map each neighbor's cluster to a local number
//            map<comm_type, comm_type>::iterator storedAlready;
//            vector<f_weight> Counter; //Number of edges in each unique cluster
            comm_type numUniqueClusters = 0;
            //Add v's current cluster:
            if(adj1 != adj2){
//                clusterLocalMap[currCommAss[i]] = 0;
//                Counter.push_back(0); //Initialize the counter to ZERO (no edges incident yet)
                comm_type sPosition = vtxPtr[i]+i; //Starting position of local map for i
                Counter[sPosition] = 0;          //Initialize the counter to ZERO (no edges incident yet)
                cid[sPosition] = currCommAss[i]; //Initialize with current community
                numUniqueClusters++; //Added the first entry

                //Find unique cluster ids and #of edges incident (eicj) to them
//                selfLoop = buildLocalMapCounter_sfp(adj1, adj2, clusterLocalMap, Counter, vtxInd, currCommAss, i);
                selfLoop = buildLocalMapCounterNoMap_SFP(i, cid, Counter, vtxPtr, head, tail, weights, currCommAss, numUniqueClusters);
                // Update delta Q calculation
//                clusterWeightInternal[i] += Counter[0]; //(e_ix)
                clusterWeightInternal[i] += Counter[sPosition]; //(e_ix)
                //Calculate the max
//                targetCommAss[i] = max_sfp(clusterLocalMap, Counter, selfLoop, cInfo, vDegree[i], currCommAss[i], constantForSecondTerm);
                targetCommAss[i] = maxNoMap_SFP(i, cid, Counter, vtxPtr, selfLoop, cInfo_size, cInfo_degree, vDegree[i],
                                                   currCommAss[i], constantForSecondTerm, numUniqueClusters);
                //assert((targetCommAss[i] >= 0)&&(targetCommAss[i] < NV));
            } else {
                targetCommAss[i] = -1;
            }

            //Update
            if(targetCommAss[i] != currCommAss[i]  && targetCommAss[i] != -1) {
                moved = true;
#pragma omp atomic update
                cUpdate_degree[targetCommAss[i]] += vDegree[i];
#pragma omp atomic update
                cUpdate_size[targetCommAss[i]] += 1;
#pragma omp atomic update
                cUpdate_degree[currCommAss[i]] -= vDegree[i];
#pragma omp atomic update
                cUpdate_size[currCommAss[i]] -=1;
                /*
                 __sync_fetch_and_add(&cUpdate[targetCommAss[i]].size, 1);
                 __sync_fetch_and_sub(&cUpdate[currCommAss[i]].degree, vDegree[i]);
                 __sync_fetch_and_sub(&cUpdate[currCommAss[i]].size, 1);*/
            }//End of If()
        }//End of for(i)
        time2 = omp_get_wtime();

        time3 = omp_get_wtime();
        f_weight e_xx = 0;
        f_weight a2_x = 0;

#pragma omp parallel for \
reduction(+:e_xx) reduction(+:a2_x)
        for (comm_type i=0; i<NV; i++) {
            e_xx += clusterWeightInternal[i];
            a2_x += (cInfo_degree[i])*(cInfo_degree[i]);
        }
        time4 = omp_get_wtime();

        currMod = (e_xx*(f_weight)constantForSecondTerm) - (a2_x*(f_weight)constantForSecondTerm*(f_weight)constantForSecondTerm);
        totItr = (time2-time1) + (time4-time3);
        total += totItr;
#ifdef PRINT_DETAILED_STATS_
        printf("%d \t %g \t %g \t %lf \t %3.3lf \t %3.3lf  \t %3.3lf\n",numItrs, e_xx, a2_x, currMod, (time2-time1), (time4-time3), totItr );
#endif
#ifdef PRINT_TERSE_STATS_
        printf("%d \t %lf \t %3.3lf  \t %3.3lf\n",numItrs, currMod, totItr, total);
#endif

        if(moved)
            *change = true;

        //Break if modularity gain is not sufficient
        if(!moved || numItrs >= 25/*(currMod - prevMod) < thresMod*/) {
            break;
        }

        //Else update information for the next iteration
        prevMod = currMod;
        if(prevMod < Lower)
            prevMod = Lower;
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            cInfo_size[i] += cUpdate_size[i];
            cInfo_degree[i] += cUpdate_degree[i];
        }

        //Do pointer swaps to reuse memory:
        comm_type* tmp;
        tmp = pastCommAss;
        pastCommAss = currCommAss; //Previous holds the current
        currCommAss = targetCommAss; //Current holds the chosen assignment
        targetCommAss = tmp;      //Reuse the vector

    }//End of while(true)
    *totTime = total; //Return back the total time for clustering
    *numItr  = numItrs;

#ifdef PRINT_DETAILED_STATS_
    printf("========================================================================================================\n");
    printf("Total time for %d iterations is: %lf\n",numItrs, total);
    printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
    printf("========================================================================================================\n");
    printf("Total time for %d iterations is: %lf\n",numItrs, total);
    printf("========================================================================================================\n");
#endif

    //Store back the community assignments in the input variable:
    //Note: No matter when the while loop exits, we are interested in the previous assignment
#pragma omp parallel for
    for (comm_type i=0; i<NV; i++) {
//        C[i] = pastCommAss[i];
        C[i] = currCommAss[i];
    }

    for (int i = 0; i < nnz; ++i) {
        vtxInd[i].head = head[i];
        vtxInd[i].tail = tail[i];
        vtxInd[i].weight = weights[i];
    }
    //Cleanup
    free(pastCommAss);
    free(currCommAss);
    free(targetCommAss);
    free(vDegree);
    free(head);
    free(tail);
    free(weights);
    free(cid);
    free(Counter);
    free(cInfo_size);
    free(cInfo_degree);
    free(cUpdate_size);
    free(cUpdate_degree);
    free(clusterWeightInternal);
    return prevMod;
}
f_weight vectorizedLouvianMethod(graph *G, comm_type *C, int nThreads, f_weight Lower,
                               f_weight thresh, double *totTime, int *numItr, bool *change) {

    cout<< "Vectorized version called" <<endl;
//    size_t alignment = sysconf(_SC_PAGESIZE);
    size_t alignment = 512;
#ifdef PRINT_DETAILED_STATS_
    printf("Within parallelLouvianMethod()\n");
#endif
    if (nThreads < 1)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(nThreads);
    int nT;
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
#ifdef PRINT_DETAILED_STATS_
    printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
#endif


    double time1, time2, time3, time4; //For timing purposes
    double total = 0, totItr = 0;

    comm_type    NV        = G->numVertices;
    comm_type    NS        = G->sVertices;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;

    comm_type nnz = 0;
    for (comm_type i=0; i<NV; i++) {
        comm_type adj1 = vtxPtr[i];        //Begin
        comm_type adj2 = vtxPtr[i + 1];    //End
        for (comm_type j = adj1; j < adj2; j++) {
            nnz++;
        }
    }
    cout << "NNZ: " << nnz << " NE: " << NE << " twice NE: " << (2*NE) << endl;
    /* Variables for computing modularity */
    comm_type totalEdgeWeightTwice;
    f_weight constantForSecondTerm;
    f_weight prevMod=-1;
    f_weight currMod=-1;
    //double thresMod = 0.000001;
    f_weight thresMod = thresh; //Input parameter
    int numItrs = 0;

    /********************** Initialization **************************/
    time1 = omp_get_wtime();
    //Store the degree of all vertices
    f_weight * vDegree; // = (f_weight *) malloc (NV * sizeof(f_weight));
    posix_memalign((void **) &vDegree, alignment, NV * sizeof(f_weight));
    assert(vDegree != 0);

    /// pointer to track existing community
    cout << "Error check pint " << endl;
    comm_type** track_cid = (comm_type **) malloc (nThreads * sizeof(comm_type *));
//    posix_memalign((void **) track_cid, alignment, nThreads * sizeof(comm_type *));
    assert(track_cid != 0);
    for (int k = 0; k < nThreads; ++k) {
        posix_memalign((void **) &track_cid[k], alignment, NV * sizeof(comm_type));
        assert(track_cid[k] != 0);
    }
    cout << "Initialization okay " << endl;

    //Community info. (ai and size)
    /*Comm *cInfo; // = (Comm *) malloc (NV * sizeof(Comm));
    posix_memalign((void **) &cInfo, alignment, NV * sizeof(Comm));
    assert(cInfo != 0);*/
    /// replace cInfo by the following cInfo_size and cInfo_degree
    comm_type* cInfo_size;
    posix_memalign((void **) &cInfo_size, alignment, NV * sizeof(comm_type));
    assert(cInfo_size != 0);
    f_weight* cInfo_degree;
    posix_memalign((void **) &cInfo_degree, alignment, NV * sizeof(f_weight));
    assert(cInfo_degree != 0);

    //use for updating Community
    /*Comm *cUpdate; // = (Comm*)malloc(NV*sizeof(Comm));
    posix_memalign((void **) &cUpdate, alignment, NV * sizeof(Comm));
    assert(cUpdate != 0);*/
    /// replace cUpdate by the following cUpdate_size and cUpdate_degree
    comm_type* cUpdate_size;
    posix_memalign((void **) &cUpdate_size, alignment, NV * sizeof(comm_type));
    assert(cUpdate_size != 0);
    f_weight* cUpdate_degree;
    posix_memalign((void **) &cUpdate_degree, alignment, NV * sizeof(f_weight));
    assert(cUpdate_degree != 0);

    //use for Modularity calculation (eii)
    f_weight* clusterWeightInternal; // = (f_weight*) malloc (NV*sizeof(f_weight));
    posix_memalign((void **) &clusterWeightInternal, alignment, NV * sizeof(f_weight));
    assert(clusterWeightInternal != 0);

    sumVertexDegreeVec_sfp(vtxInd, vtxPtr, vDegree, NV , cInfo_size, cInfo_degree);	// Sum up the vertex degree

    /*** Compute the total edge weight (2m) and 1/2m ***/
    constantForSecondTerm = calConstantForSecondTerm_sfp(vDegree, NV); // 1 over sum of the degree

    //Community assignments:
    //Store previous iteration's community assignment
    comm_type* pastCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &pastCommAss, alignment, NV * sizeof(comm_type));
    assert(pastCommAss != 0);
    //Store current community assignment
    comm_type* currCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &currCommAss, alignment, NV * sizeof(comm_type));
    assert(currCommAss != 0);
    //Store the target of community assignment
    comm_type* targetCommAss; // = (comm_type *) malloc (NV * sizeof(comm_type));
    posix_memalign((void **) &targetCommAss, alignment, NV * sizeof(comm_type));
    assert(targetCommAss != 0);

    comm_type* head;
    posix_memalign((void **) &head, alignment, (nnz * sizeof(comm_type)));
    assert(head != 0);
    comm_type* tail;
    posix_memalign((void **) &tail, alignment, (nnz * sizeof(comm_type)));
    assert(tail != 0);
    f_weight* weights;
    posix_memalign((void **) &weights, alignment, (nnz * sizeof(f_weight)));
    assert(weights != 0);
    for (int i = 0; i < nnz; ++i) {
        head[i] = vtxInd[i].head;
        tail[i] = vtxInd[i].tail;
        weights[i] = vtxInd[i].weight;
    }

    //Vectors used in place of maps: Total size = |V|+2*|E| -- The |V| part takes care of self loop
    /*mapElement* clusterLocalMap; // = (mapElement *) malloc ((NV + 2*NE) * sizeof(mapElement));
    posix_memalign((void **) &clusterLocalMap, alignment, ((NV + 2*NE) * sizeof(mapElement)));
    assert(clusterLocalMap != 0);*/
    /// Replace clusterLocalMap by the following cid and Counter
    comm_type* cid;
    posix_memalign((void **) &cid, alignment, ((NV + 2*NE) * sizeof(comm_type)));
    assert(cid != 0);
    f_weight* Counter;
    posix_memalign((void **) &Counter, alignment, ((NV + 2*NE) * sizeof(f_weight)));
    assert(Counter != 0);
    //double* Counter             = (double *)     malloc ((NV + 2*NE) * sizeof(double));     assert(Counter != 0);

    //Initialize each vertex to its own cluster
    //initCommAssOpt(pastCommAss, currCommAss, NV, clusterLocalMapX, vtxPtr, vtxInd, cInfo, constantForSecondTerm, vDegree);

    //Initialize each vertex to its own cluster
//    initCommAss_SFP(pastCommAss, currCommAss, NV);
    initCommAssOptVec_SFP(pastCommAss, currCommAss, NV, cid, Counter, vtxPtr, head, tail, weights, cInfo_size, cInfo_degree, constantForSecondTerm, vDegree);
    time2 = omp_get_wtime();


    comm_type* test_cid;
    posix_memalign((void **) &test_cid, alignment, ((NV + 2*NE) * sizeof(comm_type)));
    assert(test_cid != 0);
    f_weight* test_Counter;
    posix_memalign((void **) &test_Counter, alignment, ((NV + 2*NE) * sizeof(f_weight)));
    assert(test_Counter != 0);

    for (int i = 0; i < (NV + 2*NE); ++i) {
        test_cid[i] = cid[i];
        test_Counter[i] = Counter[i];
    }
    printf("Time to initialize: %3.3lf\n", time2-time1);

#ifdef PRINT_DETAILED_STATS_
    printf("========================================================================================================\n");
    printf("Itr      E_xx            A_x2           Curr-Mod         Time-1(s)       Time-2(s)        T/Itr(s)\n");
    printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
    printf("=====================================================\n");
    printf("Itr      Curr-Mod         T/Itr(s)      T-Cumulative\n");
    printf("=====================================================\n");
#endif
    //Start maximizing modularity
    while(true) {
        numItrs++;
        time1 = omp_get_wtime();
        /* Re-initialize datastructures */
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            clusterWeightInternal[i] = 0;
            cUpdate_degree[i] =0;
            cUpdate_size[i] =0;
        }


        bool moved = false;
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            comm_type tid = omp_get_thread_num();
//            if(numItrs == 1)
//                cout << "Tid: " << tid << endl;
            comm_type adj1 = vtxPtr[i];
            comm_type adj2 = vtxPtr[i+1];
            f_weight selfLoop = 0;
            f_weight selfLoop2 = 0;
            //Build a datastructure to hold the cluster structure of its neighbors
//            map<comm_type, comm_type> clusterLocalMap; //Map each neighbor's cluster to a local number
//            map<comm_type, comm_type>::iterator storedAlready;
//            vector<f_weight> Counter; //Number of edges in each unique cluster
            comm_type numUniqueClusters = 0;
            //Add v's current cluster:
            if(adj1 != adj2){
//                clusterLocalMap[currCommAss[i]] = 0;
//                Counter.push_back(0); //Initialize the counter to ZERO (no edges incident yet)
                comm_type sPosition = vtxPtr[i]+i; //Starting position of local map for i
                Counter[sPosition] = 0;          //Initialize the counter to ZERO (no edges incident yet)
                cid[sPosition] = currCommAss[i]; //Initialize with current community
                numUniqueClusters++; //Added the first entry

                //Find unique cluster ids and #of edges incident (eicj) to them
//                selfLoop = buildLocalMapCounter_sfp(adj1, adj2, clusterLocalMap, Counter, vtxInd, currCommAss, i);
//                selfLoop = buildLocalMapCounterVec_SFP(i, cid, Counter, vtxPtr, head, tail, weights, currCommAss, numUniqueClusters);

//                comm_type* tmp_ptr = &track_cid[tid][0];
                selfLoop2 = buildLocalMapCounterNoMap_SFP(i, test_cid, test_Counter, vtxPtr, head, tail, weights, currCommAss, numUniqueClusters);
                selfLoop = buildLocalMapCounterVec2nd_SFP(i, cid, Counter, vtxPtr, head, tail, weights, currCommAss, numUniqueClusters, &track_cid[tid][0]);
                if(selfLoop2 != selfLoop){
                    cout << "Problem: " << " parallel self-loop: " << selfLoop2 << " vectorized self-loop: " << selfLoop << endl;
//                    break;
                }
                // Update delta Q calculation
//                clusterWeightInternal[i] += Counter[0]; //(e_ix)
                clusterWeightInternal[i] += Counter[sPosition]; //(e_ix)
                //Calculate the max
//                targetCommAss[i] = max_sfp(clusterLocalMap, Counter, selfLoop, cInfo, vDegree[i], currCommAss[i], constantForSecondTerm);
//                targetCommAss[i] = maxNoMap_SFP(i, cid, Counter, vtxPtr, selfLoop, cInfo_size, cInfo_degree, vDegree[i], currCommAss[i], constantForSecondTerm, numUniqueClusters);
                targetCommAss[i] = maxNoMapVec_SFP(i, cid, Counter, vtxPtr, selfLoop, cInfo_size, cInfo_degree, vDegree[i], currCommAss[i], constantForSecondTerm, numUniqueClusters);
                //assert((targetCommAss[i] >= 0)&&(targetCommAss[i] < NV));
            } else {
                targetCommAss[i] = -1;
            }

//            cout<< "Vertex " << i  << " of " << NV << " processed" << endl;
            //Update
            if(targetCommAss[i] != currCommAss[i]  && targetCommAss[i] != -1) {
                moved = true;
#pragma omp atomic update
                cUpdate_degree[targetCommAss[i]] += vDegree[i];
#pragma omp atomic update
                cUpdate_size[targetCommAss[i]] += 1;
#pragma omp atomic update
                cUpdate_degree[currCommAss[i]] -= vDegree[i];
#pragma omp atomic update
                cUpdate_size[currCommAss[i]] -=1;
                /*
                 __sync_fetch_and_add(&cUpdate[targetCommAss[i]].size, 1);
                 __sync_fetch_and_sub(&cUpdate[currCommAss[i]].degree, vDegree[i]);
                 __sync_fetch_and_sub(&cUpdate[currCommAss[i]].size, 1);*/
            }//End of If()
        }//End of for(i)
        time2 = omp_get_wtime();

        time3 = omp_get_wtime();
        f_weight e_xx = 0;
        f_weight a2_x = 0;

#pragma omp parallel for \
reduction(+:e_xx) reduction(+:a2_x)
        for (comm_type i=0; i<NV; i++) {
            e_xx += clusterWeightInternal[i];
            a2_x += (cInfo_degree[i])*(cInfo_degree[i]);
        }
        time4 = omp_get_wtime();

        currMod = (e_xx*(f_weight)constantForSecondTerm) - (a2_x*(f_weight)constantForSecondTerm*(f_weight)constantForSecondTerm);
        totItr = (time2-time1) + (time4-time3);
        total += totItr;
#ifdef PRINT_DETAILED_STATS_
        printf("%d \t %g \t %g \t %lf \t %3.3lf \t %3.3lf  \t %3.3lf\n",numItrs, e_xx, a2_x, currMod, (time2-time1), (time4-time3), totItr );
#endif
#ifdef PRINT_TERSE_STATS_
        printf("%d \t %lf \t %3.3lf  \t %3.3lf\n",numItrs, currMod, totItr, total);
#endif

        if(moved)
            *change = true;

        //Break if modularity gain is not sufficient
        if(!moved || numItrs >= 25/*(currMod - prevMod) < thresMod*/) {
            break;
        }

        //Else update information for the next iteration
        prevMod = currMod;
        if(prevMod < Lower)
            prevMod = Lower;
#pragma omp parallel for
        for (comm_type i=0; i<NV; i++) {
            cInfo_size[i] += cUpdate_size[i];
            cInfo_degree[i] += cUpdate_degree[i];
        }

        //Do pointer swaps to reuse memory:
        comm_type* tmp;
        tmp = pastCommAss;
        pastCommAss = currCommAss; //Previous holds the current
        currCommAss = targetCommAss; //Current holds the chosen assignment
        targetCommAss = tmp;      //Reuse the vector

    }//End of while(true)
    *totTime = total; //Return back the total time for clustering
    *numItr  = numItrs;

#ifdef PRINT_DETAILED_STATS_
    printf("========================================================================================================\n");
    printf("Total time for %d iterations is: %lf\n",numItrs, total);
    printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
    printf("========================================================================================================\n");
    printf("Total time for %d iterations is: %lf\n",numItrs, total);
    printf("========================================================================================================\n");
#endif

    //Store back the community assignments in the input variable:
    //Note: No matter when the while loop exits, we are interested in the previous assignment
#pragma omp parallel for
    for (comm_type i=0; i<NV; i++) {
//        C[i] = pastCommAss[i];
        C[i] = currCommAss[i];
    }

    for (comm_type i = 0; i < nnz; ++i) {
        vtxInd[i].head = head[i];
        vtxInd[i].tail = tail[i];
        vtxInd[i].weight = weights[i];
    }
    //Cleanup
    free(pastCommAss);
    free(currCommAss);
    free(targetCommAss);
    free(vDegree);
    free(head);
    free(tail);
    free(weights);
    free(cid);
    free(Counter);
    free(cInfo_size);
    free(cInfo_degree);
    free(cUpdate_size);
    free(cUpdate_degree);
    free(clusterWeightInternal);

    return prevMod;
}
