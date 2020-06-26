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
#include "color_comm.h"
using namespace std;

double algoLouvainWithDistOneColoringNoMap(graph* G, comm_type *C, int nThreads, int* color,
			int numColor, double Lower, double thresh, double *totTime, int *numItr) {
#ifdef PRINT_DETAILED_STATS_  
	printf("Within algoLouvainWithDistOneColoring()\n");
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
	/* Indexs are vertex */
	comm_type* pastCommAss;	//Store previous iteration's community assignment
	comm_type* currCommAss;	//Store current community assignment
	//comm_type* targetCommAss;	//Store the target of community assignment
  double* vDegree;	//Store each vertex's degree
	double* clusterWeightInternal;//use for Modularity calculation (eii)
	
	/* Indexs are community */
	Comm* cInfo;	 //Community info. (ai and size)
	Comm* cUpdate; //use for updating Community
	
	/* Book keeping variables */
	comm_type    NV        = G->numVertices;
	comm_type    NS        = G->sVertices;
	comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;
	
	/* Modularity Needed variables */
	comm_type totalEdgeWeightTwice;
	double constantForSecondTerm;
	double prevMod=Lower;
	double currMod=-1;
	double thresMod = thresh;
	int numItrs = 0;
		
	/********************** Initialization **************************/
	time1 = omp_get_wtime();
	vDegree = (double *) malloc (NV * sizeof(double)); assert(vDegree != 0);
	cInfo = (Comm *) malloc (NV * sizeof(Comm)); assert(cInfo != 0);
	cUpdate = (Comm*)malloc(NV*sizeof(Comm)); assert(cUpdate != 0);

	sumVertexDegree(vtxInd, vtxPtr, vDegree, NV , cInfo);	// Sum up the vertex degree
	/*** Compute the total edge weight (2m) and 1/2m ***/
	constantForSecondTerm = calConstantForSecondTerm(vDegree, NV);	// 1 over sum of the degree
	
	pastCommAss = (comm_type *) malloc (NV * sizeof(comm_type)); assert(pastCommAss != 0);
	//Community provided as input:
	currCommAss = C; assert(currCommAss != 0);

    //Vectors used in place of maps: Total size = |V|+2*|E| -- The |V| part takes care of self loop
    mapElement* clusterLocalMap = (mapElement *) malloc ((NV + 2*NE) * sizeof(mapElement)); assert(clusterLocalMap != 0);
    
    /*** Assign each vertex to its own Community ***/
	initCommAss( pastCommAss, currCommAss, NV);

	clusterWeightInternal = (double*) malloc (NV*sizeof(double)); assert(clusterWeightInternal != 0);
	
	/*** Create a CSR-like datastructure for vertex-colors ***/
	comm_type * colorPtr = (comm_type *) malloc ((numColor+1) * sizeof(comm_type));
	comm_type * colorIndex = (comm_type *) malloc (NV * sizeof(comm_type));
	comm_type * colorAdded = (comm_type *)malloc (numColor*sizeof(comm_type));
	assert(colorPtr != 0);
        assert(colorIndex != 0);
	assert(colorAdded != 0);
	// Initialization
#pragma omp parallel for
	for(comm_type i = 0; i < numColor; i++) {
		colorPtr[i] = 0;
		colorAdded[i] = 0;
	}
	colorPtr[numColor] = 0;
	// Count the size of each color
#pragma omp parallel for
	for(comm_type i = 0; i < NV; i++) {
		__sync_fetch_and_add(&colorPtr[(comm_type)color[i]+1],1);
	}
	//Prefix sum:
	for(comm_type i=0; i<numColor; i++) {
		colorPtr[i+1] += colorPtr[i];
	}	
	//Group vertices with the same color in particular order
#pragma omp parallel for
	for (comm_type i=0; i<NV; i++) {
		comm_type tc = (comm_type)color[i];
		comm_type Where = colorPtr[tc] + __sync_fetch_and_add(&(colorAdded[tc]), 1);
		colorIndex[Where] = i;
	}
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
	while(true) {
		numItrs++;
		
		time1 = omp_get_wtime();
		for( comm_type ci = 0; ci < numColor; ci++) // Begin of color loop
		{
#pragma omp parallel for
			for (comm_type i=0; i<NV; i++) {
				clusterWeightInternal[i] = 0; //Initialize to zero
				cUpdate[i].degree =0;
				cUpdate[i].size =0;
			}
			comm_type coloradj1 = colorPtr[ci];
			comm_type coloradj2 = colorPtr[ci+1];
			
#pragma omp parallel for  
			for (comm_type K = coloradj1; K<coloradj2; K++) {
				comm_type i = colorIndex[K];
				comm_type localTarget = -1;
				comm_type adj1 = vtxPtr[i];
				comm_type adj2 = vtxPtr[i+1];
				double selfLoop = 0;
				//Build a datastructure to hold the cluster structure of its neighbors:      	
				//map<comm_type, comm_type> clusterLocalMap; //Map each neighbor's cluster to a local number
				//map<comm_type, comm_type>::iterator storedAlready;
				//vector<double> Counter; //Number of edges to each unique cluster
				comm_type numUniqueClusters = 0;
				if(adj1 != adj2) {
					//Add the current cluster of i to the local map
                    comm_type sPosition = vtxPtr[i]+i; //Starting position of local map for i
                    clusterLocalMap[sPosition].Counter = 0;          //Initialize the counter to ZERO (no edges incident yet)
                    clusterLocalMap[sPosition].cid = currCommAss[i]; //Initialize with current community
                    numUniqueClusters++; //Added the first entry
                    
					//Find unique cluster ids and #of edges incident (eicj) to them
					selfLoop = buildLocalMapCounterNoMap(i, clusterLocalMap, vtxPtr, vtxInd, currCommAss, numUniqueClusters);
					//Calculate the max
					localTarget = maxNoMap(i, clusterLocalMap, vtxPtr, selfLoop, cInfo, vDegree[i], currCommAss[i], constantForSecondTerm, numUniqueClusters);
				} else {
					localTarget = -1;
				}					
				//Update prepare
				if(localTarget != currCommAss[i] && localTarget != -1) {
          #pragma omp atomic update
          cUpdate[localTarget].degree += vDegree[i];
          #pragma omp atomic update
          cUpdate[localTarget].size += 1;
          #pragma omp atomic update
          cUpdate[currCommAss[i]].degree -= vDegree[i];
          #pragma omp atomic update
          cUpdate[currCommAss[i]].size -=1;
      /*
					__sync_fetch_and_add(&cUpdate[localTarget].degree, vDegree[i]);
	         			__sync_fetch_and_add(&cUpdate[localTarget].size, 1);
					__sync_fetch_and_sub(&cUpdate[currCommAss[i]].degree, vDegree[i]);
					__sync_fetch_and_sub(&cUpdate[currCommAss[i]].size, 1);*/
				}//End of If()
				currCommAss[i] = localTarget;      
				//clusterLocalMap.clear();
			}//End of for(i)
			
			// UPDATE
#pragma omp parallel for  
			for (comm_type i=0; i<NV; i++) {
				cInfo[i].size += cUpdate[i].size;
				cInfo[i].degree += cUpdate[i].degree;
			}
		}//End of Color loop						
		time2 = omp_get_wtime();
		
		time3 = omp_get_wtime();    
		double e_xx = 0;
		double a2_x = 0;	
		
		// CALCULATE MOD
#pragma omp parallel for  //Parallelize on each vertex
		for (comm_type i =0; i<NV;i++){
			clusterWeightInternal[i] = 0;
		}
#pragma omp parallel for  //Parallelize on each vertex
		for (comm_type i=0; i<NV; i++) {
			comm_type adj1 = vtxPtr[i];
			comm_type adj2 = vtxPtr[i+1];
			for(comm_type j=adj1; j<adj2; j++) {
				if(currCommAss[vtxInd[j].tail] == currCommAss[i]){
					clusterWeightInternal[i] += vtxInd[j].weight;
				}
			}
		}		
		
#pragma omp parallel for \
reduction(+:e_xx) reduction(+:a2_x)
		for (comm_type i=0; i<NV; i++) {
			e_xx += clusterWeightInternal[i];
			a2_x += (cInfo[i].degree)*(cInfo[i].degree);
		}
		time4 = omp_get_wtime();
		
		currMod = e_xx*(double)constantForSecondTerm  - a2_x*(double)constantForSecondTerm*(double)constantForSecondTerm;
		
		totItr = (time2-time1) + (time4-time3);
		total += totItr;

#ifdef PRINT_DETAILED_STATS_  
		printf("%d \t %g \t %g \t %lf \t %3.3lf \t %3.3lf  \t %3.3lf\n",numItrs, e_xx, a2_x, currMod, (time2-time1), (time4-time3), totItr );    
#endif
#ifdef PRINT_TERSE_STATS_
   		printf("%d \t %lf \t %3.3lf  \t %3.3lf\n",numItrs, currMod, totItr, total);
#endif
		if((currMod - prevMod) < thresMod) {	
			break;
		}
		
		prevMod = currMod;
	}//End of while(true)
	*totTime = total; //Return back the total time
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
	//Cleanup:
        free(vDegree); free(cInfo); free(cUpdate); free(clusterWeightInternal);
        free(colorPtr); free(colorIndex); free(colorAdded);
	free(pastCommAss);
    free(clusterLocalMap);
	
	return prevMod;
	
}//End of algoLouvainWithDistOneColoring()
