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
#include "coloring.h"


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  DISTANCE ONE COLORING      ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//Return the number of colors used (zero is a valid color)
int algoDistanceOneVertexColoringOpt(graph *G, int *vtxColor, int nThreads, double *totTime)
{
#ifdef PRINT_DETAILED_STATS_
  printf("Within algoDistanceOneVertexColoringOpt()\n");
#endif

  if (nThreads < 1)
		omp_set_num_threads(1); //default to one thread
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
	
  double time1=0, time2=0, totalTime=0;
  //Get the iterators for the graph:
  comm_type NVer    = G->numVertices;
  comm_type NEdge   = G->numEdges;
    comm_type *verPtr = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)

#ifdef PRINT_DETAILED_STATS_
  printf("Vertices: %ld  Edges: %ld\n", NVer, NEdge);
#endif

  //Build a vector of random numbers
  double *randValues = (double*) malloc (NVer * sizeof(double));
  assert(randValues != 0);
  generateRandomNumbers(randValues, NVer);

  comm_type *Q    = (comm_type *) malloc (NVer * sizeof(comm_type)); assert(Q != 0);
  comm_type *Qtmp = (comm_type *) malloc (NVer * sizeof(comm_type)); assert(Qtmp != 0);
  comm_type *Qswap;
  if( (Q == NULL) || (Qtmp == NULL) ) {
    printf("Not enough memory to allocate for the two queues \n");
    exit(1);
  }
  comm_type QTail=0;    //Tail of the queue
  comm_type QtmpTail=0; //Tail of the queue (implicitly will represent the size)
  comm_type realMaxDegree = 0;
	
	#pragma omp parallel for
  for (comm_type i=0; i<NVer; i++) {
      Q[i]= i;     //Natural order
      Qtmp[i]= -1; //Empty queue
  }
  QTail = NVer;	//Queue all vertices


	// Cal real Maximum degree, 2x for maxDegree to be safe
	#pragma omp parallel for reduction(max: realMaxDegree)
	for (comm_type i = 0; i < NVer; i++) {
		comm_type adj1, adj2, de;
		adj1 = verPtr[i];
		adj2 = verPtr[i+1];
		de = adj2-adj1;
		if ( de > realMaxDegree)
			realMaxDegree = de;
	}
	//realMaxDegree *= 1.5;

	ColorVector freq(MaxDegree,0);
  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// START THE WHILE LOOP ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  comm_type nConflicts = 0; //Number of conflicts
  int nLoops = 0;     //Number of rounds of conflict resolution

#ifdef PRINT_DETAILED_STATS_ 
  printf("Results from parallel coloring:\n");
  printf("***********************************************\n");
#endif
  do{
    ///////////////////////////////////////// PART 1 ////////////////////////////////////////
    //Color the vertices in parallel - do not worry about conflicts
#ifdef PRINT_DETAILED_STATS_
    printf("** Iteration : %d \n", nLoops);
#endif

    time1 = omp_get_wtime();
		#pragma omp parallel for
    for (comm_type Qi=0; Qi<QTail; Qi++) {
      comm_type v = Q[Qi]; //Q.pop_front();
			int maxColor = 0;
			BitVector mark(MaxDegree, false);
			maxColor = distanceOneMarkArray(mark,G,v,vtxColor);
				
			int myColor;
			for (myColor=0; myColor<=maxColor; myColor++) {
				if ( mark[myColor] == false )
					break;
			}     
			vtxColor[v] = myColor; //Color the vertex
		} //End of outer for loop: for each vertex
		
		time1  = omp_get_wtime() - time1;
		totalTime += time1;

#ifdef PRINT_DETAILED_STATS_
    printf("Time taken for Coloring:  %lf sec.\n", time1);
#endif
    ///////////////////////////////////////// PART 2 ////////////////////////////////////////
    //Detect Conflicts:
    //printf("Phase 2: Detect Conflicts, add to queue\n");    
    //Add the conflicting vertices into a Q:
    //Conflicts are resolved by changing the color of only one of the 
    //two conflicting vertices, based on their random values 
    time2 = omp_get_wtime();
		
		#pragma omp parallel for
		for (comm_type Qi=0; Qi<QTail; Qi++) {
			comm_type v = Q[Qi]; //Q.pop_front();
			distanceOneConfResolution(G, v, vtxColor, randValues, &QtmpTail, Qtmp, freq, 0);
		} //End of outer for loop: for each vertex
  
		time2  = omp_get_wtime() - time2;
		totalTime += time2;    
		nConflicts += QtmpTail;
		nLoops++;

#ifdef PRINT_DETAILED_STATS_
    printf("Num conflicts      : %ld \n", QtmpTail);
    printf("Time for detection : %lf sec\n", time2);
#endif

    //Swap the two queues:
    Qswap = Q;
    Q = Qtmp; //Q now points to the second vector
    Qtmp = Qswap;
    QTail = QtmpTail; //Number of elements
    QtmpTail = 0; //Symbolic emptying of the second queue    
  } while (QTail > 0);
  //Check the number of colors used
  int nColors = -1;
  for (comm_type v=0; v < NVer; v++ )
    if (vtxColor[v] > nColors) nColors = vtxColor[v];
#ifdef PRINT_DETAILED_STATS_
  printf("***********************************************\n");
  printf("Total number of colors used: %d \n", nColors);    
  printf("Number of conflicts overall: %ld \n", nConflicts);  
  printf("Number of rounds           : %d \n", nLoops);      
  printf("Total Time                 : %lf sec\n", totalTime);
  printf("***********************************************\n");
#endif  
  *totTime = totalTime;
  //////////////////////////// /////////////////////////////////////////////////////////////
  ///////////////////////////////// VERIFY THE COLORS /////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  //Verify Results and Cleanup
  int myConflicts = 0;
	#pragma omp parallel for
  for (comm_type v=0; v < NVer; v++ ) {
    comm_type adj1 = verPtr[v];
    comm_type adj2 = verPtr[v+1];
    //Browse the adjacency set of vertex v
    for(comm_type k = adj1; k < adj2; k++ ) {
      if ( v == verInd[k].tail ) //Self-loops
	continue;
      if ( vtxColor[v] == vtxColor[verInd[k].tail] ) {
	__sync_fetch_and_add(&myConflicts, 1); //increment the counter
      }
    }//End of inner for loop: w in adj(v)
  }//End of outer for loop: for each vertex
  myConflicts = myConflicts / 2; //Have counted each conflict twice
	
  if (myConflicts > 0)
    printf("Check - WARNING: Number of conflicts detected after resolution: %d \n\n", myConflicts);
  else
    printf("Check - SUCCESS: No conflicts exist\n\n");
  //Clean Up:
  free(Q);
  free(Qtmp);
  free(randValues);
  
  return nColors; //Return the number of colors used
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  DISTANCE ONE COLORING      ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//Return the number of colors used (zero is a valid color)
int algoDistanceOneVertexColoring(graph *G, int *vtxColor, int nThreads, double *totTime)
{
	printf("Within algoDistanceOneVertexColoring()\n");
	if (nThreads < 1)
		omp_set_num_threads(1); //default to one thread
	else
		omp_set_num_threads(nThreads);
	int nT;
#pragma omp parallel
	{
		nT = omp_get_num_threads();
	}
	printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
	
	
	double time1=0, time2=0, totalTime=0;
  //Get the iterators for the graph:
  comm_type NVer    = G->numVertices;
  comm_type NS      = G->sVertices;
  comm_type NT      = NVer - NS;
  comm_type NEdge           = G->numEdges;
    comm_type *verPtr         = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
  edge *verInd         = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
  printf("Vertices: %ld  Edges: %ld\n", NVer, NEdge);

  //const int MaxDegree = 4096; //Increase if number of colors is larger    

  //Build a vector of random numbers
  double *randValues = (double*) malloc (NVer * sizeof(double));
  if( randValues == NULL ) {
    printf("Not enough memory to allocate for random numbers \n");
    exit(1);
  }
  generateRandomNumbers(randValues, NVer);

  //The Queue Data Structure for the storing the vertices 
  //   the need to be colored/recolored
  //Have two queues - read from one, write into another
  //   at the end, swap the two.
  comm_type *Q    = (comm_type *) malloc (NVer * sizeof(comm_type));
  comm_type *Qtmp = (comm_type *) malloc (NVer * sizeof(comm_type));
  comm_type *Qswap;
  if( (Q == NULL) || (Qtmp == NULL) ) {
    printf("Not enough memory to allocate for the two queues \n");
    exit(1);
  }
  comm_type QTail=0;    //Tail of the queue
  comm_type QtmpTail=0; //Tail of the queue (implicitly will represent the size)
  
#pragma omp parallel for
  for (comm_type i=0; i<NVer; i++) {
      Q[i]= i;     //Natural order
      Qtmp[i]= -1; //Empty queue
  }
  QTail = NVer;	//Queue all vertices
  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// START THE WHILE LOOP ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  comm_type nConflicts = 0; //Number of conflicts
  int nLoops = 0;     //Number of rounds of conflict resolution
  int *Mark = (int *) malloc ( MaxDegree * NVer * sizeof(int) );
  if( Mark == NULL ) {
    printf("Not enough memory to allocate for Mark \n");
    exit(1);
  }
#pragma omp parallel for
  for (comm_type i=0; i<MaxDegree*NVer; i++)
     Mark[i]= -1;

  printf("Results from parallel coloring:\n");
  printf("***********************************************\n");
  do {
    ///////////////////////////////////////// PART 1 ////////////////////////////////////////
    //Color the vertices in parallel - do not worry about conflicts
    printf("** Iteration : %d \n", nLoops);
    time1 = omp_get_wtime();
#pragma omp parallel for
    for (comm_type Qi=0; Qi<QTail; Qi++) {
      comm_type v = Q[Qi]; //Q.pop_front();
      comm_type StartIndex = v*MaxDegree; //Location in Mark
      if (nLoops > 0) //Skip the first time around
	for (comm_type i=StartIndex; i<(StartIndex+MaxDegree); i++)
	  Mark[i]= -1;
      comm_type adj1 = verPtr[v];
      comm_type adj2 = verPtr[v+1];
      int maxColor = -1;
      int adjColor = -1;
      //Browse the adjacency set of vertex v
      for(comm_type k = adj1; k < adj2; k++ ) {
	//if ( v == verInd[k] ) //Skip self-loops
	//continue;
	adjColor =  vtxColor[verInd[k].tail];
	if ( adjColor >= 0 ) {
	  Mark[StartIndex+adjColor] = v;
	  //Find the largest color in the neighborhood
	  if ( adjColor > maxColor )
	    maxColor = adjColor;
	}
      } //End of for loop to traverse adjacency of v
      int myColor;
      for (myColor=0; myColor<=maxColor; myColor++) {
	if ( Mark[StartIndex+myColor] != v )
	  break;
      }
      if (myColor == maxColor)
	myColor++; /* no available color with # less than cmax */
      vtxColor[v] = myColor; //Color the vertex
    } //End of outer for loop: for each vertex
    time1  = omp_get_wtime() - time1;
    totalTime += time1;
    printf("Time taken for Coloring:  %lf sec.\n", time1);

    ///////////////////////////////////////// PART 2 ////////////////////////////////////////
    //Detect Conflicts:
    //printf("Phase 2: Detect Conflicts, add to queue\n");    
    //Add the conflicting vertices into a Q:
    //Conflicts are resolved by changing the color of only one of the 
    //two conflicting vertices, based on their random values 
    time2 = omp_get_wtime();
#pragma omp parallel for
    for (comm_type Qi=0; Qi<QTail; Qi++) {
      comm_type v = Q[Qi]; //Q.pop_front();
      comm_type adj1 = verPtr[v];
      comm_type adj2 = verPtr[v+1];
      //Browse the adjacency set of vertex v
      for(comm_type k = adj1; k < adj2; k++ ) {
	//if ( v == verInd[k] ) //Self-loops
	//continue;
	if ( vtxColor[v] == vtxColor[verInd[k].tail] ) {
	  //Q.push_back(v or w)
	  if ( (randValues[v] < randValues[verInd[k].tail]) || 
	       ((randValues[v] == randValues[verInd[k].tail])&&(v < verInd[k].tail)) ) {
	    comm_type whereInQ = __sync_fetch_and_add(&QtmpTail, 1);
	    Qtmp[whereInQ] = v;//Add to the queue
	    vtxColor[v] = -1;  //Will prevent v from being in conflict in another pairing
	    break;
	  }
	} //End of if( vtxColor[v] == vtxColor[verInd[k]] )
      } //End of inner for loop: w in adj(v)
    } //End of outer for loop: for each vertex
    time2  = omp_get_wtime() - time2;
    totalTime += time2;    
    nConflicts += QtmpTail;
    nLoops++;
    printf("Conflicts          : %ld \n", QtmpTail);
    printf("Time for detection : %lf sec\n", time2);
    //Swap the two queues:
    Qswap = Q;
    Q = Qtmp; //Q now points to the second vector
    Qtmp = Qswap;
    QTail = QtmpTail; //Number of elements
    QtmpTail = 0; //Symbolic emptying of the second queue    
  } while (QTail > 0);
  //Check the number of colors used
  int nColors = -1;
  for (comm_type v=0; v < NVer; v++ )
    if (vtxColor[v] > nColors) nColors = vtxColor[v];
  printf("***********************************************\n");
  printf("Total number of colors used: %d \n", nColors);    
  printf("Number of conflicts overall: %ld \n", nConflicts);  
  printf("Number of rounds           : %d \n", nLoops);      
  printf("Total Time                 : %lf sec\n", totalTime);
  printf("***********************************************\n");
  
  *totTime = totalTime;
  /////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// VERIFY THE COLORS /////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  //Verify Results and Cleanup
  int myConflicts = 0;
#pragma omp parallel for
  for (comm_type v=0; v < NVer; v++ ) {
    comm_type adj1 = verPtr[v];
    comm_type adj2 = verPtr[v+1];
    //Browse the adjacency set of vertex v
    for(comm_type k = adj1; k < adj2; k++ ) {
      if ( v == verInd[k].tail ) //Self-loops
	continue;
      if ( vtxColor[v] == vtxColor[verInd[k].tail] ) {
	//#pragma omp atomic
	//printf("Conflict: color[%ld]=%d AND color[%ld]=%d\n", v, vtxColor[v], verInd[k].tail, vtxColor[ verInd[k].tail]);
	__sync_fetch_and_add(&myConflicts, 1); //increment the counter
      }
    }//End of inner for loop: w in adj(v)
  }//End of outer for loop: for each vertex
  myConflicts = myConflicts / 2; //Have counted each conflict twice
  if (myConflicts > 0)
    printf("Check - WARNING: Number of conflicts detected after resolution: %d \n\n", myConflicts);
  else
    printf("Check - SUCCESS: No conflicts exist\n\n");
  //Clean Up:
  free(Q);
  free(Qtmp);
  free(Mark); 
  free(randValues);
  
  return nColors; //Return the number of colors used
}
