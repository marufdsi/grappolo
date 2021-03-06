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
#include "RngStream.h"

using namespace std;

void generateRandomNumbers(double *RandVec, comm_type size) {
    int nT;
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
#ifdef PRINT_DETAILED_STATS_
    printf("Within generateRandomNumbers() -- Number of threads: %d\n", nT);
#endif
    //Initialize parallel pseudo-random number generator
    unsigned long seed[6] = {1, 2, 3, 4, 5, 6};
    RngStream::SetPackageSeed(seed);
    RngStream RngArray[nT]; //array of RngStream Objects
    
    comm_type block = size / nT;
#ifdef PRINT_DETAILED_STATS_
    cout<<"Each thread will add "<<block<<" edges\n";
#endif
    //Each thread will generate m/nT edges each
    double start = omp_get_wtime();
#pragma omp parallel
    {
        int myRank = omp_get_thread_num();
#pragma omp for schedule(static)
        for (comm_type i=0; i<size; i++) {
            RandVec[i] =  RngArray[myRank].RandU01();
        }
    }//End of parallel region
} //End of generateRandomNumbers()

void displayGraph(graph *G) {
    comm_type    NV        = G->numVertices;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("***********************************");
    printf("|V|= %ld, |E|= %ld \n", NV, NE);
    printf("***********************************");
    for (comm_type i = 0; i < NV; i++) {
        comm_type adj1 = vtxPtr[i];
        comm_type adj2 = vtxPtr[i+1];
        printf("\nVtx: %ld [%ld]: ",i+1,adj2-adj1);
        for(comm_type j=adj1; j<adj2; j++) {
            printf("%ld (%g), ", vtxInd[j].tail+1, vtxInd[j].weight);
        }
    }
    printf("\n***********************************\n");
}

void duplicateGivenGraph(graph *Gin, graph *Gout) {
    comm_type    NV        = Gin->numVertices;
    comm_type    NS        = Gin->sVertices;
    comm_type    NE        = Gin->numEdges;
    comm_type    *vtxPtr   = Gin->edgeListPtrs;
    edge    *vtxInd   = Gin->edgeList;
#ifdef PRINT_DETAILED_STATS_
    printf("|V|= %ld, |E|= %ld \n", NV, NE);
#endif
    double time1 = omp_get_wtime();
    comm_type *edgeListPtr = (comm_type *)  malloc((NV+1) * sizeof(comm_type));
    assert(edgeListPtr != NULL);
#pragma omp parallel for
    for (comm_type i=0; i<=NV; i++) {
        edgeListPtr[i] = vtxPtr[i]; //Prefix Sum
    }
    
    //WARNING: There is a bug in edge counting when self-loops exist
    edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge));
    assert( edgeList != NULL);
#pragma omp parallel for
    for (comm_type i=0; i<NV; i++) {
        for (comm_type j=vtxPtr[i]; j<vtxPtr[i+1]; j++) {
            edgeList[j].head = vtxInd[j].head;
            edgeList[j].tail = vtxInd[j].tail;
            edgeList[j].weight = vtxInd[j].weight;
        }
    }
    
    //The last element of Cumulative will hold the total number of characters
    double time2 = omp_get_wtime();
#ifdef PRINT_DETAILED_STATS_
    printf("Done duplicating the graph:  %9.6lf sec.\n", time2 - time1);
#endif
    Gout->sVertices    = NS;
    Gout->numVertices  = NV;
    Gout->numEdges     = NE;
    Gout->edgeListPtrs = edgeListPtr;
    Gout->edgeList     = edgeList;
} //End of duplicateGivenGraph()

void displayGraphEdgeList(graph *G) {
    comm_type    NV        = G->numVertices;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("***********************************");
    printf("|V|= %ld, |E|= %ld \n", NV, NE);
    for (comm_type i = 0; i < NV; i++) {
        comm_type adj1 = vtxPtr[i];
        comm_type adj2 = vtxPtr[i+1];
        for(comm_type j=adj1; j<adj2; j++) {
            printf("%ld %ld %g\n", i+1, vtxInd[j].tail+1, vtxInd[j].weight);
        }
    }
    printf("\n***********************************\n");
}

void displayGraphEdgeList(graph *G, FILE* out) {
    comm_type    NV        = G->numVertices;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("********PRINT OUTPUT********************");
    fprintf(out,"p sp %ld %ld \n", NV, NE/2);
    for (comm_type i = 0; i < NV; i++) {
        comm_type adj1 = vtxPtr[i];
        comm_type adj2 = vtxPtr[i+1];
        for(comm_type j=adj1; j<adj2; j++) {
            if( i+1 < vtxInd[j].tail+1)
            {
                fprintf(out,"a %ld %ld %g\n", i+1, vtxInd[j].tail+1, vtxInd[j].weight);
            }
        }
    }
}

void writeEdgeListToFile(graph *G, FILE* out) {
    comm_type    NV        = G->numVertices;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    for (comm_type i = 0; i < NV; i++) {
        comm_type adj1 = vtxPtr[i];
        comm_type adj2 = vtxPtr[i+1];
        for(comm_type j=adj1; j<adj2; j++) {
            //		if( i < vtxInd[j].tail) {
            fprintf(out,"%ld %ld\n", i, vtxInd[j].tail);
            //		}
        }
    }
}

void displayGraphCharacteristics(graph *G) {
    printf("Within displayGraphCharacteristics()\n");
    long    sum = 0, sum_sq = 0;
    double  average, avg_sq, variance, std_dev;
    comm_type    maxDegree = 0;
    comm_type    isolated  = 0;
    comm_type    degreeOne = 0;
    comm_type    NS        = G->sVertices;
    comm_type    NV        = G->numVertices;
    comm_type    NT        = NV - NS;
    comm_type    NE        = G->numEdges;
    comm_type    *vtxPtr   = G->edgeListPtrs;
    comm_type    tNV       = NV; //Number of vertices
    
    if ( (NS == 0)||(NS == NV) ) {  //Nonbiparite graph
        for (comm_type i = 0; i < NV; i++) {
            comm_type degree = vtxPtr[i+1] - vtxPtr[i];
            sum_sq += degree*degree;
            sum    += degree;
            if (degree > maxDegree)
                maxDegree = degree;
            if ( degree == 0 )
                isolated++;
            if ( degree == 1 )
                degreeOne++;
        }
        average  = (double) sum / tNV;
        avg_sq   = (double) sum_sq / tNV;
        variance = avg_sq - (average*average);
        std_dev  = sqrt(variance);
        
        printf("*******************************************\n");
        printf("General Graph: Characteristics :\n");
        printf("*******************************************\n");
        printf("Number of vertices   :  %ld\n", NV);
        printf("Number of edges      :  %ld\n", NE);
        printf("Maximum out-degree is:  %ld\n", maxDegree);
        printf("Average out-degree is:  %lf\n",average);
        printf("Expected value of X^2:  %lf\n",avg_sq);
        printf("Variance is          :  %lf\n",variance);
        printf("Standard deviation   :  %lf\n",std_dev);
        printf("Isolated vertices    :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/tNV)*100);
        printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
        printf("Density              :  %lf%%\n",((double)NE/(NV*NV))*100);
        printf("*******************************************\n");
        
    }//End of nonbipartite graph
    else { //Bipartite graph
        
        //Compute characterisitcs from S side:
        for (comm_type i = 0; i < NS; i++) {
            comm_type degree = vtxPtr[i+1] - vtxPtr[i];
            sum_sq += degree*degree;
            sum    += degree;
            if (degree > maxDegree)
                maxDegree = degree;
            if ( degree == 0 )
                isolated++;
            if ( degree == 1 )
                degreeOne++;
        }
        average  = (double) sum / NS;
        avg_sq   = (double) sum_sq / NS;
        variance = avg_sq - (average*average);
        std_dev  = sqrt(variance);
        
        printf("*******************************************\n");
        printf("Bipartite Graph: Characteristics of S:\n");
        printf("*******************************************\n");
        printf("Number of S vertices :  %ld\n", NS);
        printf("Number of T vertices :  %ld\n", NT);
        printf("Number of edges      :  %ld\n", NE);
        printf("Maximum out-degree is:  %ld\n", maxDegree);
        printf("Average out-degree is:  %lf\n",average);
        printf("Expected value of X^2:  %lf\n",avg_sq);
        printf("Variance is          :  %lf\n",variance);
        printf("Standard deviation   :  %lf\n",std_dev);
        printf("Isolated (S)vertices :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/NS)*100);
        printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
        printf("Density              :  %lf%%\n",((double)NE/(NS*NS))*100);
        printf("*******************************************\n");
        
        sum = 0;
        sum_sq = 0;
        maxDegree = 0;
        isolated  = 0;
        //Compute characterisitcs from T side:
        for (comm_type i = NS; i < NV; i++) {
            comm_type degree = vtxPtr[i+1] - vtxPtr[i];
            sum_sq += degree*degree;
            sum    += degree;
            if (degree > maxDegree)
                maxDegree = degree;
            if ( degree == 0 )
                isolated++;
            if ( degree == 1 )
                degreeOne++;
        }
        
        average  = (double) sum / NT;
        avg_sq   = (double) sum_sq / NT;
        variance = avg_sq - (average*average);
        std_dev  = sqrt(variance);
        
        printf("Bipartite Graph: Characteristics of T:\n");
        printf("*******************************************\n");
        printf("Number of T vertices :  %ld\n", NT);
        printf("Number of S vertices :  %ld\n", NS);
        printf("Number of edges      :  %ld\n", NE);
        printf("Maximum out-degree is:  %ld\n", maxDegree);
        printf("Average out-degree is:  %lf\n",average);
        printf("Expected value of X^2:  %lf\n",avg_sq);
        printf("Variance is          :  %lf\n",variance);
        printf("Standard deviation   :  %lf\n",std_dev);
        printf("Isolated (T)vertices :  %ld (%3.2lf%%)\n", isolated, ((double)isolated/NT)*100);
        printf("Degree-one vertices  :  %ld (%3.2lf%%)\n", degreeOne, ((double)degreeOne/tNV)*100);
        printf("Density              :  %lf%%\n",((double)NE/(NT*NT))*100);
        printf("*******************************************\n");
    }//End of bipartite graph
}


//Convert a directed graph into an undirected graph:
//Parse through the directed graph and add edges in both directions
graph * convertDirected2Undirected(graph *G) {
    printf("Within convertDirected2Undirected()\n");
    int nthreads;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    
    double time1=0, time2=0, totalTime=0;
    //Get the iterators for the graph:
    comm_type NVer     = G->numVertices;
    comm_type NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    comm_type *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd  = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    printf("N= %ld  NE=%ld\n", NVer, NEdge);

    comm_type *degrees = (comm_type *) malloc ((NVer+1) * sizeof(comm_type));
    assert(degrees != NULL);
    
    //Count the edges from source --> sink (for sink > source)
#pragma omp parallel for
    for (comm_type v=0; v < NVer; v++ ) {
        comm_type adj1 = verPtr[v];
        comm_type adj2 = verPtr[v+1];
        for(comm_type k = adj1; k < adj2; k++ ) {
            comm_type w = verInd[k].tail;
            if (w < v)
                continue;
            __sync_fetch_and_add(&degrees[v+1], 1); //Increment by one to make space for zero
            __sync_fetch_and_add(&degrees[w+1], 1);
        }//End of for(k)
    }//End of for(i)
    
    //Build the pointer array:
    comm_type m=0;
    for (comm_type i=1; i<=NVer; i++) {
        m += degrees[i]; //Accumulate the number of edges
        degrees[i] += degrees[i-1];
    }
    //Sanity check:
    if(degrees[NVer] != 2*m) {
        printf("Number of edges added is not correct (%ld, %ld)\n", degrees[NVer], 2*m);
        exit(1);
    }
    printf("Done building pointer array\n");
    
    //Build CSR for Undirected graph:
    comm_type* counter = (comm_type *) malloc (NVer * sizeof(comm_type));
    assert(counter != NULL);
#pragma omp parallel for
    for (comm_type i=0; i<NVer; i++)
        counter[i] = 0;
    
    //Allocate memory for Edge list:
    edge *eList = (edge *) malloc ((2*m) * sizeof (edge));
    assert(eList != NULL);
    
#pragma omp parallel for
    for (comm_type v=0; v < NVer; v++ ) {
        comm_type adj1 = verPtr[v];
        comm_type adj2 = verPtr[v+1];
        for(comm_type k = adj1; k < adj2; k++ ) {
            comm_type w = verInd[k].tail;
            double weight = verInd[k].weight;
            if (w < v)
                continue;
            //Add edge v --> w
            comm_type location = degrees[v] + __sync_fetch_and_add(&counter[v], 1);
            if (location >= 2*m) {
                printf("location is out of bound: %ld \n", location);
                exit(1);
            }
            eList[location].head   = v;
            eList[location].tail   = w;
            eList[location].weight = weight;
            
            //Add edge w --> v
            location = degrees[w] + __sync_fetch_and_add(&counter[w], 1);
            if (location >= 2*m) {
                printf("location is out of bound: %ld \n", location);
                exit(1);
            }
            eList[location].head   = w;
            eList[location].tail   = v;
            eList[location].weight = weight;
        }//End of for(k)
    }//End of for(v)
    
    //Clean up:
    free(counter);
    
    //Build and return a graph data structure
    graph * Gnew = (graph *) malloc (sizeof(graph));
    Gnew->numVertices  = NVer;
    Gnew->sVertices    = NVer;
    Gnew->numEdges     = m;
    Gnew->edgeListPtrs = degrees;
    Gnew->edgeList     = eList;
    
    return Gnew;
    
}//End of convertDirected2Undirected()


comm_type removeEdges(comm_type NV, comm_type NE, edge *edgeList) {
    printf("Within removeEdges()\n");
    comm_type NGE = 0;
    comm_type *head = (comm_type *) malloc(NV * sizeof(comm_type));     /* head of linked list points to an edge */
    comm_type *next = (comm_type *) malloc(NE * sizeof(comm_type));     /* ptr to next edge in linked list       */
    
    /* Initialize linked lists */
    for (comm_type i = 0; i < NV; i++) head[i] = -1;
    for (comm_type i = 0; i < NE; i++) next[i] = -2;
    
    for (comm_type i = 0; i < NE; i++) {
        comm_type sv  = edgeList[i].head;
        comm_type ev  = edgeList[i].tail;
        if (sv == ev) continue;    /* self edge */
        comm_type * ptr = head + sv;     /* start at head of list for this key */
        while (1) {
            comm_type edgeId = *ptr;
            if (edgeId == -1) {         /* at the end of the list */
                edgeId = *ptr;             /* lock ptr               */
                if (edgeId == -1) {       /* if still end of list   */
                    comm_type newId = NGE;
                    NGE++;     /* increment number of good edges */
                    //edgeList[i].id = newId;                 /* set id of edge                 */
                    next[i] = -1;                           /* insert edge in linked list     */
                    *ptr = i;
                    break;
                }
                *ptr = edgeId;
            } else
                if (edgeList[edgeId].tail == ev) break;     /* duplicate edge */
                else  ptr = next + edgeId;
        }
    }
    /* Move good edges to front of edgeList                    */
    /* While edge i is a bad edge, swap with last edge in list */
    for (comm_type i = 0; i < NGE; i++) {
        while (next[i] == -2) {
            comm_type k = NE - 1;
            NE--;
            edgeList[i] = edgeList[k];
            next[i] = next[k];
        }
    }
    printf("About to free memory\n");
    free(head);
    free(next);
    printf("Exiting removeEdges()\n");
    return NGE;
}//End of removeEdges()

/* Since graph is undirected, sort each edge head --> tail AND tail --> head */
void SortEdgesUndirected(comm_type NV, comm_type NE, edge *list1, edge *list2, comm_type *ptrs) {
    for (comm_type i = 0; i < NV + 2; i++)
        ptrs[i] = 0;
    ptrs += 2;
    
    /* Histogram key values */
    for (comm_type i = 0; i < NE; i++) {
        ptrs[list1[i].head]++;
        ptrs[list1[i].tail]++;
    }
    /* Compute start index of each bucket */
    for (comm_type i = 1; i < NV; i++)
        ptrs[i] += ptrs[i-1];
    ptrs--;
    
    /* Move edges into its bucket's segment */
    for (comm_type i = 0; i < NE; i++) {
        comm_type head   = list1[i].head;
        comm_type index          = ptrs[head]++;
        //list2[index].id     = list1[i].id;
        list2[index].head   = list1[i].head;
        list2[index].tail   = list1[i].tail;
        list2[index].weight = list1[i].weight;
        
        comm_type tail   = list1[i].tail;
        index               = ptrs[tail]++;
        //list2[index].id     = list1[i].id;
        list2[index].head   = list1[i].tail;
        list2[index].tail   = list1[i].head;
        list2[index].weight = list1[i].weight;
    }
}//End of SortEdgesUndirected2()

/* Sort each node's neighbors by tail from smallest to largest. */
void SortNodeEdgesByIndex(comm_type NV, edge *list1, edge *list2, comm_type *ptrs) {
    for (comm_type i = 0; i < NV; i++) {
        edge *edges1 = list1 + ptrs[i];
        edge *edges2 = list2 + ptrs[i];
        comm_type size    = ptrs[i+1] - ptrs[i];
        
        /* Merge Sort */
        for (comm_type skip = 2; skip < 2 * size; skip *= 2) {
            for (comm_type sect = 0; sect < size; sect += skip)  {
                comm_type j = sect;
                comm_type l = sect;
                comm_type half_skip = skip / 2;
                comm_type k = sect + half_skip;
                
                comm_type j_limit = (j + half_skip < size) ? j + half_skip : size;
                comm_type k_limit = (k + half_skip < size) ? k + half_skip : size;
                
                while ((j < j_limit) && (k < k_limit)) {
                    if   (edges1[j].tail < edges1[k].tail) {edges2[l] = edges1[j]; j++; l++;}
                    else                                   {edges2[l] = edges1[k]; k++; l++;}
                }
                while (j < j_limit) {edges2[l] = edges1[j]; j++; l++;}
                while (k < k_limit) {edges2[l] = edges1[k]; k++; l++;}
            }
            edge *tmp = edges1;
            edges1 = edges2;
            edges2 = tmp;
        }
        // result is in list2, so move to list1
        if (edges1 == list2 + ptrs[i])
            for (comm_type j = ptrs[i]; j < ptrs[i+1]; j++) list1[j] = list2[j];
    }
}//End of SortNodeEdgesByIndex2()

//Build a reordering scheme for vertices based on clustering
//Group all communities together and number the vertices contiguously
//N = Number of vertices
//C = Community assignments for each vertex stored in an order
//old2NewMap = Stores the output of this routine
void buildOld2NewMap(comm_type N, comm_type *C, comm_type *commIndex) {
    printf("Within buildOld2NewMap(%ld) function...\n", N);
    printf("WARNING: Assumes that communities are numbered contiguously\n");
    assert(N > 0);
    //Compute number of communities:
    //Assume zero is a valid community id
    comm_type nC=-1;
    bool isZero = false;
    bool isNegative = false;
    for(comm_type i = 0; i < N; i++) {
        if(C[i] == 0)
            isZero = true; //Check if zero is a valid community
        if(C[i] < 0)
            isNegative = true; //Check if zero is a valid community
        if (C[i] > nC) {
            nC = C[i];
        }
    }
    printf("Largest community id observed: %ld\n", nC);
    if(isZero) {
        printf("Zero is a valid community id\n");
        nC++;
    }
    if(isNegative) {
        printf("Some vertices have not been assigned communities\n");
        nC++; //Place to store all the unassigned vertices
    }
    assert(nC>0);
    printf("Number of unique communities in C= %d\n", nC);
    
    //////////STEP 1: Create a CSR-like datastructure for communities in C
    comm_type * commPtr = (comm_type *) malloc ((nC+1) * sizeof(comm_type)); assert(commPtr != 0);
    //comm_type * commIndex = (comm_type *) malloc (N * sizeof(comm_type)); assert(commIndex != 0);
    comm_type * commAdded = (comm_type *) malloc (nC * sizeof(comm_type)); assert(commAdded != 0);
    
    // Initialization
#pragma omp parallel for
    for(comm_type i = 0; i < nC; i++) {
        commPtr[i] = 0;
        commAdded[i] = 0;
    }
    commPtr[nC] = 0;
    // Count the size of each community
#pragma omp parallel for
    for(comm_type i = 0; i < N; i++) {
        if(C[i] < 0) { //A negative value
            __sync_fetch_and_add(&commPtr[nC],1); //Unassigned vertices
        } else { //A positive value
            if(isZero)
                __sync_fetch_and_add(&commPtr[C[i]+1],1); //Zero-based indexing
            else
                __sync_fetch_and_add(&commPtr[C[i]],1); //One-based indexing
        }
    }//End of for(i)
    //Prefix sum:
    for(comm_type i=0; i<nC; i++) {
        commPtr[i+1] += commPtr[i];
    }
    //Group vertices with the same community in an order
#pragma omp parallel for
    for (comm_type i=0; i<N; i++) {
        comm_type tc = (comm_type)C[i];
        if(tc < 0) { //A negative value
            tc = nC-1;
            comm_type Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
            assert(Where < N);
            commIndex[Where] = i; //The vertex id
        } else {
        if(!isZero)
            tc--; //Convert to zero based index
        comm_type Where = commPtr[tc] + __sync_fetch_and_add(&(commAdded[tc]), 1);
        assert(Where < N);
        commIndex[Where] = i; //The vertex id
        }
    }
    printf("Done building structure for C...\n");
    
    //////////STEP 2: Create the old2New map:
    //This step will now be handled outside the routine
    /*
     #pragma omp parallel for
     for (comm_type i=0; i<N; i++) {
     old2NewMap[commIndex[i]] = i;
     //printf("(%ld) --> (%ld)\n", commIndex[i]+1, i+1);
     }
     */
    //Cleanup:
    free(commPtr); free(commAdded);
    //free(commIndex);
}//End of buildOld2NewMap()
