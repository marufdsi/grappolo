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

#ifndef __CLUSTERING__FUNCTIONS__
#define __CLUSTERING__FUNCTIONS__

#include "defs.h"

using namespace std;

void sumVertexDegree(edge* vtxInd, long* vtxPtr, double* vDegree, long NV, Comm* cInfo);

/// Single floating point version
void sumVertexDegree_sfp(edge* vtxInd, long* vtxPtr, f_weight* vDegree, long NV, Comm* cInfo);

/// Single floating point vectorized version
void sumVertexDegreeVec_sfp(edge* vtxInd, comm_type* vtxPtr, f_weight* vDegree, comm_type NV, comm_type* cInfo_size,
                            f_weight* cInfo_degree);
double calConstantForSecondTerm(double* vDegree, long NV);

/// Single floating point version
f_weight calConstantForSecondTerm_sfp(f_weight* vDegree, long NV);

void initCommAss(long* pastCommAss, long* currCommAss, long NV);

/// Single floating point version
void initCommAss_SFP(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV);

/// Single floating point version
void initCommAssOpt(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV,
		    mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
		    Comm* cInfo, f_weight constant, f_weight* vDegree );

void initCommAssOpt_SFP(long* pastCommAss, long* currCommAss, long NV,
		    mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
		    Comm* cInfo, double constant, double* vDegree );

void initCommAssOptVec_SFP(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV,
                           comm_type* cid, f_weight* Counter, comm_type* vtxPtr, comm_type* head, comm_type* tail,
                           f_weight* weights, comm_type* size, f_weight* degree, f_weight constant, f_weight* vDegree );

double buildLocalMapCounter(long adj1, long adj2, map<long, long> &clusterLocalMap, 
						  vector<double> &Counter, edge* vtxInd, long* currCommAss, long me);
/// Single floating point version
f_weight buildLocalMapCounter_sfp(comm_type adj1, comm_type adj2, map<comm_type, comm_type> &clusterLocalMap,
						  vector<f_weight> &Counter, edge* vtxInd, comm_type* currCommAss, comm_type me);

double buildLocalMapCounterNoMap(long v, mapElement* clusterLocalMap, long* vtxPtr, edge* vtxInd,
                               long* currCommAss, long &numUniqueClusters);
/// Single floating point version
f_weight buildLocalMapCounterNoMap_SFP(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
                                     comm_type* currCommAss, comm_type &numUniqueClusters);
/// Single floating point version
f_weight buildLocalMapCounterVec_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, comm_type* head,
                                     comm_type* tail, f_weight* weights, comm_type* currCommAss,
                                     comm_type &numUniqueClusters);

long max(map<long, long> &clusterLocalMap, vector<double> &Counter,
		 double selfLoop, Comm* cInfo, double degree, long sc, double constant ) ;

/// Single floating point version
comm_type max_sfp(map<long, long> &clusterLocalMap, vector<f_weight> &Counter,
         f_weight selfLoop, Comm* cInfo, f_weight degree, long sc, f_weight constant ) ;

long maxNoMap(long v, mapElement* clusterLocalMap, long* vtxPtr, double selfLoop, Comm* cInfo, double degree,
              long sc, double constant, long numUniqueClusters );
/// Single floating point version
comm_type maxNoMap_SFP(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, f_weight selfLoop, Comm* cInfo,
        f_weight degree, comm_type sc, f_weight constant, comm_type numUniqueClusters );
/// Single floating point vectorized version
comm_type maxNoMapVec_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, f_weight selfLoop,
                          comm_type * cInfo_size, f_weight* cInfo_degree, f_weight degree, comm_type sc, f_weight constant,
                          comm_type numUniqueClusters );

void computeCommunityComparisons(vector<long>& C1, long N1, vector<long>& C2, long N2);

double computeGiniCoefficient(long *colorSize, int numColors);
double computeMerkinMetric(long* C1, long N1, long* C2, long N2);
double computeVanDongenMetric(long* C1, long N1, long* C2, long N2);

//Sorting functions:
void merge(long* arr, long l, long m, long r);
void mergeSort(long* arr, long l, long r);
void SortNeighborListUsingInsertionAndMergeSort(graph *G);
long removeEdges(long NV, long NE, edge *edgeList);
void SortEdgesUndirected(long NV, long NE, edge *list1, edge *list2, long *ptrs);
void SortNodeEdgesByIndex(long NV, edge *list1, edge *list2, long *ptrs);

double* computeEdgeSimilarityMetrics(graph *G);
graph* buildSparifiedGraph(graph *Gin, double alpha);

void buildOld2NewMap(long N, long *C, long *commIndex); //Build the reordering map

#endif
