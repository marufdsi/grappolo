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

#include "utilityClusteringFunctions.h"
#include "../DefineStructure/defs.h"
#include <immintrin.h>
#include <math.h>
using namespace std;

void updateAxForOpt(Comm* cInfo, comm_type* currCommAss, double* vDegree, comm_type NV)
{
  //printf("NUMBER OF VERTICES: %d\n", NV);
  #pragma omp parallel for
  for(comm_type i = 0; i < NV; i++)
  {
    if(currCommAss[i] != i){
/*    __sync_fetch_and_sub(&(cInfo[i].degree),vDegree[i]);
      __sync_fetch_and_sub(&(cInfo[i].size),1);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].degree),vDegree[i]);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].size),1);
*/
      #pragma omp atomic update
      cInfo[i].degree -= vDegree[i];
      #pragma omp atomic update
      cInfo[i].size -= 1;
      #pragma omp atomic update
      cInfo[currCommAss[i]].degree += vDegree[i];
      #pragma omp atomic update
      cInfo[currCommAss[i]].size += 1;
   }
  }
}

void updateAxForOpt_SFP(Comm* cInfo, comm_type* currCommAss, f_weight* vDegree, comm_type NV)
{
  //printf("NUMBER OF VERTICES: %d\n", NV);
  #pragma omp parallel for
  for(comm_type i = 0; i < NV; i++)
  {
    if(currCommAss[i] != i){
/*    __sync_fetch_and_sub(&(cInfo[i].degree),vDegree[i]);
      __sync_fetch_and_sub(&(cInfo[i].size),1);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].degree),vDegree[i]);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].size),1);
*/
      #pragma omp atomic update
      cInfo[i].degree -= vDegree[i];
      #pragma omp atomic update
      cInfo[i].size -= 1;
      #pragma omp atomic update
      cInfo[currCommAss[i]].degree += vDegree[i];
      #pragma omp atomic update
      cInfo[currCommAss[i]].size += 1;
   }
  }
}

void updateAxForOptVec_SFP(comm_type* size, f_weight* degree, comm_type* currCommAss, f_weight* vDegree, comm_type NV)
{
  //printf("NUMBER OF VERTICES: %d\n", NV);
  #pragma omp parallel for
  for(comm_type i = 0; i < NV; i++)
  {
    if(currCommAss[i] != i){
/*    __sync_fetch_and_sub(&(cInfo[i].degree),vDegree[i]);
      __sync_fetch_and_sub(&(cInfo[i].size),1);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].degree),vDegree[i]);
      __sync_fetch_and_add(&(cInfo[currCommAss[i]].size),1);
*/
      #pragma omp atomic update
      degree[i] -= vDegree[i];
      #pragma omp atomic update
      size[i] -= 1;
      #pragma omp atomic update
      degree[currCommAss[i]] += vDegree[i];
      #pragma omp atomic update
      size[currCommAss[i]] += 1;
   }
  }
}
void sumVertexDegree(edge* vtxInd, comm_type* vtxPtr, double* vDegree, comm_type NV, Comm* cInfo) {
#pragma omp parallel for
  for (comm_type i=0; i<NV; i++) {
    comm_type adj1 = vtxPtr[i];	    //Begin
    comm_type adj2 = vtxPtr[i+1];	//End
    double totalWt = 0;
    for(comm_type j=adj1; j<adj2; j++) {
      totalWt += vtxInd[j].weight;
    }
    vDegree[i] = totalWt;	//Degree of each node
    cInfo[i].degree = totalWt;	//Initialize the community
    cInfo[i].size = 1;
  }
}//End of sumVertexDegree()

/// Single floating point version
void sumVertexDegree_sfp(edge* vtxInd, comm_type* vtxPtr, f_weight* vDegree, comm_type NV, Comm* cInfo) {
#pragma omp parallel for
  for (comm_type i=0; i<NV; i++) {
      comm_type adj1 = vtxPtr[i];	    //Begin
      comm_type adj2 = vtxPtr[i+1];	//End
      f_weight totalWt = 0;
    for(comm_type j=adj1; j<adj2; j++) {
      totalWt += vtxInd[j].weight;
    }
    vDegree[i] = totalWt;	//Degree of each node
    cInfo[i].degree = totalWt;	//Initialize the community
    cInfo[i].size = 1;
  }
}//End of sumVertexDegree()

/// Single floating point vectorized version
void sumVertexDegreeVec_sfp(edge* vtxInd, comm_type* vtxPtr, f_weight* vDegree, comm_type NV, comm_type* cInfo_size,
        f_weight* cInfo_degree) {
#pragma omp parallel for
  for (comm_type i=0; i<NV; i++) {
      comm_type adj1 = vtxPtr[i];	    //Begin
      comm_type adj2 = vtxPtr[i+1];	//End
      f_weight totalWt = 0;
    for(comm_type j=adj1; j<adj2; j++) {
      totalWt += vtxInd[j].weight;
    }
    vDegree[i] = totalWt;	//Degree of each node
    cInfo_degree[i] = totalWt;	//Initialize the community
    cInfo_size[i] = 1;
  }
}//End of sumVertexDegree()

double calConstantForSecondTerm(double* vDegree, comm_type NV) {
  double totalEdgeWeightTwice = 0;
  #pragma omp parallel for reduction(+:totalEdgeWeightTwice)
  for (comm_type i=0; i<NV; i++) {
      totalEdgeWeightTwice += vDegree[i];
  }
  return (double)1/totalEdgeWeightTwice;
}//End of calConstantForSecondTerm()

/// Single floating point version
f_weight calConstantForSecondTerm_sfp(f_weight* vDegree, comm_type NV) {
    f_weight totalEdgeWeightTwice = 0;
  #pragma omp parallel for reduction(+:totalEdgeWeightTwice)
  for (comm_type i=0; i<NV; i++) {
      totalEdgeWeightTwice += vDegree[i];
  }
  return (f_weight)1/totalEdgeWeightTwice;
}//End of calConstantForSecondTerm()

void initCommAss(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV) {
#pragma omp parallel for
  for (comm_type i=0; i<NV; i++) {
    pastCommAss[i] = i; //Initialize each vertex to its cluster
    currCommAss[i] = i;
  }
}//End of initCommAss()

/// Single floating point version
void initCommAss_SFP(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV) {
#pragma omp parallel for
  for (comm_type i=0; i<NV; i++) {
    pastCommAss[i] = i; //Initialize each vertex to its cluster
    currCommAss[i] = i;
  }
}//End of initCommAss()

//Smart initialization assuming that each vertex is assigned to its own cluster
//WARNING: Will ignore duplicate edge entries (multi-graph)
void initCommAssOpt(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV,
		    mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
		    Comm* cInfo, double constant, double* vDegree ) {

#pragma omp parallel for
  for (comm_type v=0; v<NV; v++) {
    comm_type adj1  = vtxPtr[v];
    comm_type adj2  = vtxPtr[v+1];
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v
    
    pastCommAss[v] = v; //Initialize each vertex to its own cluster
    //currCommAss[v] = v; //Initialize with a self cluster
    
    //Step-1: Build local map counter (without a map):
    comm_type numUniqueClusters = 0;
    double selfLoop = 0;
    clusterLocalMap[sPosition].cid     = v; //Add itself
    clusterLocalMap[sPosition].Counter = 0; //Initialize the count
    numUniqueClusters++;
    //Parse through the neighbors
    for(comm_type j=adj1; j<adj2; j++) {
      if(vtxInd[j].tail == v) {	// SelfLoop need to be recorded
	      selfLoop += (comm_type)vtxInd[j].weight;
        clusterLocalMap[sPosition].Counter = vtxInd[j].weight; //Initialize the count
        continue;
      }
      //Assume each neighbor is assigned to a separate cluster
      //Assume no duplicates (only way to improve performance at this step)
      clusterLocalMap[sPosition + numUniqueClusters].cid     = vtxInd[j].tail; //Add the cluster id (initialized to itself)
      clusterLocalMap[sPosition + numUniqueClusters].Counter = vtxInd[j].weight; //Initialize the count
      numUniqueClusters++;
    }//End of for(j)
    
    //Step 2: Find max:
    comm_type maxIndex = v;	//Assign the initial value as the current community
    double curGain = 0;
    double maxGain = 0;
    double eix = clusterLocalMap[sPosition].Counter - selfLoop; //NOT SURE ABOUT THIS.
    double ax  = cInfo[v].degree - vDegree[v];
    double eiy = 0;
    double ay  = 0;    
    for(comm_type k=0; k<numUniqueClusters; k++) {
      if(v != clusterLocalMap[sPosition + k].cid) {
        ay = cInfo[clusterLocalMap[sPosition + k].cid].degree; // degree of cluster y
        eiy = clusterLocalMap[sPosition + k].Counter; 	//Total edges incident on cluster y
        curGain = 2*(eiy - eix) - 2*vDegree[v]*(ay - ax)*constant;
	
        if( (curGain > maxGain) || ((curGain==maxGain) && (curGain != 0) && (clusterLocalMap[sPosition + k].cid < maxIndex)) ) {
          maxGain  = curGain;
          maxIndex = clusterLocalMap[sPosition + k].cid;
        }
      }
    }//End of for()
    
    if(cInfo[maxIndex].size == 1 && cInfo[v].size == 1 && maxIndex > v) { //Swap protection
      maxIndex = v;
    }    
    currCommAss[v] = maxIndex; //Assign the new community
  }

  updateAxForOpt(cInfo,currCommAss,vDegree,NV);
}//End of initCommAssOpt()

//Smart initialization assuming that each vertex is assigned to its own cluster
//WARNING: Will ignore duplicate edge entries (multi-graph)
void initCommAssOpt_SFP(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV,
		    mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
		    Comm* cInfo, f_weight constant, f_weight* vDegree ) {

#pragma omp parallel for
  for (comm_type v=0; v<NV; v++) {
      comm_type adj1  = vtxPtr[v];
      comm_type adj2  = vtxPtr[v+1];
      comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v

    pastCommAss[v] = v; //Initialize each vertex to its own cluster
    //currCommAss[v] = v; //Initialize with a self cluster

    //Step-1: Build local map counter (without a map):
      comm_type numUniqueClusters = 0;
      f_weight selfLoop = 0;
    clusterLocalMap[sPosition].cid     = v; //Add itself
    clusterLocalMap[sPosition].Counter = 0; //Initialize the count
    numUniqueClusters++;
    //Parse through the neighbors
    for(comm_type j=adj1; j<adj2; j++) {
      if(vtxInd[j].tail == v) {	// SelfLoop need to be recorded
	      selfLoop += (comm_type)vtxInd[j].weight;
        clusterLocalMap[sPosition].Counter = vtxInd[j].weight; //Initialize the count
        continue;
      }
      //Assume each neighbor is assigned to a separate cluster
      //Assume no duplicates (only way to improve performance at this step)
      clusterLocalMap[sPosition + numUniqueClusters].cid     = vtxInd[j].tail; //Add the cluster id (initialized to itself)
      clusterLocalMap[sPosition + numUniqueClusters].Counter = vtxInd[j].weight; //Initialize the count
      numUniqueClusters++;
    }//End of for(j)

    //Step 2: Find max:
      comm_type maxIndex = v;	//Assign the initial value as the current community
      f_weight curGain = 0;
      f_weight maxGain = 0;
      f_weight eix = clusterLocalMap[sPosition].Counter - selfLoop; //NOT SURE ABOUT THIS.
      f_weight ax  = cInfo[v].degree - vDegree[v];
      f_weight eiy = 0;
      f_weight ay  = 0;
    for(comm_type k=0; k<numUniqueClusters; k++) {
      if(v != clusterLocalMap[sPosition + k].cid) {
        ay = cInfo[clusterLocalMap[sPosition + k].cid].degree; // degree of cluster y
        eiy = clusterLocalMap[sPosition + k].Counter; 	//Total edges incident on cluster y
        curGain = 2*(eiy - eix) - 2*vDegree[v]*(ay - ax)*constant;

        if( (curGain > maxGain) || ((curGain==maxGain) && (curGain != 0) && (clusterLocalMap[sPosition + k].cid < maxIndex)) ) {
          maxGain  = curGain;
          maxIndex = clusterLocalMap[sPosition + k].cid;
        }
      }
    }//End of for()

    if(cInfo[maxIndex].size == 1 && cInfo[v].size == 1 && maxIndex > v) { //Swap protection
      maxIndex = v;
    }
    currCommAss[v] = maxIndex; //Assign the new community
  }

    updateAxForOpt_SFP(cInfo,currCommAss,vDegree,NV);
}//End of initCommAssOpt()

//Smart initialization assuming that each vertex is assigned to its own cluster
//WARNING: Will ignore duplicate edge entries (multi-graph)
void initCommAssOptVec_SFP(comm_type* pastCommAss, comm_type* currCommAss, comm_type NV,
                           comm_type* cid, f_weight* Counter, comm_type* vtxPtr, comm_type* head, comm_type* tail,
                           f_weight* weights, comm_type* size, f_weight* degree, f_weight constant, f_weight* vDegree ) {

#pragma omp parallel for
  for (comm_type v=0; v<NV; v++) {
      comm_type adj1  = vtxPtr[v];
      comm_type adj2  = vtxPtr[v+1];
      comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v

    pastCommAss[v] = v; //Initialize each vertex to its own cluster
    //currCommAss[v] = v; //Initialize with a self cluster

    //Step-1: Build local map counter (without a map):
      comm_type numUniqueClusters = 0;
      f_weight selfLoop = 0;
    cid[sPosition]     = v; //Add itself
    Counter[sPosition] = 0; //Initialize the count
    numUniqueClusters++;
    //Parse through the neighbors
    for(comm_type j=adj1; j<adj2; j++) {
      if(tail[j] == v) {	// SelfLoop need to be recorded
	      selfLoop += (comm_type)weights[j];
        Counter[sPosition] = weights[j]; //Initialize the count
        continue;
      }
      //Assume each neighbor is assigned to a separate cluster
      //Assume no duplicates (only way to improve performance at this step)
      cid[sPosition + numUniqueClusters] = tail[j]; //Add the cluster id (initialized to itself)
      Counter[sPosition + numUniqueClusters] = weights[j]; //Initialize the count
      numUniqueClusters++;
    }//End of for(j)

    //Step 2: Find max:
      comm_type maxIndex = v;	//Assign the initial value as the current community
      f_weight curGain = 0;
      f_weight maxGain = 0;
      f_weight eix = Counter[sPosition] - selfLoop; //NOT SURE ABOUT THIS.
      f_weight ax  = degree[v] - vDegree[v];
      f_weight eiy = 0;
      f_weight ay  = 0;
    for(comm_type k=0; k<numUniqueClusters; k++) {
      if(v != cid[sPosition + k]) {
        ay = degree[cid[sPosition + k]]; // degree of cluster y
        eiy = Counter[sPosition + k]; 	//Total edges incident on cluster y
        curGain = 2*(eiy - eix) - 2*vDegree[v]*(ay - ax)*constant;

        if( (curGain > maxGain) || ((curGain==maxGain) && (curGain != 0) && (cid[sPosition + k] < maxIndex)) ) {
          maxGain  = curGain;
          maxIndex = cid[sPosition + k];
        }
      }
    }//End of for()

    if(size[maxIndex] == 1 && size[v] == 1 && maxIndex > v) { //Swap protection
      maxIndex = v;
    }
    currCommAss[v] = maxIndex; //Assign the new community
  }

    updateAxForOptVec_SFP(size, degree, currCommAss, vDegree, NV);
}//End of initCommAssOpt()


double buildLocalMapCounter(comm_type adj1, comm_type adj2, map<comm_type, comm_type> &clusterLocalMap,
			 vector<double> &Counter, edge* vtxInd, comm_type* currCommAss, comm_type me) {
  
  map<comm_type, comm_type>::iterator storedAlready;
  comm_type numUniqueClusters = 1;
  double selfLoop = 0;
  for(comm_type j=adj1; j<adj2; j++) {
    if(vtxInd[j].tail == me) {	// SelfLoop need to be recorded
      selfLoop += vtxInd[j].weight;
    }
    
    storedAlready = clusterLocalMap.find(currCommAss[vtxInd[j].tail]); //Check if it already exists
    if( storedAlready != clusterLocalMap.end() ) {	//Already exists
      Counter[storedAlready->second]+= vtxInd[j].weight; //Increment the counter with weight
    } else {
      clusterLocalMap[currCommAss[vtxInd[j].tail]] = numUniqueClusters; //Does not exist, add to the map
      Counter.push_back(vtxInd[j].weight); //Initialize the count
      numUniqueClusters++;
    }
  }//End of for(j)

  return selfLoop;
}//End of buildLocalMapCounter()


f_weight buildLocalMapCounter_sfp(comm_type adj1, comm_type adj2, map<comm_type, comm_type> &clusterLocalMap,
			 vector<f_weight> &Counter, edge* vtxInd, comm_type* currCommAss, comm_type me) {

  map<comm_type, comm_type>::iterator storedAlready;
  comm_type numUniqueClusters = 1;
    f_weight selfLoop = 0;
  for(comm_type j=adj1; j<adj2; j++) {
    if(vtxInd[j].tail == me) {	// SelfLoop need to be recorded
      selfLoop += vtxInd[j].weight;
    }

    storedAlready = clusterLocalMap.find(currCommAss[vtxInd[j].tail]); //Check if it already exists
    if( storedAlready != clusterLocalMap.end() ) {	//Already exists
      Counter[storedAlready->second]+= vtxInd[j].weight; //Increment the counter with weight
    } else {
      clusterLocalMap[currCommAss[vtxInd[j].tail]] = numUniqueClusters; //Does not exist, add to the map
      Counter.push_back(vtxInd[j].weight); //Initialize the count
      numUniqueClusters++;
    }
  }//End of for(j)

  return selfLoop;
}//End of buildLocalMapCounter()

//Build the local-map data structure using vectors
double buildLocalMapCounterNoMap(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
                               comm_type* currCommAss, comm_type &numUniqueClusters) {
    comm_type adj1  = vtxPtr[v];
    comm_type adj2  = vtxPtr[v+1];
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v

    comm_type storedAlready = 0;
    double selfLoop = 0;
    for(comm_type j=adj1; j<adj2; j++) {
        if(vtxInd[j].tail == v) {	// SelfLoop need to be recorded
            selfLoop += vtxInd[j].weight;
        }
        bool storedAlready = false; //Initialize to zero
        for(comm_type k=0; k<numUniqueClusters; k++) { //Check if it already exists
            if(currCommAss[vtxInd[j].tail] ==  clusterLocalMap[sPosition+k].cid) {
                storedAlready = true;
                clusterLocalMap[sPosition + k].Counter += vtxInd[j].weight; //Increment the counter with weight
                break;
            }
        }
        if( storedAlready == false ) {	//Does not exist, add to the map
            clusterLocalMap[sPosition + numUniqueClusters].cid     = currCommAss[vtxInd[j].tail];
            clusterLocalMap[sPosition + numUniqueClusters].Counter = vtxInd[j].weight; //Initialize the count
            numUniqueClusters++;
        }
    }//End of for(j)
    return selfLoop;
}//End of buildLocalMapCounter()

//Build the local-map data structure using vectors
f_weight buildLocalMapCounterNoMap_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, comm_type* head,
                                       comm_type* tail, f_weight* weights, comm_type* currCommAss,
                                       comm_type &numUniqueClusters) {
    comm_type adj1  = vtxPtr[v];
    comm_type adj2  = vtxPtr[v+1];
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v

    comm_type storedAlready = 0;
    f_weight selfLoop = 0;
    for(comm_type j=adj1; j<adj2; j++) {
        if(tail[j] == v) {	// SelfLoop need to be recorded
            selfLoop += weights[j];
        }
        bool storedAlready = false; //Initialize to zero
        for(comm_type k=0; k<numUniqueClusters; k++) { //Check if it already exists
            if(currCommAss[tail[j]] ==  cid[sPosition+k]) {
                storedAlready = true;
                Counter[sPosition + k] += weights[j]; //Increment the counter with weight
                break;
            }
        }
        if( storedAlready == false ) {	//Does not exist, add to the map
            cid[sPosition + numUniqueClusters]     = currCommAss[tail[j]];
            Counter[sPosition + numUniqueClusters] = weights[j]; //Initialize the count
            numUniqueClusters++;
        }
    }//End of for(j)
    return selfLoop;
}//End of buildLocalMapCounter()

//Build the local-map data structure using vectors
f_weight buildLocalMapCounterVec_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, comm_type* head,
                                     comm_type* tail, f_weight* weights, comm_type* currCommAss,
                                     comm_type *numUniqueClusters) {
    comm_type adj1  = vtxPtr[v];
    comm_type adj2  = vtxPtr[v+1];
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v
    /// 512 bit integer register initialize by all -1
    const   __m512 fl_set_neg_1 = _mm512_set1_ps(-1.0);
    /// 512 bit integer register initialize by all 1
    const   __m512 fl_set_plus_1 = _mm512_set1_ps(1.0);
    /// 512 bit integer register initialize by all 1
    const   __m512i set_plus_1 = _mm512_set1_epi32(1);
    /// 512 bit integer register initialize by all -1
    const   __m512i set1 = _mm512_set1_epi32(0xFFFFFFFF);
    /// 512 bit integer register initialize by all 0
    const   __m512i set0 = _mm512_set1_epi32(0x00000000);
    comm_type storedAlready = 0;
    f_weight selfLoop = 0;
    comm_type vector_op = (adj2-adj1)/16;
    cout << "vector_op : " << vector_op << " adj1: " << adj1 << " adj2: " << adj2 << " start: " << vector_op *16 << endl;
//    cout << "Neighbors: " << adj2-adj1  << " numUniqueClusters: " << (*numUniqueClusters) << endl;
    /// perform intrinsic on the neighbors that are multiple of 16
    const   __m512i check_self_loop = _mm512_set1_epi32(v);
//    cout << "Perform vector operation" << endl;
    for(comm_type j=adj1; j<(vector_op * 16); j+=16) {
        /// Load at most 16 tail of the neighbors.
        __m512i tail_vec = _mm512_loadu_si512((__m512i *) &tail[j]);
        /// Load at most 16 neighbors edge weight.
        __m512 w_vec = _mm512_loadu_ps((__m512 *) &weights[j]);
        /// Mask to find u != v
        __mmask16 self_loop_mask = _mm512_cmpeq_epi32_mask(check_self_loop, tail_vec);
        selfLoop += _mm512_mask_reduce_add_ps(self_loop_mask, w_vec);
        /*if(tail[j] == v) {	// SelfLoop need to be recorded
            selfLoop += weights[j];
        }*/
        __m512i currCommAss_vec = _mm512_i32gather_epi32(tail_vec, &currCommAss[0], 4);
        bool storedAlready = false; //Initialize to zero
        int count_existing_cluster = 0;
        __mmask16 comm_mask = pow(2, 16) - 1;
        for(comm_type k=0; k<(*numUniqueClusters); k++) { //Check if it already exists
            const   __m512i check_existing_cluster = _mm512_set1_epi32(cid[sPosition+k]);
            __mmask16 existing_cluster_mask = _mm512_cmpeq_epi32_mask(check_existing_cluster, currCommAss_vec);
            comm_mask = _mm512_kand(comm_mask, _mm512_knot(existing_cluster_mask));
            Counter[sPosition + k] += _mm512_mask_reduce_add_ps(existing_cluster_mask, w_vec);
            count_existing_cluster += _mm_popcnt_u32((unsigned) existing_cluster_mask);
            if(count_existing_cluster >= 16){
                break;
            }
            /*if(currCommAss[tail[j]] ==  cid[sPosition+k]) {
                storedAlready = true;
                Counter[sPosition + k] += weights[j]; //Increment the counter with weight
                break;
            }*/
        }

        if(count_existing_cluster < 16) {
            __m512i C_conflict = _mm512_mask_conflict_epi32(set_plus_1, comm_mask, currCommAss_vec);
            /// Calculate mask using NAND of C_conflict and set1
            const __mmask16 mask = _mm512_kand(_mm512_testn_epi32_mask(C_conflict, set1), comm_mask);
            __m512i distinct_comm;
            /// It will find out the distinct community.
            distinct_comm = _mm512_mask_compress_epi32(set0, mask, currCommAss_vec);
            comm_type *remaining_comm = (comm_type *) &distinct_comm;
            for (int k = 0; k < _mm_popcnt_u32((unsigned) mask); ++k) {
                const __m512i comm = _mm512_set1_epi32(remaining_comm[k]);
                __mmask16 comm_mask = _mm512_cmpeq_epi32_mask(comm, currCommAss_vec);
                Counter[sPosition + (*numUniqueClusters)] += _mm512_mask_reduce_add_ps(comm_mask, w_vec);
                cid[sPosition + (*numUniqueClusters)] = remaining_comm[k];
                (*numUniqueClusters)++;
            }
            /*if (storedAlready == false) {    //Does not exist, add to the map
                cid[sPosition + numUniqueClusters] = currCommAss[tail[j]];
                Counter[sPosition + numUniqueClusters] = weights[j]; //Initialize the count
                numUniqueClusters++;
            }*/
        }
    }//End of for(j)
//    cout << "End vector operation" << endl;
    for(comm_type j=(vector_op * 16); j<adj2; j++) {
        if(tail[j] == v) {	// SelfLoop need to be recorded
            selfLoop += weights[j];
        }
        bool storedAlready = false; //Initialize to zero
        for(comm_type k=0; k<(*numUniqueClusters); k++) { //Check if it already exists
            if(currCommAss[tail[j]] ==  cid[sPosition+k]) {
                storedAlready = true;
                Counter[sPosition + k] += weights[j]; //Increment the counter with weight
                break;
            }
        }
        if( storedAlready == false ) {	//Does not exist, add to the map
            cid[sPosition + (*numUniqueClusters)]     = currCommAss[tail[j]];
            Counter[sPosition + (*numUniqueClusters)] = weights[j]; //Initialize the count
            (*numUniqueClusters)++;
        }
    }//End of for(j)
    cout << "[" << (adj2-adj1) << "] numUniqueClusters: " << (*numUniqueClusters) << endl;
    return selfLoop;
}//End of buildLocalMapCounter()
                                                                                
comm_type max(map<comm_type, comm_type> &clusterLocalMap, vector<double> &Counter,
         double selfLoop, Comm* cInfo, double degree, comm_type sc, double constant ) {
    
    map<comm_type, comm_type>::iterator storedAlready;
    comm_type maxIndex = sc;	//Assign the initial value as self community
    double curGain = 0;
    double maxGain = 0;
    double eix = Counter[0] - selfLoop;
    double ax = cInfo[sc].degree - degree;
    double eiy = 0;
    double ay = 0;
    
    storedAlready = clusterLocalMap.begin();
    do {
        if(sc != storedAlready->first) {
            ay = cInfo[storedAlready->first].degree; // degree of cluster y
            eiy = Counter[storedAlready->second]; 	//Total edges incident on cluster y
            curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;
            
            if( (curGain > maxGain) ||
               ((curGain==maxGain) && (curGain != 0) && (storedAlready->first < maxIndex)) ) {
                maxGain = curGain;
                maxIndex = storedAlready->first;
            }
        }
        storedAlready++; //Go to the next cluster
    } while ( storedAlready != clusterLocalMap.end() );
    
    if(cInfo[maxIndex].size == 1 && cInfo[sc].size ==1 && maxIndex > sc) { //Swap protection
        maxIndex = sc;
    }
    
    return maxIndex;		
}//End max()

comm_type max_sfp(map<comm_type, comm_type> &clusterLocalMap, vector<f_weight> &Counter,
             f_weight selfLoop, Comm* cInfo, f_weight degree, comm_type sc, f_weight constant ) {

    map<comm_type, comm_type>::iterator storedAlready;
    comm_type maxIndex = sc;	//Assign the initial value as self community
    f_weight curGain = 0;
    f_weight maxGain = 0;
    f_weight eix = Counter[0] - selfLoop;
    f_weight ax = cInfo[sc].degree - degree;
    f_weight eiy = 0;
    f_weight ay = 0;

    storedAlready = clusterLocalMap.begin();
    do {
        if(sc != storedAlready->first) {
            ay = cInfo[storedAlready->first].degree; // degree of cluster y
            eiy = Counter[storedAlready->second]; 	//Total edges incident on cluster y
            curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;

            if( (curGain > maxGain) ||
               ((curGain==maxGain) && (curGain != 0) && (storedAlready->first < maxIndex)) ) {
                maxGain = curGain;
                maxIndex = storedAlready->first;
            }
        }
        storedAlready++; //Go to the next cluster
    } while ( storedAlready != clusterLocalMap.end() );

    if(cInfo[maxIndex].size == 1 && cInfo[sc].size ==1 && maxIndex > sc) { //Swap protection
        maxIndex = sc;
    }

    return maxIndex;
}//End max()



comm_type maxNoMap(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, double selfLoop, Comm* cInfo, double degree,
              comm_type sc, double constant, comm_type numUniqueClusters ) {
                                                                                
    comm_type maxIndex = sc;	//Assign the initial value as the current community
    double curGain = 0;
    double maxGain = 0;
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v
    double eix = clusterLocalMap[sPosition].Counter - selfLoop;
    double ax  = cInfo[sc].degree - degree;
    double eiy = 0;
    double ay  = 0;
    
    for(comm_type k=0; k<numUniqueClusters; k++) {
        if(sc != clusterLocalMap[sPosition + k].cid) {
            ay = cInfo[clusterLocalMap[sPosition + k].cid].degree; // degree of cluster y
            eiy = clusterLocalMap[sPosition + k].Counter; 	//Total edges incident on cluster y
            curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;

            if( (curGain > maxGain) ||
               ((curGain==maxGain) && (curGain != 0) && (clusterLocalMap[sPosition + k].cid < maxIndex)) ) {
                maxGain  = curGain;
                maxIndex = clusterLocalMap[sPosition + k].cid;
            }
        }
    }//End of for()

    if(cInfo[maxIndex].size == 1 && cInfo[sc].size ==1 && maxIndex > sc) { //Swap protection
        maxIndex = sc;
    }

    return maxIndex;
}//End maxNoMap()

comm_type maxNoMap_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, f_weight selfLoop,
                       comm_type * cInfo_size, f_weight* cInfo_degree, f_weight degree, comm_type sc, f_weight constant,
                       comm_type numUniqueClusters ) {

    comm_type maxIndex = sc;	//Assign the initial value as the current community
    f_weight curGain = 0;
    f_weight maxGain = 0;
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v
    f_weight eix = Counter[sPosition] - selfLoop;
    f_weight ax  = cInfo_degree[sc] - degree;
    f_weight eiy = 0;
    f_weight ay  = 0;

    for(comm_type k=0; k<numUniqueClusters; k++) {
        if(sc != cid[sPosition + k]) {
            ay = cInfo_degree[cid[sPosition + k]]; // degree of cluster y
            eiy = Counter[sPosition + k]; 	//Total edges incident on cluster y
            curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;

            if( (curGain > maxGain) ||
                ((curGain==maxGain) && (curGain != 0) && (cid[sPosition + k] < maxIndex)) ) {
                maxGain  = curGain;
                maxIndex = cid[sPosition + k];
            }
        }
    }//End of for()

    if(cInfo_size[maxIndex] == 1 && cInfo_size[sc] ==1 && maxIndex > sc) { //Swap protection
        maxIndex = sc;
    }

    return maxIndex;
}//End maxNoMap()

comm_type maxNoMapVec_SFP(comm_type v, comm_type *cid, f_weight *Counter, comm_type* vtxPtr, f_weight selfLoop,
        comm_type * cInfo_size, f_weight* cInfo_degree, f_weight degree, comm_type sc, f_weight constant,
        comm_type numUniqueClusters ) {

    comm_type maxIndex = sc;	//Assign the initial value as the current community
    f_weight curGain = 0;
    f_weight maxGain = 0;
    comm_type sPosition = vtxPtr[v]+v; //Starting position of local map for v
    f_weight eix = Counter[sPosition] - selfLoop;
    f_weight ax  = cInfo_degree[sc] - degree;
    f_weight eiy = 0;
    f_weight ay  = 0;

    for(comm_type k=0; k<numUniqueClusters; k++) {
        if(sc != cid[sPosition + k]) {
            ay = cInfo_degree[cid[sPosition + k]]; // degree of cluster y
            eiy = Counter[sPosition + k]; 	//Total edges incident on cluster y
            curGain = 2*(eiy - eix) - 2*degree*(ay - ax)*constant;

            if( (curGain > maxGain) ||
               ((curGain==maxGain) && (curGain != 0) && (cid[sPosition + k] < maxIndex)) ) {
                maxGain  = curGain;
                maxIndex = cid[sPosition + k];
            }
        }
    }//End of for()

    if(cInfo_size[maxIndex] == 1 && cInfo_size[sc] ==1 && maxIndex > sc) { //Swap protection
        maxIndex = sc;
    }

    return maxIndex;
}//End maxNoMap()
            
    
  
