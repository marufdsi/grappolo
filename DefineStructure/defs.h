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

#ifndef _DEFS_H
#define _DEFS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <map>
#include <vector>
#include <unistd.h> //For getopts()

#include <immintrin.h>
#define MilanRealMax HUGE_VAL       // +INFINITY
#define MilanRealMin -MilanRealMax  // -INFINITY

#define PRINT_DETAILED_STATS_
//#define PRINT_TERSE_STATS_

typedef float f_weight;
typedef int comm_type;

/// 512 bit integer register initialize by all -1
const __m512 fl_set_neg_1 = _mm512_set1_ps(-1.0);
/// 512 bit integer register initialize by all 1
const __m512 fl_set_plus_1 = _mm512_set1_ps(1.0);
/// 512 bit integer register initialize by all 1
const __m512i set_plus_1 = _mm512_set1_epi32(1);
/// 512 bit integer register initialize by all -1
const __m512i set_minus_1 = _mm512_set1_epi32(-1);
/// 512 bit integer register initialize by bit wise 1
const __m512i set1 = _mm512_set1_epi32(0xFFFFFFFF);
/// 512 bit integer register initialize by all 0
const __m512i set0 = _mm512_set1_epi32(0x00000000);
/// 512 bit floating register initialize by all 0.0
const   __m512 fl_set0 = _mm512_set1_ps(0.0);

typedef struct comm
{
    comm_type size;
    f_weight degree;
}Comm;

typedef struct
{
    comm_type cid;       //community ID
    f_weight Counter; //Weight relative to that community
} mapElement;

typedef struct /* the edge data structure */
{
    comm_type head;
    comm_type tail;
    f_weight weight;
} edge;

typedef struct /* the graph data structure */
{
    comm_type numVertices;        /* Number of columns                                */
    comm_type sVertices;          /* Number of rows: Bipartite graph: number of S vertices; T = N - S */
    comm_type numEdges;           /* Each edge stored twice, but counted once        */
    comm_type * edgeListPtrs;     /* start vertex of edge, sorted, primary key        */
    edge * edgeList;         /* end   vertex of edge, sorted, secondary key      */
} graph;

typedef struct /* the graph data structure for directed graph */
{
    comm_type numVertices;        /* Number of vertices                                */
    comm_type numEdges;           /* Each edge is stored only once (u --> v)           */
    //Outgoing edges
    comm_type *edgeListPtrsOut;            /* Edge pointer vector O(|V|)    */
    edge * edgeListOut;        /* Edge weight vector O(|E|)     */
    //Incoming edges
    comm_type *edgeListPtrsIn;             /* Edge pointer vector O(|V|)    */
    edge * edgeListIn;         /* Edge weight vector O(|E|)     */
} dGraph;

struct clustering_parameters
{
    const char *inFile; //Input file
    int ftype;  //File type

    bool strongScaling; //Enable strong scaling
    bool output; //Printout the clustering data
    bool VF; //Vertex following turned on
    int coloring; // Type of coloring
    bool replaceMap; //If map data structure is replaced with a vector
    int numColors; // Type of coloring
    int syncType; // Type of synchronization method
    int basicOpt; //If map data structure is replaced with a vector
    bool threadsOpt;
    f_weight C_thresh; //Threshold with coloring on
    comm_type minGraphSize; //Min |V| to enable coloring
    f_weight threshold; //Value of threshold
    int percentage;
    bool compute_metrics;
    bool intrinsic;

    clustering_parameters();
    void usage();

    //Define in parseInputParameter.cpp
    bool parse(int argc, char *argv[]);
};

//Reverse Cuthill-McKee Algorithm
void algoReverseCuthillMcKee( graph *G, comm_type *pOrder, int nThreads );
void algoReverseCuthillMcKeeStrict( graph *G, comm_type *pOrder, int nThreads );

#endif
