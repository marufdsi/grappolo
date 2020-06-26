// ***********************************************************************
//
//     Rundemanen: CUDA C++ parallel program for community detection
//   Md Naim (naim.md@gmail.com), Fredrik Manne (Fredrik.Manne@uib.no)
//                       University of Bergen
//
// ***********************************************************************
//
//       Copyright (2016) University of Bergen
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

#ifndef GRAPHGPU_H
#define	GRAPHGPU_H

#include "thrust/device_vector.h"
#include "iostream"
#include"commonconstants.h"

struct GraphGPU {
    // Graph 
    unsigned int nb_nodes;
    unsigned comm_type nb_links;
    double total_weight;

    int type;

    //thrust::device_vector<unsigned comm_type> degrees;

    thrust::device_vector<int> indices;

    thrust::device_vector<unsigned int> links;
    thrust::device_vector<float> weights;
    thrust::device_vector<int> colors;
    void greedyColoring(unsigned int wrpSz);

    //unsigned int nb_neighbors(unsigned int node);
    //double weighted_degree(unsigned int node);


};


#endif	/* GRAPHGPU_H */

