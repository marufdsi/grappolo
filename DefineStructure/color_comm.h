#ifndef __color_comm__
#define __color_comm__

#include "basic_comm.h"
#include "coloring.h"

void runMultiPhaseColoring(graph *G, comm_type *C_orig, int coloring, int numColors, int replaceMap, comm_type minGraphSize,
            double threshold, double C_threshold, int numThreads, int threadsOpt);

double algoLouvainWithDistOneColoring(graph* G, comm_type *C, int nThreads, int* color,
			int numColor, double Lower, double thresh, double *totTime, int *numItr);

double algoLouvainWithDistOneColoringNoMap(graph* G, comm_type *C, int nThreads, int* color,
			int numColor, double Lower, double thresh, double *totTime, int *numItr);

#endif
