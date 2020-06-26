#ifndef __sync_comm__
#define __sync_comm__

#include "basic_util.h"
#include "utilityClusteringFunctions.h"

void runMultiPhaseSyncType(graph *G, comm_type *C_orig, int syncType, comm_type minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt);

double parallelLouvainMethodFullSyncEarly(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr,int ytype, int freedom);
				
double parallelLouvainMethodFullSync(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr,int ytype, int freedom);
				
double parallelLouvianMethodEarlyTerminate(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);
				
// Define in fullSyncUtility.cpp
double buildAndLockLocalMapCounter(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd,
                               comm_type* currCommAss, comm_type &numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double& eix, int freedom);

void maxAndFree(comm_type v, mapElement* clusterLocalMap, comm_type* vtxPtr, edge* vtxInd, double selfLoop, Comm* cInfo, comm_type* CA,
							double constant, comm_type numUniqueClusters, omp_lock_t* vlocks, omp_lock_t* clocks, int ytype, double eix, double* vDegree);

#endif
