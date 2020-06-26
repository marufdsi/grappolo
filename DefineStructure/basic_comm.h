#ifndef __BASIC__COMM__
#define __BASIC__COMM__

// Define in louvainMultiPhaseRun.cpp
void runMultiPhaseBasic(graph *G, comm_type *C_orig, int basicOpt, comm_type minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt, char * graphName = "");

/// Single floating point version
void runMultiPhaseBasic_sfp(graph *G, comm_type *C_orig, int basicOpt, comm_type minGraphSize,
                            f_weight threshold, f_weight C_threshold, int numThreads, int threadsOpt, char * graphName = "");

// same as above, but runs exactly one phase
void runMultiPhaseBasicOnce(graph *G, comm_type *C_orig, int basicOpt, comm_type minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt);

// uses Granell, Arenas, et al. Fast track resistance
void runMultiPhaseBasicFastTrackResistance(graph *G, comm_type *C_orig, int basicOpt, comm_type minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt);


void runMultiPhaseBasicApprox(graph *G, comm_type *C_orig, int basicOpt, comm_type minGraphSize,
			double threshold, double C_threshold, int numThreads, int threadsOpt, int percentage);

// Define in parallelLouvianMethod.cpp
double parallelLouvianMethod(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr, bool *change);

/// Single floating point vectorized version
f_weight vectorizedLouvianMethod(graph *G, comm_type *C, int nThreads, f_weight Lower,
                               f_weight thresh, double *totTime, int *numItr, bool *change);
f_weight parallelLouvianMethod_SFP(graph *G, comm_type *C, int nThreads, f_weight Lower,
                               f_weight thresh, double *totTime, int *numItr, bool *change);

// Define in parallelLouvianMethodApprox.cpp
double parallelLouvianMethodApprox(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr, int percentage);

double parallelLouvianMethodNoMap(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);
				
double parallelLouvianMethodScale(graph *G, comm_type *C, int nThreads, double Lower,
				double thresh, double *totTime, int *numItr);

// implements Granell, Arenas, et al. Fast track resistance
// Granell, Clara, Sergio Gomez, and Alex Arenas. "Hierarchical multiresolution method to 
// overcome the resolution limit in complex networks." International Journal of Bifurcation 
// and Chaos 22, no. 07 (2012): 1250171.

// Define in parallelLouvianMethodFastTrackResistance.cpp
double parallelLouvianMethodFastTrackResistance(graph *G, comm_type *C, int nThreads, double Lower,
        double thresh, double *totTime, int *numItr, int phase, double* rmin, double* finMod);

// Define in parallelLouvianMethodNoMapFastTrackResistance.cpp
double parallelLouvianMethodNoMapFastTrackResistance(graph *G, comm_type *C, int nThreads, double Lower,
        double thresh, double *totTime, int *numItr, int phase, double* rmin, double* finMod);
				
// Define in parallelLouvianMethodScaleFastTrackResistance.cpp
double parallelLouvianMethodScaleFastTrackResistance(graph *G, comm_type *C, int nThreads, double Lower,
        double thresh, double *totTime, int *numItr, int phase, double* rmin, double* finMod);

#endif
