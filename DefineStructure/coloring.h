#ifndef __coloring__
#define __coloring__

#include "basic_util.h"
#include "coloringUtils.h"

// In coloringDistanceOne.cpp
int algoDistanceOneVertexColoringOpt(graph *G, int *vtxColor, int nThreads, double *totTime);
int algoDistanceOneVertexColoring(graph *G, int *vtxColor, int nThreads, double *totTime);

// In ColoringMultiHasMaxMin.cpp
int algoColoringMultiHashMaxMin(graph *G, int *vtxColor, int nThreads, double *totTime, int nHash, int nItrs);

// In vBase.cpp
int vBaseRedistribution(graph* G, int* vtxColor, int ncolors, int type);

// In equtiableColoringDistanceOne.cpp
void buildColorSize(comm_type NVer, int *vtxColor, int numColors, comm_type *colorSize);
void computeVariance(comm_type NVer, int numColors, comm_type *colorSize);

void equitableDistanceOneColorBased(graph *G, int *vtxColor, int numColors, comm_type *colorSize,
				    int nThreads, double *totTime, int type);

int algoColoringMultiHashMaxMin(graph *G, int *vtxColor, int nThreads, double *totTime, int nHash, int nItrs);


#endif
