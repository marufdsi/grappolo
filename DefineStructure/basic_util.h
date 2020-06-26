#ifndef __UTILITLY__
#define __UTILITLY__

// Define in buildNextPhase.cpp
comm_type renumberClustersContiguously(comm_type *C, comm_type size);
double buildNextLevelGraphOpt(graph *Gin, graph *Gout, comm_type *C, comm_type numUniqueClusters, int nThreads);
double buildNextLevelGraphOpt_SFP(graph *Gin, graph *Gout, comm_type *C, comm_type numUniqueClusters, int nThreads);
void buildNextLevelGraph(graph *Gin, graph *Gout, comm_type *C, comm_type numUniqueClusters);
comm_type buildCommunityBasedOnVoltages(graph *G, comm_type *Volts, comm_type *C, comm_type *Cvolts);
void segregateEdgesBasedOnVoltages(graph *G, comm_type *Volts);
inline void Visit(comm_type v, comm_type myCommunity, short *Visited, comm_type *Volts,
				  comm_type* vtxPtr, edge* vtxInd, comm_type *C);
inline void Visit(comm_type v, comm_type myCommunity, short *Visited, comm_type *Volts,
                  comm_type* vtxPtr, edge* vtxInd, comm_type *C);
inline void Visit(comm_type v, comm_type myCommunity, comm_type *Visited, comm_type *Volts,
                  comm_type* vtxPtr, edge* vtxInd, comm_type *C);
inline void Visit(comm_type v, comm_type myCommunity, short *Visited, comm_type *Volts,
                  comm_type* vtxPtr, edge* vtxInd, comm_type *C);
				  
// Define in vertexFollowing.cpp
comm_type vertexFollowing(graph *G, comm_type *C);
double buildNewGraphVF(graph *Gin, graph *Gout, comm_type *C, comm_type numUniqueClusters);

// Define in utilityFunctions.cpp
double computeGiniCoefficient(comm_type *colorSize, int numColors);
void generateRandomNumbers(double *RandVec, comm_type size);
void displayGraph(graph *G);
void duplicateGivenGraph(graph *Gin, graph *Gout);
void displayGraphEdgeList(graph *G);
void writeEdgeListToFile(graph *G, FILE* out);
void displayGraphCharacteristics(graph *G);


#endif
