#include "input_output.h"
#include "defs.h"
#include "sstream"
#include "utilityStringTokenizer.hpp"

void parse_EdgeListBinaryNew(graph * G, char *fileName) {
  printf("Parsing a file in binary format...\n");
  printf("WARNING: Assumes that the graph is undirected -- every edge is stored twice.\n");
  int nthreads = 0;

  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  
  double time1, time2;

  std::ifstream ifs;  
  ifs.open(fileName, std::ifstream::in | std::ifstream::binary);
  if (!ifs) {
    std::cerr << "Error opening binary format file: " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  comm_type NV, NE, weighted;
  //Parse line-1: #Vertices #Edges
  ifs.read(reinterpret_cast<char*>(&NV), sizeof(NV));
  ifs.read(reinterpret_cast<char*>(&NE), sizeof(NE));
  ifs.read(reinterpret_cast<char*>(&weighted), sizeof(weighted));

  comm_type* verPtrRaw = (comm_type*) malloc( (NV+1)*sizeof(comm_type)); assert(verPtrRaw != 0);
  edge* edgeListRaw = (edge*) malloc(2*NE*sizeof(edge)); assert(edgeListRaw != 0);

  ifs.read(reinterpret_cast<char*>(verPtrRaw), sizeof(comm_type) * (NV+1));
  ifs.read(reinterpret_cast<char*>(edgeListRaw), sizeof(edge) * (2*NE));
 
  ifs.close(); //Close the file

  G->sVertices    = NV;
  G->numVertices  = NV;
  G->numEdges     = NE;
  G->edgeListPtrs = verPtrRaw;
  G->edgeList     = edgeListRaw;
  
  //Clean up

  //displayGraph(G);
}//End of parse_Dimacs9FormatDirectedNewD()
