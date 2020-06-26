#include "input_output.h"
void writeGraphMetisSimpleFormat(graph* G, char *filename) {
    //Get the iterators for the graph:
    comm_type NVer     = G->numVertices;
    comm_type NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    comm_type *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    
    printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);
    printf("Writing graph in Metis format - each edge represented twice -- no weights; 1-based indices\n");
    printf("Graph will be stored in file: %s\n", filename);
    
    FILE *fout;
    fout = fopen(filename, "w");
    if (!fout) {
        printf("Could not open the file \n");
        exit(1);
    }
    //First Line: #Vertices #Edges
    fprintf(fout, "%ld %ld\n", NVer, NEdge);
    //Write the edges:
    for (comm_type v=0; v<NVer; v++) {
        comm_type adj1 = verPtr[v];
        comm_type adj2 = verPtr[v+1];
        for(comm_type k = adj1; k < adj2; k++ ) {
            fprintf(fout, "%ld ", (verInd[k].tail+1) );
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    printf("Graph has been stored in file: %s\n",filename);
}//End of writeGraphPajekFormat()

//Output the graph in Pajek format:
//Each edge is represented twice
void writeGraphPajekFormat(graph* G, char *filename) {
    //Get the iterators for the graph:
    comm_type NVer     = G->numVertices;
    comm_type NS       = G->sVertices;
    comm_type NT       = NVer - NS;
    comm_type NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    comm_type *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);
    
    printf("Writing graph in Pajek format - Undirected graph - each edge represented ONLY ONCE!\n");
    printf("Graph will be stored in file: %s\n", filename);
    
    
    FILE *fout;
    fout = fopen(filename, "w");
    if (!fout) {
        printf("Could not open the file \n");
        exit(1);
    }
    //First Line: Vertices
    fprintf(fout, "*Vertices %ld\n", NVer);
    for (comm_type i=0; i<NVer; i++) {
        fprintf(fout, "%ld\n", i+1);
    }
    
    //Write the edges:
    fprintf(fout, "*Edges %ld\n", NEdge);
    for (comm_type v=0; v<NVer; v++) {
        comm_type adj1 = verPtr[v];
        comm_type adj2 = verPtr[v+1];
        //Edge lines: <adjacent> <weight>
        for(comm_type k = adj1; k < adj2; k++ ) {
            if (v <= verInd[k].tail) { //Print only once
                fprintf(fout, "%ld %ld %g\n", v+1, (verInd[k].tail+1), (verInd[k].weight) );
            }
        }
    }
    fclose(fout);
    printf("Graph has been stored in file: %s\n",filename);
}//End of writeGraphPajekFormat()

//Output the graph in Pajek format:
//Each edge is represented twice
void writeGraphPajekFormatWithCommunityInfo(graph* G, char *filename, comm_type *C) {
    //Get the iterators for the graph:
    comm_type NVer     = G->numVertices;
    comm_type NS       = G->sVertices;
    comm_type NT       = NVer - NS;
    comm_type NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    comm_type *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    
    printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);
    printf("Writing graph in Pajek format - Undirected graph - each edge represented ONLY ONCE!\n");
    printf("Graph will be stored in file: %s\n", filename);
    
    FILE *fout;
    fout = fopen(filename, "w");
    if (!fout) {
        printf("Could not open the file \n");
        exit(1);
    }
    //First Line: Vertices
    fprintf(fout, "*Vertices %ld\n", NVer);
    for (comm_type i=0; i<NVer; i++) {
        fprintf(fout, "%ld  \"%ld\"\n", i+1, C[i]);
    }
    //Write the edges:
    fprintf(fout, "*Edges %ld\n", NEdge);
    for (comm_type v=0; v<NVer; v++) {
        comm_type adj1 = verPtr[v];
        comm_type adj2 = verPtr[v+1];
        //Edge lines: <adjacent> <weight>
        for(comm_type k = adj1; k < adj2; k++ ) {
            if (v <= verInd[k].tail) { //Print only once
                fprintf(fout, "%ld %ld %g\n", v+1, (verInd[k].tail+1), (verInd[k].weight) );
            }
        }
    }
    fclose(fout);
    printf("Graph has been stored in file: %s\n", filename);
}//End of writeGraphPajekFormat()
