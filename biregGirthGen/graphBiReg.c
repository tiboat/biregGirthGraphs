#include "../nauty/nauty.h"
#include "bitset/bitset.h"


/*
 * This file contains code related to a graph datastructure.
 * It involves
 *  - creation, deletion of the graph datastructure,
 *  - printing different parts of it,
 *  - checking conditions on it,
 *  - adding, deleting edges,
 *  - adapting other parts of the datastructure.
 */



// Maximum number of edges with graph of MAXN vertices
#define MAXEDGES (MAXN*(MAXN-1) >> 1)

// Returns the minimum of two values
#define MIN(a,b) (((a)<(b))?(a):(b))



/************************************************
 * Data structures
 ***********************************************/

// Struct representing an edge between two vertices v1 and v2
struct edge { int v1, v2; };

// Data structure for a graph
// Based on GenK2
struct graph {
    graph nautyGraph[MAXN*MAXM];
    graph gCan[MAXN*MAXM];
    int numberOfVertices;
    int numEdgesAfterTree; // Number of edges added after construction of initial tree
    bitset* adjacencyList;
    bitset* verticesOfDeg;
    int minDeg; // Minimum degree of the degree set
    int maxDeg; // Maximum degree of the degree set
    bitset indexMinDeg; // Vertices of degree < minDeg
    bitset indexMaxDeg; // Vertices of degree < maxDeg (yet vertices of the Moore tree that are
    int distMatrix[MAXN][MAXN]; // 2D |V|x|V| distance matrix
                                // Only lower triangular part is used: distMatrix[i][j] with i < j
    struct edge addedEdges[MAXEDGES]; // edges added after initial construction
    bitset* legalChoices; // legalChoices[i] contains the vertices i is allowed to connect to
    int heightMooreTree;
    bitset* verticesAtLevel; // verticesAtLevel[l] contains the vertices at level l in the Moore tree
    int levelOfVertex[MAXN]; // levelOfVertex[i] equals the level of vertex i  in the Moore tree
};


// Data structures to do breadth-first-search and hold old legalChoices
// Source: EGR (https://github.com/JorikJooken/edgeGirthRegularGraphs)
int distBFS1[MAXN];
int distBFS2[MAXN];

int bfsQueue[MAXN];
int startQueue;
int endQueue;

#define MAXR 15
bitset oldLegalChoices[MAXR*MAXN][MAXN];
bitset oldLegalChoicesAddEdge[MAXR*MAXN][MAXN]; // separate bitsets, because otherwise might copy wrong things
int oldDist[MAXR*MAXN][MAXN][MAXN];
int endpoint1[MAXR*MAXN][MAXR*MAXN];
int endpoint2[MAXR*MAXN][MAXR*MAXN];



/************************************************
 * Macros for creating, checking and changing graph
 ***********************************************/


//  Initializer for empty graph.
// Based on GenK2
#define emptyGraph(g,minDegree,maxDegree) EMPTYGRAPH((g)->nautyGraph, (g)->numberOfVertices, MAXM);\
 (g)->verticesOfDeg[0] = compl(EMPTY,(g)->numberOfVertices);\
 for(int i = 1; i <= maxDegree; i++) { (g)->verticesOfDeg[i] = EMPTY; }\
 (g)->minDeg = minDegree;\
 (g)->maxDeg = maxDegree;\
 (g)->indexMinDeg = compl(EMPTY, numVertices);\
 (g)->indexMaxDeg = compl(EMPTY, numVertices);\
 for(int i = 0; i < (g)->numberOfVertices; i++) {(g)->adjacencyList[i] = EMPTY;\
 }

//  Add one edge.
// Based on GenK2
#define addEdge(g,i,j) {ADDONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);       \
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) - 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) - 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);                   \
 if(deg((g), i) == (g)->minDeg) removeElement((g)->indexMinDeg, i);\
 else if(deg((g), i) == (g)->maxDeg) removeElement((g)->indexMaxDeg, i);\
 if(deg((g), j) == (g)->minDeg) removeElement((g)->indexMinDeg, j);\
 else if(deg((g), j) == (g)->maxDeg) removeElement((g)->indexMaxDeg, j);\
 }

//  Remove one edge.
// Based on GenK2
#define removeEdge(g,i,j) {REMOVEONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 removeElement((g)->adjacencyList[i], j); removeElement((g)->adjacencyList[j],i); \
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) + 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) + 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j); \
 if(deg((g), i) == (g)->minDeg-1) add((g)->indexMinDeg, i);\
 else if(deg((g), i) == (g)->maxDeg-1) add((g)->indexMaxDeg, i);\
 if(deg((g), j) == (g)->minDeg-1) add((g)->indexMinDeg, j);\
 else if(deg((g), j) == (g)->maxDeg-1) add((g)->indexMaxDeg, j);\
}

// From GenK2
#define REMOVEONEEDGE(g,i,j,MAXM)\
 DELELEMENT(GRAPHROW(g,i,MAXM),j); DELELEMENT(GRAPHROW(g,j,MAXM),i)

// deg(g,v) is the degree of vertex v in graph g
#define deg(g,v) size(g->adjacencyList[v])

// Use this macro when you are not sure which of the two vertices i or j is smaller
#define dist(g,i,j) (i < j ? (g)->distMatrix[i][j] : (g)->distMatrix[j][i])



/************************************************
 * Printing functions
 ***********************************************/

// From GenK2
void printAdjacencyList(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr,"%d:", i);
        forEach(neighbour, g->adjacencyList[i]) {
            fprintf(stderr," %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");
}

void printLegalChoicesMaxDeg(struct graph *g) {
    fprintf(stderr,"Legal choices\n");
    forEach(v, g->indexMaxDeg) {
        fprintf(stderr,"%d:", v);
        forEach(vChoice, g->legalChoices[v]) {
            fprintf(stderr," %d", vChoice);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

void printLegalChoicesMinDeg(struct graph *g) {
    fprintf(stderr,"Legal choices\n");
    forEach(v, g->indexMinDeg) {
        fprintf(stderr,"%d:", v);
        forEach(vChoice, g->legalChoices[v]) {
            fprintf(stderr," %d", vChoice);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

void printIndexMaxDeg(struct graph *g) {
    fprintf(stderr,"indexMaxDeg\n");
    forEach(v, g->indexMaxDeg) {
        fprintf(stderr,"%d ", v);
    }
    fprintf(stderr,"\n");
}

void printIndexMinDeg(struct graph *g) {
    fprintf(stderr,"indexMinDeg\n");
    forEach(v, g->indexMinDeg) {
        fprintf(stderr,"%d ", v);
    }
    fprintf(stderr,"\n");
}

void printVerticesOfDeg(struct graph *g) {
    fprintf(stderr,"verticesOfDeg\n");
    for(int i = 0; i <= g->maxDeg; i++) {
        fprintf(stderr,"%d:", i);
        forEach(v, g->verticesOfDeg[i]) {
            fprintf(stderr," %d", v);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

void printVerticesAtLevel(struct graph *g) {
    fprintf(stderr,"verticesAtLevel\n");
    for(int i = 0; i < g->heightMooreTree; i++) {
        fprintf(stderr,"%d:", i);
        forEach(v, g->verticesAtLevel[i]) {
            fprintf(stderr," %d", v);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

void printDistMatrixIndexMaxDeg(struct graph *g) {
    fprintf(stderr,  "Indexed Dist matrix \n");
    forEach(j, g->indexMaxDeg) {
        for(int i = j+1; i < g->numberOfVertices; i++) {
            if (g->distMatrix[j][i] != g->numberOfVertices) {
                fprintf(stderr, " (%d,%d) dist %d\n", j, i, g->distMatrix[j][i]);
            }
        }
    }
    fprintf(stderr, "\n");
}



/************************************************
 * Nauty-like functions
 ***********************************************/

// Based on GenK2
int gCanToG6(char* graphString, graph gCan[], int numberOfVertices) {
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) ((numberOfVertices >> i*6) & 0x3f) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) ((numberOfVertices >> i*6) & 0x3f) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(ISELEMENT(GRAPHROW(gCan, i, MAXM), j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline character.
    graphString[pointer++] = '\n';
    return pointer;
}


double timeCreateCanonical = 0;

// From Genk2
//  Uses Nauty to get canonical form.
void createCanonicalForm(struct graph *g, graph gCan[]) {
    clock_t start = clock();
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    statsblk stats;

    densenauty(g->nautyGraph, lab, ptn, orbits, &options, &stats, MAXM, g->numberOfVertices, gCan);

    clock_t end = clock();
    timeCreateCanonical += (double)(end - start) / CLOCKS_PER_SEC;
}


//  Write gCan to stdout in graph6 format.
void writeToG6(graph gCan[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int ptr = gCanToG6(graphString, gCan, numberOfVertices);
    graphString[ptr] = '\0';
    printf("%s", graphString);
}


/************************************************
 * Boolean functions checking graph properties
 ***********************************************/

// Returns true if all vertices have degree d1 or d2; otherwise returns false
bool isBiRegular(struct graph *g, int d1, int d2) {
    return size(g->verticesOfDeg[d1]) + size(g->verticesOfDeg[d2]) == g->numberOfVertices;
}


// Pre-condition: v1 < v2
// Returns true if
//  - v1 and v2 are both no isolated vertices
//  - v1 is an isolated vertex, v2 not and v1 is the smallest isolated vertex
//  - v1 is not an isolated vertex, v2 is and v2 is the smallest isolated vertex
//  - v1 and v2 are both isolated vertices and v1 is the smallest isolated vertex and v2 is the
//    second smallest isolated vertex
bool isolatedLexMin(struct graph *g, int v1, int v2) {
    bitset isolatedVertices = g->verticesOfDeg[0];
    // if v1 is an isolated vertex, check if all other isolated vertices
    // are larger than v1
    if(contains(isolatedVertices, v1)) {
        removeElement(isolatedVertices, v1);
        forEach(v, isolatedVertices) {
            if(v < v1)
                return false;
        }
    }
    // if v2 is an isolated vertex, check if all other isolated vertices
    // are larger than v2, except for v1.
    if(contains(isolatedVertices, v2)) {
        forEach(v, isolatedVertices) {
            if(v < v2)
                return false;
        }
    }
    return true;
}


// Checks if the vertices have enough legal choices, more specifically:
//  - if deg(v) < minDeg, then need at least minDeg - deg(v) legal choices
//  - if deg(v) > minDeg and deg(v) < maxDeg, then need at least maxDeg - deg(v) legal choices
//  - else, nothing is checked
bool enoughChoicesLeftVertex(struct graph *g, int v) {
    if ((deg(g, v) < g->minDeg && size(g->legalChoices[v]) < g->minDeg - deg(g, v)) ||
        (deg(g, v) > g->minDeg && deg(g, v) < g->maxDeg && size(g->legalChoices[v]) < g->maxDeg - deg(g, v))) {
        return false;
    }
    return true;
}

// Returns true if all vertices of graph g have at least degree g->minDeg
bool allVerticesAtLeastMinDeg(struct graph *g) {
    return isEmpty(g->indexMinDeg);
}



/************************************************
 * Functions for changing the graph
 ***********************************************/

// Removes vertex v from all vertices in indexMaxDeg
// and leaves v with no legal choices
void makeVertexIllegal(struct graph *g, int v) {
    forEach(vIndex, g->indexMaxDeg) {
        removeElement(g->legalChoices[vIndex], v);
    }
    g->legalChoices[v] = EMPTY;
}

// Removes vertex v from all vertices in indexMaxDeg, leaves v with no legal choices,
// Registers which vertices' legalChoices changed in legalChoicesChanged and
// sets oldLegalChoices to it's legalChoices if the vertex was not yet registered in
// legalChoicesChanged.
// Returns false if removing one of the legalChoices leads to a vertex not being able
// to complete to its necessary degree. Otherwise, returns true.
bool makeVertexIllegal2(struct graph *g, int v, int searchDepth, bitset *legalChoicesChanged) {
    forEach(vLegal, g->legalChoices[v]) {
        if(!contains(*legalChoicesChanged, vLegal)) {
            add(*legalChoicesChanged, vLegal);
            oldLegalChoices[searchDepth][vLegal] = g->legalChoices[vLegal];
        }
        removeElement(g->legalChoices[vLegal], v);
        if(!enoughChoicesLeftVertex(g, vLegal))
            return false;
    }
    if(!contains(*legalChoicesChanged, v)) {
        add(*legalChoicesChanged, v);
        oldLegalChoices[searchDepth][v] = g->legalChoices[v];
    }
    g->legalChoices[v] = EMPTY;
    return true;
}

// Removes vertex v from all vertices in indexMaxDeg, leaves v with no legal choices,
// Registers which vertices' legalChoices changed in legalChoicesChanged and
// sets oldLegalChoicesAddEdge to it's legalChoices if the vertex was not yet registered in
// legalChoicesChanged.
// Returns false if removing one of the legalChoices leads to a vertex not being able
// to complete to its necessary degree. Otherwise, returns true.
// This function is used in comparison to makeVertexIllegal2, when making
// vertices illegal when in the process of adding an edge.
bool makeVertexIllegal2AddEdge(struct graph *g, int v, int searchDepth, bitset *legalChoicesChanged) {
    forEach(vLegal, g->legalChoices[v]) {
        if(!contains(*legalChoicesChanged, vLegal)) {
            add(*legalChoicesChanged, vLegal);
            oldLegalChoicesAddEdge[searchDepth][vLegal] = g->legalChoices[vLegal];
        }
        removeElement(g->legalChoices[vLegal], v);
        if(!enoughChoicesLeftVertex(g, vLegal))
            return false;
    }
    if(!contains(*legalChoicesChanged, v)) {
        add(*legalChoicesChanged, v);
        oldLegalChoicesAddEdge[searchDepth][v] = g->legalChoices[v];
    }
    g->legalChoices[v] = EMPTY;
    return true;
}

// Removes vertex v1 as a legal choice for v2 and vice versa.
// Registers which vertices' legalChoices changed in legalChoicesChanged and
// sets oldLegalChoicesAddEdge to it's legalChoices if the vertex was not yet registered in
// legalChoicesChanged.
// Returns false if removing one of the legalChoices leads to a vertex not being able
// to complete to its necessary degree. Otherwise, returns true.
bool makeEdgeIllegal(struct graph *g, int v1, int v2, int searchDepth, bitset *legalChoicesChanged) {
    if(!contains(*legalChoicesChanged, v1)) {
        add(*legalChoicesChanged, v1);
        oldLegalChoices[searchDepth][v1] = g->legalChoices[v1];
    }
    removeElement(g->legalChoices[v1], v2);
    if(!enoughChoicesLeftVertex(g, v1))
        return false;

    if(!contains(*legalChoicesChanged, v2)) {
        add(*legalChoicesChanged, v2);
        oldLegalChoices[searchDepth][v2] = g->legalChoices[v2];
    }
    removeElement(g->legalChoices[v2], v1);
    return enoughChoicesLeftVertex(g, v2);
}


// Returns a vertex in indexMinDeg with the least amount of legalChoices larger than zero.
// If there is no such vertex, -1 is returned.
int getVertexMinLegalChoiceMinDeg(struct graph *g) {
    int vMin = -1;
    int minChoices = g->numberOfVertices;
    forEach(v, g->indexMinDeg) {
        if (size(g->legalChoices[v]) < minChoices && size(g->legalChoices[v]) != 0) {
            minChoices = size(g->legalChoices[v]);
            vMin = v;
        }
    }
    return vMin;
}

// Returns a vertex in indexMaxDeg with the least amount of legalChoices larger than zero.
// If there is no such vertex, -1 is returned.
int getVertexMinLegalChoiceMaxDeg(struct graph *g) {
    int vMin = -1;
    int minChoices = g->numberOfVertices;
    forEach(v, g->indexMaxDeg) {
        if (size(g->legalChoices[v]) < minChoices && size(g->legalChoices[v]) != 0) {
            minChoices = size(g->legalChoices[v]);
            vMin = v;
        }
    }
    return vMin;
}

// Returns a vertex vMin in indexMaxDeg with degree larger than minDeg with the least amount of legal choices and
// the amount of legal choices must be at least the amount of edges it needs to complete vMin to degree maxDeg.
// If there is no such vertex, -1 is returned.
int getVertexBetweenMinMaxDeg(struct graph *g) {
    // Also, gives the one with the least amount of legal choices and
    // doens't give one that has no legal choices left
    int vMin = -1;
    int minChoices = g->numberOfVertices;
    forEach(v, g->indexMaxDeg) {
        if (deg(g, v) > g->minDeg && size(g->legalChoices[v]) >= g->maxDeg - deg(g,v) && size(g->legalChoices[v]) < minChoices) {
            minChoices = size(g->legalChoices[v]);
            vMin = v;
        }
    }
    return vMin;
}


// Pre-condition: deg(v1) < g->maxDeg and deg(v2) < g->maxDeg
// Returns true if adding the edge was successful;
// returns false if adding the edge led to a vertex having not enough
// legalChoices left to complete to its full degree or to having two vertices of maxDeg
// being at a distance smaller than minDistMaxDeg
// This function is used when constructing the Moore tree.
void addEdgeWithDistMatrixLegal(struct graph *g, int minDistMaxDeg, int girth, int v1, int v2) {
    int potDist1, potDist2, minPotDist;

    addEdge(g, v1, v2);

    for(int vi1 = 0; vi1 < g->numberOfVertices; vi1++) {
        for(int vi2 = vi1+1; vi2 < g->numberOfVertices; vi2++) {
            potDist1 = dist(g,vi1,v1) + 1 + dist(g,vi2,v2);
            potDist2 = dist(g,vi1,v2) + 1 + dist(g,vi2,v1);
            minPotDist = MIN(potDist1,potDist2);
            if(g->distMatrix[vi1][vi2] > minPotDist) {
                // if previous distance was larger than or equal to girth-1
                // and now it's smaller than girth-1, then we can decrement
                // the amount of legal choices for both vertices
                if(g->distMatrix[vi1][vi2] >= girth - 1 && minPotDist < girth - 1) {
                    removeElement(g->legalChoices[vi1], vi2);
                    removeElement(g->legalChoices[vi2], vi1);
                }
                g->distMatrix[vi1][vi2] = minPotDist;
            }
        }
    }

    if (deg(g,v1) == g->maxDeg)
        makeVertexIllegal(g, v1);
    if (deg(g,v2) == g->maxDeg)
        makeVertexIllegal(g, v2);

    if (g->distMatrix[0][v1] < minDistMaxDeg && v1 != 0) {
        removeElement(g->indexMaxDeg, v1);
        makeVertexIllegal(g, v1);
    }

    if (g->distMatrix[0][v2] < minDistMaxDeg && v2 != 0) {
        removeElement(g->indexMaxDeg, v2);
        makeVertexIllegal(g, v2);
    }
}

// Returns false if adding an edge between v1 and v2 led to having two vertices of maxDeg
// being at a distance smaller than minDistMaxDeg. Otherwise, returns false.
bool minDistMaxDegSatisfied(struct graph *g, int minDistMaxDeg, int v1, int v2, int searchDepth,
                            bitset *legalChoicesChanged) {

    // Only check if minDistMaxDeg > 1. If minDistMaxDeg == 1, checking is redundant
    if(minDistMaxDeg <= 1)
        return true;

    // based on https://stackoverflow.com/questions/10258305/how-to-implement-a-breadth-first-search-to-a-certain-depth
    if (deg(g, v1) == g->minDeg + 1) {
        int dist = 0;
        int amtToDoToDepthIncrease = 1;
        int counterToDepthIncrease = 0;
        startQueue = 0;
        endQueue = 0;
        bfsQueue[0] = v1;
        while (startQueue <= endQueue) {
            int currVertex = bfsQueue[startQueue];
            startQueue++;
            if (currVertex != v1) {
                if (deg(g, currVertex) > g->minDeg)
                    return false;
                else if (deg(g, currVertex) == g->minDeg &&
                         !makeVertexIllegal2AddEdge(g, currVertex, searchDepth, legalChoicesChanged))
                    return false;
            }
            forEach(neigh, g->adjacencyList[currVertex]) {
                if (dist(g, neigh, v1) > dist) {
                    bfsQueue[++endQueue] = neigh;
                    counterToDepthIncrease++;
                }
            }
            if (--amtToDoToDepthIncrease == 0) {
                if (++dist >= minDistMaxDeg) break;
                amtToDoToDepthIncrease = counterToDepthIncrease;
                counterToDepthIncrease = 0;
            }
        }
    } else if (deg(g, v1) == g->minDeg) {
        int dist = 0;
        int amtToDoToDepthIncrease = 1;
        int counterToDepthIncrease = 0;
        startQueue = 0;
        endQueue = 0;
        bfsQueue[0] = v1;
        while (startQueue <= endQueue) {
            int currVertex = bfsQueue[startQueue];
            startQueue++;
            if (deg(g, currVertex) > g->minDeg && currVertex != v1 &&
                !makeVertexIllegal2AddEdge(g, v1, searchDepth, legalChoicesChanged))
                return false;
            forEach(neigh, g->adjacencyList[currVertex]) {
                if (dist(g, neigh, v1) > dist) {
                    bfsQueue[++endQueue] = neigh;
                    counterToDepthIncrease++;
                }
            }
            if (--amtToDoToDepthIncrease == 0) {
                if (++dist >= minDistMaxDeg) break;
                amtToDoToDepthIncrease = counterToDepthIncrease;
                counterToDepthIncrease = 0;
            }
        }
    }

    if (deg(g, v2) == g->minDeg + 1) {
        int dist = 0;
        int amtToDoToDepthIncrease = 1;
        int counterToDepthIncrease = 0;
        startQueue = 0;
        endQueue = 0;
        bfsQueue[0] = v2;
        while (startQueue <= endQueue) {
            int currVertex = bfsQueue[startQueue];
            startQueue++;
            if (currVertex != v2) {
                if (deg(g, currVertex) > g->minDeg)
                    return false;
                else if (deg(g, currVertex) == g->minDeg &&
                         !makeVertexIllegal2AddEdge(g, currVertex, searchDepth, legalChoicesChanged))
                    return false;
            }
            forEach(neigh, g->adjacencyList[currVertex]) {
                if (dist(g, neigh, v2) > dist) {
                    bfsQueue[++endQueue] = neigh;
                    counterToDepthIncrease++;
                }
            }
            if (--amtToDoToDepthIncrease == 0) {
                if (++dist >= minDistMaxDeg) break;
                amtToDoToDepthIncrease = counterToDepthIncrease;
                counterToDepthIncrease = 0;
            }
        }
    } else if (deg(g, v2) == g->minDeg) {
        int dist = 0;
        int amtToDoToDepthIncrease = 1;
        int counterToDepthIncrease = 0;
        startQueue = 0;
        endQueue = 0;
        bfsQueue[0] = v2;
        while (startQueue <= endQueue) {
            int currVertex = bfsQueue[startQueue];
            startQueue++;
            if (deg(g, currVertex) > g->minDeg && currVertex != v2 &&
                !makeVertexIllegal2AddEdge(g, v2, searchDepth, legalChoicesChanged))
                return false;
            forEach(neigh, g->adjacencyList[currVertex]) {
                if (dist(g, neigh, v2) > dist) {
                    bfsQueue[++endQueue] = neigh;
                    counterToDepthIncrease++;
                }
            }
            if (--amtToDoToDepthIncrease == 0) {
                if (++dist >= minDistMaxDeg) break;
                amtToDoToDepthIncrease = counterToDepthIncrease;
                counterToDepthIncrease = 0;
            }
        }
    }
    return true;
}

// If the level l of vertex v contains comb[l] vertices of degree larger than g->minDeg,
// then the vertices of degree minDeg at level l are made illegal.
// Returns false if making this illegal leads to not having enough choices: oterhwise retursn true.
bool checkLevelComb(struct graph *g, int v, int searchDepth, bitset *legalChoicesChanged, int comb[]) {
    int levelV = g->levelOfVertex[v];
    if (levelV < g->heightMooreTree && levelV > 0 &&
        size(difference(g->verticesAtLevel[levelV], g->verticesOfDeg[g->minDeg])) == comb[levelV]) {
        forEach(vSameLevel, intersection(g->verticesAtLevel[levelV],g->verticesOfDeg[g->minDeg])) {
            if (vSameLevel != v && !makeVertexIllegal2AddEdge(g, vSameLevel, searchDepth, legalChoicesChanged)) {
                return false;
            }
        }
    }
    return true;
}


// Pre-condition: deg(v1) < g->maxDeg and deg(v2) < g->maxDeg
// Returns true if adding the edge was successful;
// returns false if adding the edge led to a vertex having not enough
// legalChoices left to complete to its full degree.
// This function is used when recursively adding edges to the graph after the Moore tree is constructed.
// Strongly based on EGR
// (https://github.com/JorikJooken/edgeGirthRegularGraphs/blob/master/Code/generateRGLambdaGraphs.c)
bool addEdgeWithDistMatrixLegal2(struct graph *g, int girth, int v1, int v2, int searchDepth,
                                 bitset *legalChoicesChanged, int* numDistChanged, int comb[]) {
    addEdge(g, v1, v2);

    if (deg(g, v1) == g->minDeg+1 && !checkLevelComb(g, v1, searchDepth, legalChoicesChanged, comb))
        return false;
    if (deg(g, v2) == g->minDeg+1 && !checkLevelComb(g, v2, searchDepth, legalChoicesChanged, comb))
        return false;


    if (deg(g, v1) == g->maxDeg)
        if (!makeVertexIllegal2AddEdge(g, v1, searchDepth, legalChoicesChanged))
            return false;

    if (deg(g, v2) == g->maxDeg)
        if (!makeVertexIllegal2AddEdge(g, v2, searchDepth, legalChoicesChanged))
            return false;

    *numDistChanged = -1;

    // Set distBFS to numVertices for all vertices
    for (int i = 0; i < g->numberOfVertices; i++)
        distBFS1[i] = distBFS2[i] = g->numberOfVertices;

    // BFS for v1
    startQueue = 0;
    endQueue = 0;
    bfsQueue[0] = v1;
    distBFS1[v1] = 0;
    while (startQueue <= endQueue) {
        int currVertex = bfsQueue[startQueue];
        startQueue++;
        forEach(neigh, g->adjacencyList[currVertex]) {
            if (distBFS1[neigh] == g->numberOfVertices) {
                distBFS1[neigh] = distBFS1[currVertex] + 1;
                bfsQueue[++endQueue] = neigh;
            }
        }
    }

    // BFS for v2
    startQueue = 0;
    endQueue = 0;
    bfsQueue[0] = v2;
    distBFS2[v2] = 0;
    while (startQueue <= endQueue) {
        int currVertex = bfsQueue[startQueue];
        startQueue++;
        forEach(neigh, g->adjacencyList[currVertex]) {
            if (distBFS2[neigh] == g->numberOfVertices) {
                distBFS2[neigh] = distBFS2[currVertex] + 1;
                bfsQueue[++endQueue] = neigh;
            }
        }
    }

    // update distances and keep checking if enough choices left
    int newDist;
    forEach(v, g->indexMaxDeg) {
        forEachAfterIndex(neigh, g->legalChoices[v], v) {
            newDist = g->distMatrix[v][neigh];
            if (1 + distBFS1[v] + distBFS2[neigh] < newDist)
                newDist = 1 + distBFS1[v] + distBFS2[neigh];
            if (1 + distBFS2[v] + distBFS1[neigh] < newDist)
                newDist = 1 + distBFS2[v] + distBFS1[neigh];

            if (newDist < g->distMatrix[v][neigh]) {
                if (g->distMatrix[v][neigh] >= girth - 1 && newDist < girth - 1) {
                    if (!contains(*legalChoicesChanged, v)) {
                        add(*legalChoicesChanged, v);
                        oldLegalChoicesAddEdge[searchDepth][v] = g->legalChoices[v];
                    }
                    removeElement(g->legalChoices[v], neigh);
                    if (!enoughChoicesLeftVertex(g, v))
                        return false;

                    if (!contains(*legalChoicesChanged, neigh)) {
                        add(*legalChoicesChanged, neigh);
                        oldLegalChoicesAddEdge[searchDepth][neigh] = g->legalChoices[neigh];
                    }
                    removeElement(g->legalChoices[neigh], v);
                    if (!enoughChoicesLeftVertex(g, neigh))
                        return false;
                }
                *numDistChanged = *numDistChanged + 1;
                endpoint1[searchDepth][*numDistChanged] = v;
                endpoint2[searchDepth][*numDistChanged] = neigh;
                oldDist[searchDepth][v][neigh] = g->distMatrix[v][neigh];
                g->distMatrix[v][neigh] = newDist;
            }
        }
    }

    return true;
}

// Resets the legalChoices of the vertices in legalChoicesChanged to these in oldLegalChoicesAddEdge.
// Also resets the distance matrix of the vertices in legalChoicesChanged and removes the edge v1,v2 in graph g.
void removeEdgeWithDistMatrixAndReset(struct graph *g, int v1, int v2, int searchDepth, bitset *legalChoicesChanged,
                                      int *numDistChanged) {
    // Rest legal choices
    forEach(v, *legalChoicesChanged) {
        g->legalChoices[v] = oldLegalChoicesAddEdge[searchDepth][v];
    }

    // Reset distance matrix
    int u,v;
    for(int i = 0; i <= *numDistChanged; i++) {
        u = endpoint1[searchDepth][i];
        v = endpoint2[searchDepth][i];
        g->distMatrix[u][v] = oldDist[searchDepth][u][v];
    }

    // Remove edge
    removeEdge(g, v1, v2);
}

// Resets the legalChoices of the vertices in legalChoicesChanged to these in oldLegalChoices.
void resetLegalChoices(struct graph *g, int searchDepth, bitset *legalChoicesChanged) {
    forEach(v, *legalChoicesChanged) {
        g->legalChoices[v] = oldLegalChoices[searchDepth][v];
    }
}




/************************************************
 * Initializing and deleting graph
 ***********************************************/

// Returns a graph of numVertices
struct graph makeEmptyGraph(int numVertices, int minDegree, int maxDegree, int h) {
    // From GenK2
    struct graph g = {.numberOfVertices = numVertices};

    g.heightMooreTree = h;

    // Allocate space for graph
    g.adjacencyList = malloc(sizeof(bitset)*numVertices);
    g.verticesOfDeg = malloc(sizeof(bitset)*(maxDegree+1));
    g.legalChoices = malloc(sizeof(bitset)*numVertices);
    g.verticesAtLevel = malloc(sizeof(bitset)*h);

    // Initialize distance matrix
    for(int i = 0; i < numVertices; i++) {
        g.distMatrix[i][i] = 0;
        // Initialize values to numVertices, which means they are not connected
        for (int j = i+1; j < numVertices; j++) {
                g.distMatrix[i][j] = numVertices;
        }
    }

    // Initialize the rest to obtain a valid empty graph
    emptyGraph(&g, minDegree, maxDegree);

    for (int i = 0; i < h; i++) {
        g.verticesAtLevel[i] = EMPTY;
    }

    for (int i = 0; i < numVertices; i++) {
        g.levelOfVertex[i] = numVertices;
    }

    // Initialize every vertex as being legal to connect to every other vertex
    for(int i = 0; i < numVertices; i++)
        g.legalChoices[i] = compl(singleton(i), numVertices);

    return g;
}

// Frees allocated memory of the graph g
void deleteGraph(struct graph *g) {
    free(g->adjacencyList);
    free(g->verticesOfDeg);
    free(g->legalChoices);
    free(g->verticesAtLevel);
}
