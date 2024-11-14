#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../nauty/nauty.h"
#include "../nauty/gtools.h"
#include "../nauty/gutils.h"
#include "bitset/bitset.h"


/************************************************
 *  - Regular to biregular graph construction algorithm -
 ***********************************************/


// Based on GenK2
#define USAGE \
"\nUsage:./regToBiregExec r m girth minN maxN\n\n"


#define MIN(a,b) (((a)<(b))?(a):(b))
#define EVEN(n) (n % 2 == 0)
#define ODD(n) (n % 2 == 1)

struct edge { int v1,v2; };


/************************************************
 *  Global variables
 ***********************************************/

struct edge *edges;

int **distEdgeToVertex;
bitset_t **possibleEdges;

int leastGirthWanted;
bool foundBiregGraph = false;
int maxm;

//  Macro's for nauty representation
#define FOREACH(element,nautySet,maxm)\
 for(int element = nextelement((nautySet),maxm,-1); (element) >= 0;\
 (element) = nextelement((nautySet),maxm,(element)))
#define REMOVEONEEDGE(g,i,j,maxm)\
 DELELEMENT(GRAPHROW(g,i,maxm),j); DELELEMENT(GRAPHROW(g,j,maxm),i)



/************************************************
*  Helper functions related to counting vertices
***********************************************/

// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
int ipow(int base, int exp) {
    int result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

// Pre-condition: base != 1
int powerSum(int base, int max) {
    return (ipow(base, max+1) - 1) / (base - 1);
}


int getNumVerticesLevel0(int g) {
    return (ODD(g) ? 1 : 2);
}

int getNumVerticesLevel1(int r, int m, int g, int config[]) {
    if(ODD(g))
        return m;
    return m + (config[0] == 1 ? r : m) - 2;
}

int getNumVerticesComb(int r, int m, int g, int t, int config[]) {
    int sum = getNumVerticesLevel0(g);
    sum += getNumVerticesLevel1(r, m, g, config) * powerSum(r-1, t-1);
    for(int level = 1; level <= t-1; level++) {
        // Remember: for even g=2d: t=d-1
        sum += (m-r) * config[level] * powerSum(r-1, t-1-level);
    }
    return sum;
}

// Sets the edges global variables to the edges of graph g with order n.
// Returns the number of edges in g.
int setEdges(graph *g, int n) {
    int count = 0;
    for(int v = 0; v < n; v++) {
        FOREACH(nbr, GRAPHROW(g, v, maxm),maxm) {
            if (nbr > v) {
                struct edge newEdge = { .v1 = v, .v2 = nbr };
                edges[count++] = newEdge;
            }
        }
    }
    return count;
}



/************************************************
*  Printing functions
***********************************************/

void printDist(int n, int nOfEdges) {
    for(int i = 0; i < nOfEdges; i++) {
        for(int j = 0; j < n; j++) {
            fprintf(stderr, "%d ", distEdgeToVertex[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void bitset_myprint(bitset_t* b) {
    fprintf(stderr,"{");
    for(size_t i = 0; bitset_next_set_bit(b,&i) ; i++) {
        fprintf(stderr, "%zu, ",i);
    }
    fprintf(stderr,"}\n");
}



/************************************************
*  Construction functions
***********************************************/

bool noVertexInCommon(int e1, int e2) {
    return edges[e1].v1 != edges[e2].v1 && edges[e1].v1 != edges[e2].v2 &&
           edges[e1].v2 != edges[e2].v1 && edges[e1].v2 != edges[e2].v2;
}

int distTwoEdges(int e1, int e2) {
    return MIN(distEdgeToVertex[e1][edges[e2].v1], distEdgeToVertex[e1][edges[e2].v2]);
}

bool edgesTooClose(int e1, int e2, int minDist) {
    return distTwoEdges(e1, e2) < minDist;
}


// Removes the chosenEdges in the graph and adds new edges to the new vertex newV1,
// which will be of degree m. If it succeeds to generate a graph with girth >= leastGirthWanted,
// then it sets foundBiregGraph tot true.
void tryAddingEdgesEvenM(graph *g, int* chosenEdges, int numEdgesToAdd, int newV1, int newNumVertices) {
    int newGirth;

    for (int i = 0; i < numEdgesToAdd; i++) {
        REMOVEONEEDGE(g, edges[chosenEdges[i]].v1, edges[chosenEdges[i]].v2, maxm);
        ADDONEEDGE(g, newV1, edges[chosenEdges[i]].v1, maxm);
        ADDONEEDGE(g, newV1, edges[chosenEdges[i]].v2, maxm);
    }

    newGirth = girth(g, maxm, newNumVertices);

    if(newGirth >= leastGirthWanted) {
        foundBiregGraph = true;
        printf("n, new girth: %d, %d\n", newNumVertices, newGirth);
        char *g6String = ntos6(g,maxm,newNumVertices);
        printf("%s", g6String);
    }

    for (int i = 0; i < numEdgesToAdd; i++) {
        ADDONEEDGE(g, edges[chosenEdges[i]].v1, edges[chosenEdges[i]].v2, maxm);
        REMOVEONEEDGE(g, newV1, edges[chosenEdges[i]].v1, maxm);
        REMOVEONEEDGE(g, newV1, edges[chosenEdges[i]].v2, maxm);
    }
}

// Removes the chosenEdges in the graph and adds new edges to the new vertices newV1 and newV2,
// which will both be of degree m. If it succeeds to generate a graph with girth >= leastGirthWanted,
// then it sets foundBiregGraph tot true.
void tryAddingEdgesOddM(graph *g, int* chosenEdges, int numEdgesToAdd, int newV1, int newV2, int newNumVertices) {
    int newGirth;

    ADDONEEDGE(g, newV1, newV2, maxm);

    for (int i = 0; i < numEdgesToAdd; i++) {
        REMOVEONEEDGE(g, edges[chosenEdges[i]].v1, edges[chosenEdges[i]].v2, maxm);
    }

    // Finding all subsets based on ChatGPT
    int numSubsets = 1 << numEdgesToAdd; // Number of subsets (2^n)
    for (int subsetMask = 0; subsetMask < numSubsets && !foundBiregGraph; subsetMask++) {
        for (int j = 0; j < numEdgesToAdd; j++) {
            if ((subsetMask & (1 << j)) == 0) {
                ADDONEEDGE(g, newV1, edges[chosenEdges[j]].v1, maxm);
                ADDONEEDGE(g, newV2, edges[chosenEdges[j]].v2, maxm);
            } else {
                ADDONEEDGE(g, newV1, edges[chosenEdges[j]].v2, maxm);
                ADDONEEDGE(g, newV2, edges[chosenEdges[j]].v1, maxm);
            }
        }

        newGirth = girth(g, maxm, newNumVertices);

        if (newGirth >= leastGirthWanted) {
            foundBiregGraph = true;
            printf("n, new girth: %d, %d\n", newNumVertices, newGirth);
            char *g6String = ntos6(g,maxm,newNumVertices);
            printf("%s", g6String);
        }

        for (int j = 0; j < numEdgesToAdd; j++) {
            if ((subsetMask & (1 << j)) == 0) {
                REMOVEONEEDGE(g, newV1, edges[chosenEdges[j]].v1, maxm);
                REMOVEONEEDGE(g, newV2, edges[chosenEdges[j]].v2, maxm);
            } else {
                REMOVEONEEDGE(g, newV1, edges[chosenEdges[j]].v2, maxm);
                REMOVEONEEDGE(g, newV2, edges[chosenEdges[j]].v1, maxm);
            }
        }
    }

    for (int i = 0; i < numEdgesToAdd; i++) {
        ADDONEEDGE(g, edges[chosenEdges[i]].v1, edges[chosenEdges[i]].v2, maxm);
    }

    REMOVEONEEDGE(g, newV1, newV2, maxm);
}



// Try the construction
void makeBiregular(graph *g, int n, int r, int m) {
    if (r >= m)
        return;

    int nOfEdges = setEdges(g,n);
    int minDist = (ODD(m)) ? leastGirthWanted-3 : leastGirthWanted-2;
    int numEdgesToAdd = (ODD(m)) ? m-1 : m/2;
    int newV1 = n;
    int newV2 = n+1;
    int newNumVertices = (ODD(m)) ? n+2 : n+1;
    int chosenEdges[numEdgesToAdd];
    bool valid = true;

    bitset_t* enoughPossibleEdges = bitset_create_with_capacity(nOfEdges);

    // Calculate the distances between the edges
    for (int e1 = 0; e1 < nOfEdges; e1++) {
        find_dist2(g, maxm, n, edges[e1].v1, edges[e1].v2, distEdgeToVertex[e1]);
        possibleEdges[e1] = bitset_create_with_capacity(nOfEdges);
        int ctr = 0;
        for (int e2 = 0; e2 < nOfEdges; e2++) {
            if (!edgesTooClose(e1,e2,minDist)) {
                bitset_set(possibleEdges[e1], e2);
                ctr++;
            }
        }
        if (ctr >= numEdgesToAdd-1) {
            bitset_set(enoughPossibleEdges, e1);
        }
    }

    // Calculate the set of edges that have enough edges at minDist or further.
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t e1 = 0; bitset_next_set_bit(enoughPossibleEdges, &e1); e1++) {
            if (!bitset_contains_all(enoughPossibleEdges, possibleEdges[e1])) {
                bitset_inplace_intersection(possibleEdges[e1], enoughPossibleEdges);
                changed = true;
            }
            if (bitset_count(possibleEdges[e1]) < numEdgesToAdd-1)
                bitset_set_to_value(enoughPossibleEdges,  e1, false);
        }
    }

    if (bitset_count(enoughPossibleEdges) < numEdgesToAdd) {
        for (int e = 0; e < nOfEdges; e++) {
            bitset_free(possibleEdges[e]);
        }
        bitset_free(enoughPossibleEdges);
        return;
    }

    // Get first possible edge
    size_t firstEdge = 0;
    bitset_next_set_bit(enoughPossibleEdges, &firstEdge);
    chosenEdges[0] = firstEdge;

    int numChosenEdges = 1;

    bitset_t* intersecPossibleEdges[numEdgesToAdd];
    for (int i = 1; i < numEdgesToAdd; i++) {
        intersecPossibleEdges[i] = bitset_create_with_capacity(nOfEdges);
    }
    intersecPossibleEdges[0] = bitset_copy(possibleEdges[firstEdge]);


    while(!foundBiregGraph) {
        valid = true;

        if (numChosenEdges == numEdgesToAdd) {
            if (ODD(m))
                tryAddingEdgesOddM(g, chosenEdges, numEdgesToAdd, newV1, newV2, newNumVertices);
            else
                tryAddingEdgesEvenM(g, chosenEdges, numEdgesToAdd, newV1, newNumVertices);
        }

        if (bitset_empty(intersecPossibleEdges[numChosenEdges - 1]))
            valid = false;

        // Logic to iterate over the possible edges
        bool goBack = true;
        if (!foundBiregGraph) {
            // Change an edge
            size_t nextEdge;
            if (valid && numChosenEdges < numEdgesToAdd) {
                nextEdge = chosenEdges[numChosenEdges - 1]+1;
                if (bitset_next_set_bit(intersecPossibleEdges[numChosenEdges - 1], &nextEdge)) {
                    chosenEdges[numChosenEdges] = nextEdge;
                    bitset_intersection(intersecPossibleEdges[numChosenEdges - 1], possibleEdges[nextEdge],
                                        intersecPossibleEdges[numChosenEdges]);
                    numChosenEdges++;
                    goBack = false;
                }
            }
            if (goBack) {
                bitset_t* bsetToIterate;
                do {
                    numChosenEdges--;
                    if (numChosenEdges < 0) {
                        for (int e = 0; e < nOfEdges; e++) {
                            bitset_free(possibleEdges[e]);
                        }
                        for (int e = 0; e < numEdgesToAdd; e++) {
                            bitset_free(intersecPossibleEdges[e]);
                        }
                        bitset_free(enoughPossibleEdges);
                        return;
                    } else if (numChosenEdges == 0) {
                        bsetToIterate = enoughPossibleEdges;
                    } else {
                        bsetToIterate = intersecPossibleEdges[numChosenEdges - 1];
                    }
                    nextEdge = chosenEdges[numChosenEdges]+1;
                } while (!bitset_next_set_bit(bsetToIterate, &nextEdge));
                chosenEdges[numChosenEdges] = nextEdge;
                bitset_intersection(bsetToIterate, possibleEdges[nextEdge],
                                    intersecPossibleEdges[numChosenEdges]);
                numChosenEdges++;
            }
        } else {
            for (int e = 0; e < nOfEdges; e++) {
                bitset_free(possibleEdges[e]);
            }
            for (int e = 0; e < numEdgesToAdd; e++) {
                bitset_free(intersecPossibleEdges[e]);
            }
            bitset_free(enoughPossibleEdges);
            return;
        }
    }
}



int main(int argc, char ** argv) {
    if (argc != 6) {
        // Print the usage information and exit with an error code
        printf(USAGE);
        return 1;
    }

    int r = atoi(argv[1]);
    int m = atoi(argv[2]);
    leastGirthWanted = atoi(argv[3]);
    int minVertices = atoi(argv[4]);
    int maxVertices = atoi(argv[5]);

    int n,initialGirth;
    int radius, diameter;
    int lowerBound, upperBound;

    lowerBound = m % 2 ? minVertices-2 : minVertices-1;
    upperBound = m % 2 ? maxVertices-3 : maxVertices-2;

    //  Start looping over lines of stdin.
    char *graphString = NULL;
    size_t size;
    int ctrUnfiltered = 0;
    int ctrAllGraphs = 0;

    int maxn = maxVertices;
    int maxEdges = (maxn*r >> 1) + m;

    // Do mallocs
    distEdgeToVertex = (int **)malloc(maxEdges * sizeof(int *));
    if (distEdgeToVertex == NULL) {
        printf("Memory allocation dist failed\n");
        return -1;
    }

    for (int i = 0; i < maxEdges; i++) {
        distEdgeToVertex[i] = (int *)malloc(maxn * sizeof(int));
        if (distEdgeToVertex[i] == NULL) {
            printf("Memory allocation failed for row %d!\n", i);
            return -1;
        }
    }

    edges = (struct edge*)malloc(maxEdges * sizeof(struct edge));
    if (edges == NULL) {
        printf("Memory allocation edges failed\n");
        return -1;
    }

    possibleEdges = (bitset_t**)malloc(maxEdges*sizeof(bitset_t*));
    if (possibleEdges == NULL) {
        printf("Memory allocation possibleEdges failed\n");
        return -1;
    }

    maxm = SETWORDSNEEDED(maxn);

    DYNALLSTAT(graph,g,g_sz);
    DYNALLOC2(graph,g,g_sz,maxm,maxn,"malloc");


    while(getline(&graphString, &size, stdin) != -1) {
        ctrAllGraphs++;
        n = graphsize(graphString);

        if (n > upperBound)
            break;

        if (EVEN(leastGirthWanted) && leastGirthWanted >= 8 && EVEN(m)) {
            int t = leastGirthWanted/2;
            int potentialNewLowerBound = m * ceil((double) ipow(r-1,t-1)-1 / (double) (r-2) +
                                                  (double) ipow(r-1,t-1) / (double) r);
            if (potentialNewLowerBound > lowerBound)
                lowerBound = potentialNewLowerBound;
        }

        if (ODD(m)) {
            int t, potentialNewLowerBound;
            if (EVEN(leastGirthWanted)) {
                t = leastGirthWanted/2-1;
                int *config = (int *)malloc(t * sizeof(int));
                config[0] = 2 ;
                for (int i = 1; i < t; i++) {
                    config[i] = 0;
                }
                potentialNewLowerBound = getNumVerticesComb(r, m, leastGirthWanted, t, config) - 2;
                free(config);
            } else {
                t = leastGirthWanted/2;
                int *config = (int *)malloc(t * sizeof(int));
                config[0] = 1;
                config[1] = 1;
                for (int i = 2; i < t; i++) {
                    config[i] = 0;
                }
                potentialNewLowerBound = getNumVerticesComb(r, m, leastGirthWanted, t, config) - 2;
                free(config);
            }
            if (potentialNewLowerBound > lowerBound)
                lowerBound = potentialNewLowerBound;
        }

        if (n >= lowerBound && n <= upperBound) {

            stringtograph(graphString,g,maxm);
            initialGirth = girth(g, maxm, n);
            if (initialGirth >= leastGirthWanted) {

                // Compute diameter
                diamstats(g, maxm, n, &radius, &diameter);
                int diameterBound = diameter+2;
                if (ODD(m))
                    diameterBound++;

                if (leastGirthWanted <= diameterBound) {
                    ctrUnfiltered++;

                    printf("diameter=%d\n", diameter);
                    printf("initial girth=%d\n", initialGirth);
                    printf("initial n=%d\n", n);

                    makeBiregular(g, n, r, m);

                    if (foundBiregGraph)
                        break;
                }
            }
        }
    }

    printf("filtered %d graphs out of %d\n", ctrAllGraphs-ctrUnfiltered, ctrAllGraphs);

    free(graphString);

    for (int i = 0; i < maxEdges; i++) {
        free(distEdgeToVertex[i]);
    }
    free(distEdgeToVertex);
    free(edges);
    free(possibleEdges);

    return 0;
}






