#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../nauty/nauty.h"
#include "../nauty/gtools.h"
#include "../nauty/gutils.h"


#define MIN(a,b) (((a)<(b))?(a):(b))

struct edge { int v1,v2; };


#define USAGE \
"\nUsage:./calcMaxAtMinDistExec [-v|-e|-g] \n\n"

#define HELPTEXT \
"Computes the maximum amount of vertices and/or edges at minimum pairwise distance ceil(g/2).\n\
of an (r,g)-graph. It reads these (r,g)-graphs in graph6 or sparse6 format from stdin.\n\
All output is sent to stderr. This implementation uses parallelisation using openmp, so \n\
make sure it is installed.\n\
\n\
There are three optional arguments.\n\
    -v,     output the most amount of vertices at the minimum pairwise distance\n\
            found so far.\n\
    -e,     output the most amount of edges at the minimum pairwise distance\n\
            found so far.\n\
    -g,     applies a greedy search for finding a number of vertices and edges\n\
            at the minimum pairwise distance. This number is not necessarily optimal.\n\
When no argument is given, it computes the maximum amount of vertices and edges\n\
exactly and only shows the the output when the computation is finished.\n"


/************************************************
 *  Global variables
 ***********************************************/


int maxm;

int numEdges;
struct edge *edges;

int **dist;

int **currVerticesAll;
int **currEdgesAll;

bool doApprox = false;
int bestMaxVerticesSoFar = 0;
int bestMaxEdgesSoFar = 0;


//  Macro's for nauty representation
#define FOREACH(element,nautySet,maxm)\
for(int element = nextelement((nautySet),maxm,-1); (element) >= 0;\
(element) = nextelement((nautySet),maxm,(element)))



bool allocateArrays(int n, int r){
    dist = (int **)malloc(n * sizeof(int *));

    if (dist == NULL) {
        fprintf(stderr,"Memory allocation dist failed\n");
        return false;
    }

    for (int i = 0; i < n; i++) {
        dist[i] = (int *)malloc(n * sizeof(int));
        if (dist[i] == NULL) {
            fprintf(stderr,"Memory allocation failed for row %d!\n", i);
            return false;
        }
    }

    currVerticesAll = (int **)malloc(n * sizeof(int *));

    if (currVerticesAll == NULL) {
        fprintf(stderr,"Memory allocation currVerticesAll failed\n");
        return false;
    }

    for (int i = 0; i < n; i++) {
        currVerticesAll[i] = (int *)malloc(n * sizeof(int));
        if (currVerticesAll[i] == NULL) {
            fprintf(stderr,"Memory allocation failed for row %d!\n", i);
            return false;
        }
    }

    numEdges = (n*r) >> 1;

    edges = (struct edge*)malloc(numEdges * sizeof(struct edge));

    if (edges == NULL) {
        fprintf(stderr,"Memory allocation edges failed\n");
        return false;
    }

    currEdgesAll = (int **)malloc(numEdges * sizeof(int *));

    if (currEdgesAll == NULL) {
        fprintf(stderr,"Memory allocation currEdgesAll failed\n");
        return false;
    }

    for (int i = 0; i < numEdges; i++) {
        currEdgesAll[i] = (int *)malloc(n * sizeof(int));
        if (currEdgesAll[i] == NULL) {
            fprintf(stderr,"Memory allocation failed for row %d!\n", i);
            return false;
        }
    }

    return true;
}


void freeArrays(int n){
    for (int i = 0; i < n; i++)
        free(dist[i]);
    free(dist);

    for (int i = 0; i < n; i++)
        free(currVerticesAll[i]);
    free(currVerticesAll);

    for (int i = 0; i < numEdges; i++)
        free(currEdgesAll[i]);
    free(currEdgesAll);

    free(edges);
}


// Sets the edges global variables to the edges of graph g with order n.
// Returns the number of edges in g.
void setEdges(graph *g, int n) {
    int count = 0;
    for(int v = 0; v < n; v++) {
        FOREACH(nbr, GRAPHROW(g, v, maxm),maxm) {
            if (nbr > v) {
                struct edge newEdge = { .v1 = v, .v2 = nbr };
                edges[count++] = newEdge;
            }
        }
    }
}

int distTwoEdges(int e1, int e2) {
    return MIN(
            MIN(dist[edges[e1].v1][edges[e2].v1], dist[edges[e1].v1][edges[e2].v2]),
            MIN(dist[edges[e1].v2][edges[e2].v1], dist[edges[e1].v2][edges[e2].v2])
            );
}



// Try the construction
int checkDistances(graph *g, int n, int minDist, int startVertex, int ctrVertices) {

    int maxAtMinDist = ctrVertices;

    for (int v = currVerticesAll[startVertex][ctrVertices-1]+1; v < n; v++) {

        bool distantEnough = true;

        for (int u = 0; u < ctrVertices; u++) {
            if (dist[currVerticesAll[startVertex][u]][v] < minDist) {
                distantEnough = false;
                break;
            }
        }

        if (distantEnough) {
            currVerticesAll[startVertex][ctrVertices] = v;
            ctrVertices += 1;
            int numAtMinDist = checkDistances(g, n, minDist, startVertex, ctrVertices);
            if (maxAtMinDist < numAtMinDist) {
                maxAtMinDist = numAtMinDist;

                if (doApprox) {
                    # pragma omp critical
                    {
                        if (bestMaxVerticesSoFar < maxAtMinDist) {
                            bestMaxVerticesSoFar = maxAtMinDist;
                            fprintf(stderr,"Best maxVerticesSoFar = %d\n", bestMaxVerticesSoFar);
                        }
                    }
                }
            }
            ctrVertices -= 1;
        }
    }

    return maxAtMinDist;
}


int checkDistancesEdges(graph *g, int numEdges, int minDist, int startEdge, int ctrEdges) {
    int maxAtMinDist = ctrEdges;
    for (int e1 = currEdgesAll[startEdge][ctrEdges-1]+1; e1 < numEdges; e1++) {
        bool distantEnough = true;
        for (int e2 = 0; e2 < ctrEdges; e2++) {
            if (distTwoEdges(currEdgesAll[startEdge][e2],e1) < minDist) {
                distantEnough = false;
                break;
            }
        }
        if (distantEnough) {
            currEdgesAll[startEdge][ctrEdges] = e1;
            ctrEdges += 1;
            int numAtMinDist = checkDistancesEdges(g, numEdges, minDist, startEdge, ctrEdges);
            if (maxAtMinDist < numAtMinDist)
                maxAtMinDist = numAtMinDist;

            if (doApprox) {
                # pragma omp critical
                {
                    if (bestMaxEdgesSoFar < maxAtMinDist) {
                        bestMaxEdgesSoFar = maxAtMinDist;
                        fprintf(stderr,"Best maxEdgesSoFar = %d\n", bestMaxEdgesSoFar);
                    }
                }
            }
            ctrEdges -= 1;
        }
    }
    return maxAtMinDist;
}


// Try the construction greedily
int checkDistancesGreedy(graph *g, int n, int minDist, int startVertex) {
    int ctrVertices = 1;
    for (int v = startVertex+1; v < n; v++) {
        bool distantEnough = true;
        for (int u = 0; u < ctrVertices; u++) {
            if (dist[currVerticesAll[startVertex][u]][v] < minDist) {
                distantEnough = false;
                break;
            }
        }
        if (distantEnough) {
            currVerticesAll[startVertex][ctrVertices] = v;
            ctrVertices += 1;
        }
    }
    return ctrVertices;
}


// Try the construction greedily
int checkDistancesEdgesGreedy(graph *g, int numEdges, int minDist, int startEdge) {
    int ctrEdges = 1;
    for (int e1 = startEdge+1; e1 < numEdges; e1++) {
        bool distantEnough = true;
        for (int e2 = 0; e2 < ctrEdges; e2++) {
            if (distTwoEdges(currEdgesAll[startEdge][e2],e1) < minDist) {
                distantEnough = false;
                break;
            }
        }
        if (distantEnough) {
            currEdgesAll[startEdge][ctrEdges] = e1;
            ctrEdges += 1;
        }
    }
    return ctrEdges;
}


int main(int argc, char ** argv) {
    int n;
    int radius, diameter, theGirth;
    unsigned long es;
    int r, mincount, maxdeg, maxcount;
    boolean eulerian;

    bool doVertices;
    bool doGreedy = false;

    // Loop through command-line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            fprintf(stderr, "%s", USAGE);
            fprintf(stderr, "%s", HELPTEXT);
            return 0;
        }
        if (strcmp(argv[i], "-v") == 0) {
            doApprox = true;
            doVertices = true;
            fprintf(stderr,"Printing the maximum amount of vertices at min dist found so far\n");
            break;
        }
        if (strcmp(argv[i], "-e") == 0) {
            doApprox = true;
            doVertices = false;
            fprintf(stderr,"Printing the maximum amount of edges at min dist found so far\n");
            break;
        }
        if (strcmp(argv[i], "-g") == 0) {
            doGreedy = true;
            fprintf(stderr,"Greedily (not exact) searching for vertices and edges at min dist\n");
            break;
        }
        fprintf(stderr,"Only showing the exact results\n");
    }


    //  Start looping over lines of stdin.
    char *graphString = NULL;
    size_t size;
    int ctrAllGraphs = 0;
    while(getline(&graphString, &size, stdin) != -1) {

        // reset bestMaxVerticesSoFar and bestMaxEdgesSoFar
        bestMaxVerticesSoFar = 0;
        bestMaxEdgesSoFar = 0;

        ctrAllGraphs++;

        n = graphsize(graphString);
        maxm = SETWORDSNEEDED(n);

        // Allocations and graph information
        DYNALLSTAT(graph,g,g_sz);
        DYNALLOC2(graph,g,g_sz,maxm,n,"malloc");
        stringtograph(graphString,g,maxm);
        theGirth = girth(g, maxm, n);
        diamstats(g, maxm, n, &radius, &diameter);
        degstats(g, maxm, n, &es, &r, &mincount, &maxdeg, &maxcount, &eulerian);

        if (!allocateArrays(n,r))
            return -1;
        setEdges(g,n);

        // Compute distances
        for (int v = 0; v < n; v++) {
            find_dist(g, maxm, n , v, dist[v]);
        }
        fprintf(stderr,"n=%d (%d,%d) diam=%d\n", n, r,theGirth, diameter);

        int minDist = (theGirth + 1) / 2;

        if (!doApprox) {

            int maxNumVerticesDistantEnough = 0;

            int* numVerticesDistEnough = (int*)malloc(n * sizeof(int));
            if (numVerticesDistEnough == NULL) {
                fprintf(stderr,"Memory allocation numVerticesDistEnough failed\n");
                exit(-1);
            }

            # pragma omp parallel for
            for (int v = 0; v < n; v++) {
                currVerticesAll[v][0] = v;
                if (doGreedy)
                    numVerticesDistEnough[v] = checkDistancesGreedy(g, n, minDist, v);
                else
                    numVerticesDistEnough[v] = checkDistances(g, n, minDist, v, 1);
            }

            for (int v = 0; v < n; v++) {
                if (numVerticesDistEnough[v] > maxNumVerticesDistantEnough)
                    maxNumVerticesDistantEnough = numVerticesDistEnough[v];
            }

            fprintf(stderr,"num vertices at pairwise distance at least %d: %d\n",minDist,maxNumVerticesDistantEnough);


            int* numEdgesDistEnough = (int*)malloc(numEdges * sizeof(int));
            if (numEdgesDistEnough == NULL) {
                fprintf(stderr,"Memory allocation numEdgesDistEnough failed\n");
                exit(-1);
            }

            int maxNumEdgesDistantEnough = 0;

            # pragma omp parallel for
            for (int e = 0; e < numEdges; e++) {
                currEdgesAll[e][0] = e;
                if (doGreedy) {
                    numEdgesDistEnough[e] = checkDistancesEdgesGreedy(g, numEdges, minDist, e);
                } else
                    numEdgesDistEnough[e] = checkDistancesEdges(g, numEdges, minDist, e, 1);
            }

            for (int e = 0; e < numEdges; e++) {
                if (numEdgesDistEnough[e] > maxNumEdgesDistantEnough)
                    maxNumEdgesDistantEnough = numEdgesDistEnough[e];
            }

            fprintf(stderr,"num edges at pairwise distance at least %d: %d\n",minDist,maxNumEdgesDistantEnough);

            // Free allocations
            freeArrays(n);
            free(numVerticesDistEnough);
            free(numEdgesDistEnough);


        } else if (doVertices) {
            fprintf(stderr,"Best maxVerticesSoFar = %d\n", 1);
            #pragma omp parallel
            {
                # pragma omp for schedule(dynamic)
                for (int v = 0; v < n; v++) {
                    currVerticesAll[v][0] = v;
                    checkDistances(g, n, minDist, v, 1);
                }

                #pragma omp single
                {
                    freeArrays(n);
                }
            }
        } else {
            fprintf(stderr,"Best maxEdgesSoFar = %d\n", 1);
            # pragma omp parallel
            {
                # pragma omp for schedule(dynamic)
                for (int e = 0; e < numEdges; e++) {
                    currEdgesAll[e][0] = e;
                    checkDistancesEdges(g, numEdges, minDist, e, 1);
                }

                #pragma omp single
                {
                    freeArrays(n);
                }
            }
        }
    }

    fprintf(stderr,"Processed %d graphs\n", ctrAllGraphs);

    free(graphString);

    return 0;
}






;
