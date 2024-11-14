#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "../nauty/nauty.h"
#include "../nauty/gtools.h"
#include "../nauty/gutils.h"


/************************************************
 *  - Regular to biregular graph construction algorithm -
 ***********************************************/


// Based on GenK2
#define USAGE \
"\nUsage:./regToBiregAddEdgeExec r girth minN maxN\n\n"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define EVEN(n) (n % 2 == 0)
#define ODD(n) (n % 2 == 1)

struct edge { int v1,v2; };


/************************************************
 *  Global variables
 ***********************************************/

int leastGirthWanted;
bool foundBiregGraph = false;
int maxm;


//  Macro's for nauty representation
#define FOREACH(element,nautySet)\
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



/************************************************
*  Construction functions
***********************************************/


bool addEdgeCheckGirth(graph *g, int v1, int v2, int n) {
    bool success = false;
    ADDONEEDGE(g, v1, v2, maxm);
    int newGirth = girth(g, maxm, n);
    if(newGirth >= leastGirthWanted) {
        success = true;
        printf("n, new girth: %d, %d\n", n, newGirth);
        foundBiregGraph = true;
        char *g6String = ntos6(g,maxm,n);
        printf("%s", g6String);
    } else {
        REMOVEONEEDGE(g, v1, v2, maxm);
    }
    return success;
}


// Try the construction
void tryAddEdge(graph *g, int n, int r, int m, int **dist) {
    int minDist = leastGirthWanted-1;

    for (int v = 0; v < n; v++) {
        find_dist(g, maxm, n , v, dist[v]);
    }

    for (int v1 = 0; v1 < n && !foundBiregGraph; v1++) {
        for (int v2 = 0; v2 < n && !foundBiregGraph; v2++) {
            if (dist[v1][v2] >= minDist) {
                addEdgeCheckGirth(g, v1, v2, n);
            }
        }
    }
}




int main(int argc, char ** argv) {
    if (argc != 5) {
        // Print the usage information and exit with an error code
        printf(USAGE);
        return 1;
    }

    // Ask for r because of performance improvement
    int r = atoi(argv[1]);
    int m = r+1;

    leastGirthWanted = atoi(argv[2]);
    int minVertices = atoi(argv[3]);
    int maxVertices = atoi(argv[4]);

    int n,initialGirth;
    int radius, diameter;

    int lowerBound = minVertices;
    int upperBound = maxVertices;

    //  Start looping over lines of stdin.
    char *graphString = NULL;
    size_t size;
    int ctrUnfiltered = 0;
    int ctrAllGraphs = 0;
    while(getline(&graphString, &size, stdin) != -1) {
        ctrAllGraphs++;

        n = graphsize(graphString);

        if (n > upperBound)
            break;

        int t, potentialNewLowerBound;
        if (EVEN(leastGirthWanted)) {
            t = leastGirthWanted/2-1;
            int *config = (int *)malloc(t * sizeof(int));
            config[0] = 2 ;
            for (int i = 1; i < t; i++) {
                config[i] = 0;
            }
            potentialNewLowerBound = getNumVerticesComb(r, m, leastGirthWanted, t, config);
            free(config);
        } else {
            t = leastGirthWanted/2;
            int *config = (int *)malloc(t * sizeof(int));
            config[0] = 1;
            config[1] = 1;
            for (int i = 2; i < t; i++) {
                config[i] = 0;
            }
            potentialNewLowerBound = getNumVerticesComb(r, m, leastGirthWanted, t, config);
            free(config);
        }
        if (potentialNewLowerBound > lowerBound)
            lowerBound = potentialNewLowerBound;


        if (n >= lowerBound && n <= upperBound) {
            maxm = SETWORDSNEEDED(n);
            DYNALLSTAT(graph,g,g_sz);
            DYNALLOC2(graph,g,g_sz,maxm,n,"malloc");

            stringtograph(graphString,g,maxm);

            initialGirth = girth(g, maxm, n);
            if (initialGirth >= leastGirthWanted) {
                // Compute diameter
                diamstats(g, maxm, n, &radius, &diameter);
                int diameterBound = leastGirthWanted-1;
                if (diameter >= diameterBound) {
                    ctrUnfiltered++;

                    int **dist = (int **)malloc(n * sizeof(int *));
                    if (dist == NULL) {
                        printf("Memory allocation dist failed\n");
                        return -1;
                    }

                    for (int i = 0; i < n; i++) {
                        dist[i] = (int *)malloc(n * sizeof(int));
                        if (dist[i] == NULL) {
                            printf("Memory allocation failed for row %d!\n", i);
                            return -1;
                        }
                    }

                    printf("diameter=%d\n", diameter);
                    printf("initial girth=%d\n", initialGirth);
                    printf("initial n=%d\n", n);
                    tryAddEdge(g, n, r, m, dist);

                    for (int i = 0; i < n; i++) {
                        free(dist[i]);
                    }

                    free(dist);

                    if (foundBiregGraph)
                        break;
                }
            }
        }
    }

    printf("filtered %d graphs out of %d\n", ctrAllGraphs-ctrUnfiltered, ctrAllGraphs);

    free(graphString);

    return 0;
}






