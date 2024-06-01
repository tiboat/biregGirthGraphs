#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <unistd.h> // For sleep function

#include "biregBounds.c"
#include "graphBiReg.c"
#include "../nauty/gutils.h"


/************************************************
 *  - Biregular graph with given girth generation algorithm -
 ***********************************************/

// Based on GenK2
#define USAGE \
"\nUsage:./biregGirthGen n r m g l0 l1 l2 phase res mod\n\n"

// A lot of code is taken from: https://github.com/JarneRenders/GenK2Hypohamiltonian/blob/main/genK2Hypohamiltonian.c
// Parts that are literally taken from there or modified are shown with GenK2

// Name of the audit files produced in the first phase, to read in the second phase
#define AUDIT_FILE_STR "audit/audit%d.txt"


/************************************************
 * Data structures
 ***********************************************/

//  Struct keeping counters
// based on GenK2
struct counters {
    long long unsigned int nOfGraphsVisited;
    long long unsigned int nOfGraphsVisitedPerLevel[2*MAXN];
    long long unsigned int nOfBiRegularGraphsFound;
    long long unsigned int nOfPrunedAddEdge;
    long long unsigned int nOfPrunedAddEdgePerLevel[2*MAXN];
    long long unsigned int nOfPrunedPruningRule;
    long long unsigned int nOfPrunedPruningRulePerLevel[2*MAXN];
    long long unsigned int nOfPrunedAlreadyVisited;
    long long unsigned int nOfPrunedPruningRuleMadeIllegal;
};


//  Struct for passing options.
// based on GenK2
struct options {
    // minDeg < maxDeg
    int minDeg, maxDeg;
    int girth;
    // levels for two-phase isomorphic pruning
    // level0 > level1 > level2
    int level0, level1, level2;
    // Minimal distance that two vertices of degree maxDeg should have
    int minDistMaxDeg;
    bool minDistCheckEnabled;
    bool l2PruningFirst;
    bool firstPhase;
    int res, mod;
};




/************************************************
 * Definition's for Nauty's splay tree (from GenK2)
 ***********************************************/

typedef struct SPLAYNODE {
    graph* canonForm;
    struct SPLAYNODE *left, *right, *parent;
} SPLAYNODE;

// Graphs in the search tree for which adding the edges
// to the current vertex all have been tried and
// the amount of added edges is not larger
// then level1
SPLAYNODE *completedGraphs[MAXEDGES] = {NULL};
int nrCompletedGraphs;

// Graphs called with the backtrack method
SPLAYNODE *visitedGraphs[MAXEDGES] = {NULL};
int nrVisitedGraphs;

#define SCAN_ARGS

#define ACTION(p)

#define INSERT_ARGS , graph gCan[], int numberOfVertices, bool *isPresent

int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices);

#define COMPARE(p) compareSplayNodeToGraph(p, gCan, numberOfVertices);

#define PRESENT(p) {(*isPresent) = true;}

#define NOT_PRESENT(p) {p->canonForm = gCan; (*isPresent) = false;}

#define LOOKUP_ARGS , graph gCan[], int numberOfVertices

#include "../nauty/splay.c"


int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices) {
    return memcmp(p->canonForm, gCan, MAXM * numberOfVertices * sizeof(graph));
}


bool addToCompletedGraphs(struct graph *g, int numEdges) {
    bool isPresent;
    graph* gCan = malloc(sizeof(graph)*g->numberOfVertices*MAXM);
    createCanonicalForm(g, gCan);
    splay_insert(&completedGraphs[numEdges], gCan, g->numberOfVertices, &isPresent);
    if(isPresent) {
        free(gCan);
    } else {
        nrCompletedGraphs++;
    }
    return isPresent;
}

bool addToVisitedGraphs(struct graph *g) {
    bool isPresent;
    graph* gCan = malloc(sizeof(graph)*g->numberOfVertices*MAXM);
    createCanonicalForm(g, gCan);
    splay_insert(&visitedGraphs[g->numEdgesAfterTree], gCan, g->numberOfVertices, &isPresent);
    if(isPresent) {
        free(gCan);
    } else {
        nrVisitedGraphs++;
    }
    return isPresent;
}


int freeSplayNode(SPLAYNODE *splayNode) {
    if(splayNode != NULL) {
        int l = freeSplayNode(splayNode->left);
        int r = freeSplayNode(splayNode->right);
        free(splayNode->canonForm);
        free(splayNode);
        return r+l+1;
    }
    return 0;
}

int freeSplayTrees() {
    int sum = 0;
    for(int i = 0; i < MAXEDGES; i++) {
        sum += freeSplayNode(completedGraphs[i]);
        completedGraphs[i] = NULL;
        freeSplayNode(visitedGraphs[i]);
        visitedGraphs[i] = NULL;
    }
    return sum;
}


/************************************************
 * Functions for writing auditing files
 ***********************************************/

#define DELIMITER "--"
#define LEN_DELIMITER 2  \
// Change the length of the DELIMITER accordingly if you would change the DELIMITER!

#define MAX_RETRIES 5
#define RETRY_DELAY_SECONDS 2

int nOfSubSearches = 0;

// Based on output of Chat-GPT
char* getIndexMaxDegString(bitset indexMaxDeg) {
    // Calculate the length of the string
    int strLength = 0;
    forEach(v, indexMaxDeg) {
        strLength += snprintf(NULL, 0, "%d ", v);
    }

    // Allocate memory for the string
    char *str = (char *)malloc((strLength + 2) * sizeof(char));
    if (str == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // Build the string
    int offset = 0;
    int ctr = 0;
    forEach(v, indexMaxDeg) {
        int written = snprintf(str + offset, strLength + 1 - offset, "%d", v);
        if (written < 0) {
            printf("Error occurred while formatting string.\n");
            free(str);
            exit(1);
        }
        offset += written;
        ctr++;
        if (ctr < size(indexMaxDeg)) {
            str[offset++] = ' '; // Add space delimiter
        }
    }
    str[offset++] = '\n';
    str[offset] = '\0'; // Null terminator

    return str;
}


// Based on output of Chat-GPT
void writeAdjacencyList(FILE *file, struct graph *g) {
    for (int i = 0; i < g->numberOfVertices; i++) {
        bool first = true;
        forEach(neigh, g->adjacencyList[i]) {
            if (first) {
                fprintf(file, "%d", neigh);
                first = false;
            } else {
                fprintf(file, " %d", neigh);
            }
        }
        fprintf(file, "\n");
    }
}


void writeAuditFile(struct graph *g, int currentVertex, int currentDeg, int searchDepth, int mod) {
    int res = nOfSubSearches % mod;

    char legalChoices[7*MAXEDGES];
    int currLen = 0;
    for (int i = 0; i < g->numberOfVertices; i++) {
        forEach(j, g->legalChoices[i]) {
            // Based on output of ChatGPT
            int sizeLine = snprintf(legalChoices + currLen,
                                    sizeof(legalChoices) - currLen, "%d %d  ", i, j);
            currLen += sizeLine;
        }
    }
    legalChoices[currLen] = '\0';

    char* indexMaxDegStr = getIndexMaxDegString(g->indexMaxDeg);

    // Calculate the length of the fileName
    int fileNameLength = snprintf(NULL, 0, AUDIT_FILE_STR, res) + 1;
    // Dynamically allocate memory for fileName
    char *fileName = (char *)malloc(fileNameLength * sizeof(char));
    if (fileName == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
    // Using sprintf to format the string
    sprintf(fileName, AUDIT_FILE_STR, res);

    // Retrying code based on Chat-GPT
    // Why retrying? When running the code for larger graphs, it might be that you write more than a million
    // of audit files. For this quantity of files, the chance of something going wrong here becomes not negligible
    // anymore, especially on complex hardware systems, like supercomputers.
    FILE *file;
    int attempt = 0;

    while (attempt < MAX_RETRIES) {
        // Write file if first time writing to it, otherwise append
        if(nOfSubSearches < mod)
            file = fopen(fileName, "w");
        else
            file = fopen(fileName, "a");

        if (file != NULL) {
            break;
        } else {
            fprintf(stderr, "Error opening audit%d.txt for writing. Attempt %d/%d\n",
                    res, attempt+1, MAX_RETRIES);
            attempt++;
            sleep(RETRY_DELAY_SECONDS); // Wait before retrying
        }
    }

    if (file == NULL) {
        // If the loop completes without opening the file, print error and return NULL
        fprintf(stderr, "Failed to open audit%d.txt file after %d attempts.\n",
                res, MAX_RETRIES);
        exit(-1);
    }

    writeAdjacencyList(file, g);
    fprintf(file, "%d %d %d %d\n", g->numEdgesAfterTree, currentVertex, currentDeg, searchDepth);
    fprintf(file, "%s", indexMaxDegStr);
    fprintf(file, "%s\n", legalChoices);
    fprintf(file, "%s\n", DELIMITER); // to separate different sub-branches
    fclose(file);

    // free the allocated memory
    free(fileName);
    free(indexMaxDegStr);

    nOfSubSearches++;
}



/************************************************
 * Main graph generator functions
 ***********************************************/


double timeSpentPrunable = 0;

// Returns true if graph g with next edge (v1,v2) is prunable or not according to the pruning rule
// of 'Computational determination of (3,11) and (4,7) cages'
// (https://arxiv.org/pdf/1009.3608.pdf)
bool prunable(struct graph *g, int v1, int v2) {

    // Only adds and deletes edges of nauty graph for efficiency
    ADDONEEDGE(g->nautyGraph, v1, v2, MAXM);

    // Look at all subsets of addedEdges of graph g, with an edge between v1 and v2 added
    // Finding all subsets based on ChatGPT
    int numSubsets = 1 << g->numEdgesAfterTree; // Number of subsets (2^n)
    int edgesToRemove[g->numEdgesAfterTree];
    int numEdgesToRemove;
    for (int subsetMask = 0; subsetMask < numSubsets; subsetMask++) {
        numEdgesToRemove = 0;
        for (int j = 0; j < g->numEdgesAfterTree; j++) {
            if ((subsetMask & (1 << j)) == 0) {
                // Remove edge from g
                REMOVEONEEDGE(g->nautyGraph, g->addedEdges[j].v1, g->addedEdges[j].v2, MAXM);
                // Keep track of which edge from g is removed and how many
                edgesToRemove[numEdgesToRemove] = j;
                numEdgesToRemove++;
            }
        }

        // Get canonical form
        createCanonicalForm(g, g->gCan);

        // Returns the splaynode of the graph that is isomorphic to g if it exists in te splay tree
        // otherwise returns NULL
        SPLAYNODE* isomorphicG = splay_lookup(&completedGraphs[g->numEdgesAfterTree + 1 - numEdgesToRemove],
                                              g->gCan, g->numberOfVertices);

        // Add edges back again
        for (int j = 0; j < numEdgesToRemove; j++) {
            ADDONEEDGE(g->nautyGraph, g->addedEdges[edgesToRemove[j]].v1, g->addedEdges[edgesToRemove[j]].v2, MAXM);
        }

        // If an isomorphic graph is found in the splay tree, then edge (v1,v2) is prunable
        if (isomorphicG != NULL) {
            REMOVEONEEDGE(g->nautyGraph, v1, v2, MAXM);
            return true;
        }
    }

    // No isomorphic graph found -> edge (v1,v2) is not prunable
    REMOVEONEEDGE(g->nautyGraph, v1, v2, MAXM);
    return false;
}



// Recursive backtracking function
void backtrack(struct graph *g, struct options *options, struct counters *counters, int currentVertex,
               int currentDeg, bool allHaveAtLeastMinDegPhase, int searchDepth, int comb[]) {
    counters->nOfGraphsVisited++;
    counters->nOfGraphsVisitedPerLevel[searchDepth]++;

    /* Might be useful to uncomment underneath lines when debugging or experimenting with the code :-) */
//    fprintf(stderr, "# visited graphs = %llu\n", counters->nOfGraphsVisited);
//    printAdjacencyList(g);
//    printLegalChoicesMaxDeg(g);

    // If we are in the allHaveAtLeastMinDegPhase, we can call the backtrack function again without
    // adding an edge. To prevent reprinting the graph and not checking the visited graphs
    // we perform a check.
    if(!allHaveAtLeastMinDegPhase || deg(g,currentVertex) != g->minDeg) {
        // Returns true if this graph has already been visited
        if(addToVisitedGraphs(g)) {
            counters->nOfPrunedAlreadyVisited++;
            return;
        }

        if (isBiRegular(g, options->minDeg, options->maxDeg)) {
            counters->nOfBiRegularGraphsFound++;
            createCanonicalForm(g, g->gCan);
            writeToG6(g->gCan, g->numberOfVertices);
        }
    }

    if(options->firstPhase && g->numEdgesAfterTree == options->level0) {
        // At this point, you are at the end of phase 1 => write the current
        // graph to a file.
        writeAuditFile(g, currentVertex, currentDeg, searchDepth, options->mod);
        return;
    }

    if(!allHaveAtLeastMinDegPhase) {
        // Go to allHaveAtLeastMinDegPhase if all vertices have at least deg minDeg
        if (allVerticesAtLeastMinDeg(g)) {
            allHaveAtLeastMinDegPhase = true;
        } else {
            if (deg(g,currentVertex) == g->minDeg) {
                // Only choose vertex with degree < minDeg
                currentVertex = getVertexMinLegalChoiceMinDeg(g);

                if (currentVertex == -1) {
                    // At this point there exists a vertex, not yet of degree minDeg,
                    // that has no legal choices left. We abort.
                    return;
                }
            }
        }
    }

    bitset legalChoicesChanged = EMPTY;

    if(allHaveAtLeastMinDegPhase) {
        if(deg(g, currentVertex) == currentDeg) {
            // Try to find a vertex that has degree > minDeg, but < maxDeg
            currentVertex = getVertexBetweenMinMaxDeg(g);

            if (currentVertex == -1) {
                // If there exists a vertex between degree minDeg and maxDeg, but returned -1,
                // we know that the vertex or vertices don't have enough choices left to complete
                for(int i = g->minDeg+1; i < g->maxDeg; i++) {
                    if (size(g->verticesOfDeg[i]) > 0)
                        return;
                }

                // At this point there is no vertex of degree between minDeg and maxDeg
                // We are going to finish the vertices of deg minDeg
                while(true) {
                    currentVertex = getVertexMinLegalChoiceMaxDeg(g);

                    if (currentVertex == -1) {
                        // It is not possible to finish any vertex of degree minDeg, we abort.
                        resetLegalChoices(g, searchDepth, &legalChoicesChanged);
                        return;
                    }

                    // Try to complete the graph where this vertex has degree either minDeg or maxDeg

                    // Only try to complete to maxDeg if there are at least maxDeg-minDeg legal choices
                    if (size(g->legalChoices[currentVertex]) >= g->maxDeg - g->minDeg)
                        backtrack(g, options, counters, currentVertex, g->maxDeg, allHaveAtLeastMinDegPhase,
                                  searchDepth+1, comb);

                    if (!makeVertexIllegal2(g, currentVertex, searchDepth, &legalChoicesChanged)) {
                        resetLegalChoices(g, searchDepth, &legalChoicesChanged);
                        return;
                    }
                }
            }
            // Set the degree for the currentVertex to maxDeg, because we necessarily need to finish it to that
            // degree
            currentDeg = g->maxDeg;
        }
    }

    // Extra sanity check
    if (currentVertex == -1) {
        fprintf(stderr, "currentVertex is -1");
        exit(-1);
    }

    // Extra check to see if tree is not empty, otherwise looping is useless
    if(options->l2PruningFirst && g->numEdgesAfterTree <= options->level2 && nrCompletedGraphs != 0) {
        forEach(v1, g->indexMaxDeg) {
            if(v1 != currentVertex) {
                forEachAfterIndex(v2, g->legalChoices[v1], v1) {
                    clock_t start = clock();
                    bool isPrunable = prunable(g, v1, v2);
                    clock_t end = clock();
                    timeSpentPrunable += (double)(end - start) / CLOCKS_PER_SEC;
                    if (isPrunable) {
                        counters->nOfPrunedPruningRule++;
                        ADDONEEDGE(g, v1, v2, MAXM);
                        addToCompletedGraphs(g, g->numEdgesAfterTree + 1);
                        REMOVEONEEDGE(g, v1, v2, MAXM);
                        if (!makeEdgeIllegal(g, v1, v2, searchDepth, &legalChoicesChanged)) {
                            counters->nOfPrunedPruningRuleMadeIllegal++;
                            resetLegalChoices(g, searchDepth, &legalChoicesChanged);
                            return;
                        }
                    }
                }
            }
        }
    }

    bitset potentialNeighbrs = g->legalChoices[currentVertex];

    forEach(otherV, potentialNeighbrs) {
        // otherV might have been made illegal in the meanwhile
        if(!contains(g->legalChoices[currentVertex], otherV)) {
            continue;
        }

        if(isolatedLexMin(g, currentVertex, otherV)) {
            bitset legalChoicesChangedAddEdge = EMPTY;
            int numDistChanged = -1;

            // Add edge
            bool stillLegal = addEdgeWithDistMatrixLegal2(g,  options->girth, currentVertex, otherV,
                                                          searchDepth, &legalChoicesChangedAddEdge,
                                                          &numDistChanged, comb);

            if(stillLegal && options->minDistCheckEnabled &&
               !minDistMaxDegSatisfied(g, options->minDistMaxDeg, currentVertex, otherV,
                                       searchDepth, &legalChoicesChangedAddEdge)) {
                stillLegal = false;
            }

            struct edge newEdge = { .v1 = currentVertex, .v2 = otherV };
            g->addedEdges[g->numEdgesAfterTree] = newEdge;
            g->numEdgesAfterTree++;

            // 4. Early abortion -> don't go further if not legal anymore
            if(stillLegal) {

                if (!allHaveAtLeastMinDegPhase && g->numEdgesAfterTree <= options->level1 && nrCompletedGraphs != 0) {
                    REMOVEONEEDGE(g, currentVertex, otherV, MAXM);
                    g->numEdgesAfterTree--;
                    clock_t start = clock();
                    bool isPrunable = prunable(g, currentVertex, otherV);
                    clock_t end = clock();
                    timeSpentPrunable += (double)(end - start) / CLOCKS_PER_SEC;
                    ADDONEEDGE(g, currentVertex, otherV, MAXM);
                    g->numEdgesAfterTree++;
                    if (isPrunable) {
                        stillLegal = false;
                        counters->nOfPrunedPruningRule++;
                        counters->nOfPrunedPruningRulePerLevel[searchDepth]++;
                    }
                }

                if(stillLegal) {
                    backtrack(g, options, counters, currentVertex, currentDeg, allHaveAtLeastMinDegPhase,
                              searchDepth+1, comb);
                }

            } else {
                counters->nOfPrunedAddEdge++;
                counters->nOfPrunedAddEdgePerLevel[searchDepth]++;
            }

            if (g->numEdgesAfterTree <= options->level1) {
                // Subsearch is completed -> add graph to splay tree
                addToCompletedGraphs(g, g->numEdgesAfterTree);
            }

            // Remove edge
            removeEdgeWithDistMatrixAndReset(g, currentVertex, otherV, searchDepth,
                                             &legalChoicesChangedAddEdge, &numDistChanged);
            g->numEdgesAfterTree--;

            // Still leave added edge to be illegal
            if (!makeEdgeIllegal(g, currentVertex, otherV, searchDepth, &legalChoicesChanged)) {
                resetLegalChoices(g, searchDepth, &legalChoicesChanged);
                return;
            }
         }
    }

    // Extra check to see if tree is not empty, otherwise looping is useless
    if(!options->l2PruningFirst && g->numEdgesAfterTree <= options->level2 && nrCompletedGraphs != 0) {
        forEach(v1, g->indexMaxDeg) {
            if(v1 != currentVertex) {
                forEachAfterIndex(v2, g->legalChoices[v1], v1) {
                    clock_t start = clock();
                    bool isPrunable = prunable(g, v1, v2);
                    clock_t end = clock();
                    timeSpentPrunable += (double)(end - start) / CLOCKS_PER_SEC;
                    if (isPrunable) {
                        counters->nOfPrunedPruningRule++;
                        ADDONEEDGE(g, v1, v2, MAXM);
                        addToCompletedGraphs(g, g->numEdgesAfterTree + 1);
                        REMOVEONEEDGE(g, v1, v2, MAXM);
                        if (!makeEdgeIllegal(g, v1, v2, searchDepth, &legalChoicesChanged)) {
                            counters->nOfPrunedPruningRuleMadeIllegal++;
                            resetLegalChoices(g, searchDepth, &legalChoicesChanged);
                            return;
                        }
                    }
                }
            }
        }
    }

    // Reset legalChoice changes from pruning rule
    resetLegalChoices(g, searchDepth, &legalChoicesChanged);
}





// Helper function for createMooreTreeGivenConfig
// Special case: returns the highest level if only the top (root) level
// has non-zero in config
int getHighestNonZeroLevel(int* config, int girth) {
    int highestLevel = EVEN(girth) ? girth/2-1 : girth/2;
    int highestNonZeroLevel = highestLevel;
    for (int i = highestLevel-1; i > 0; i--) {
        if (config[i] != 0) {
            highestNonZeroLevel = i;
        }
    }
    return highestNonZeroLevel;
}

// Create the Moore tree for the given Moore tree config
void createMooreTreeGivenConfig(struct graph *g, struct options *options, int* config) {
    int minDeg = options->minDeg;
    int maxDeg = options->maxDeg;
    int girth = options->girth;
    int minDistMaxDeg = options->minDistMaxDeg;

    int treeHeight = floor((girth - 1) >> 1);
    int startVertexCurrentLevel, startVertexNextLevel;

    // Create root and children of root in the tree

    add(g->verticesAtLevel[0],0);
    g->levelOfVertex[0] = 0;

    if(ODD(girth)) {
        for (int i = 1; i <= maxDeg; i++) {
            addEdgeWithDistMatrixLegal(g, minDistMaxDeg, girth, 0, i);
        }
        startVertexCurrentLevel = 1;
        startVertexNextLevel = maxDeg + startVertexCurrentLevel;

    } else {
        addEdgeWithDistMatrixLegal(g, minDistMaxDeg, girth, 0, 1);
        for (int i = 2; i <= maxDeg; i++) {
            addEdgeWithDistMatrixLegal(g, minDistMaxDeg, girth, 0, i);
        }
        for (int i = maxDeg + 1; i <= minDeg + maxDeg - 1; i++) {
            addEdgeWithDistMatrixLegal(g, minDistMaxDeg, girth, 1, i);
        }

        if (config[0] == 1) {
            removeElement(g->indexMaxDeg, 1);
            makeVertexIllegal(g, 1);
        }

        add(g->verticesAtLevel[0],1);
        g->levelOfVertex[1] = 0;

        startVertexCurrentLevel = 2;
        startVertexNextLevel = maxDeg + minDeg;
    }

    int beginVertexNextLevel, numVerticesDegMLevel, numEdgesUnderneath;
    int highestLevelDeg2 = getHighestNonZeroLevel(config, girth);

    // create rest of tree
    for (int h = 1; h < treeHeight; h++) {
        beginVertexNextLevel = startVertexNextLevel;
        for (int i = startVertexCurrentLevel; i < beginVertexNextLevel; i++) {
            add(g->verticesAtLevel[h],i);
            g->levelOfVertex[i] = h;
            numEdgesUnderneath = minDeg - 1;
            for (int j = startVertexNextLevel; j < startVertexNextLevel + numEdgesUnderneath; j++) {
                addEdgeWithDistMatrixLegal(g, minDistMaxDeg, girth, i, j);
            }
            if (h < highestLevelDeg2) {
                removeElement(g->indexMaxDeg, i);
                makeVertexIllegal(g, i);
            }
            startVertexNextLevel += numEdgesUnderneath;
        }
        startVertexCurrentLevel = beginVertexNextLevel;
    }

    int startVertexLastToCheckLevel, lastVertexLastToCheckLevel;

    for (int h = highestLevelDeg2; h < treeHeight; h++) {
        if (ODD(girth)) {
            startVertexLastToCheckLevel = 1 + maxDeg * powerSum(minDeg - 1, h - 2);
            lastVertexLastToCheckLevel = maxDeg * powerSum(minDeg - 1, h - 1);
        } else {
            startVertexLastToCheckLevel = 2 + (maxDeg + minDeg - 2) * powerSum(minDeg - 1, h - 2);
            lastVertexLastToCheckLevel = 1 + (maxDeg + minDeg - 2) * powerSum(minDeg - 1, h - 1);
        }
        int numDeg2 = config[h];

        if (numDeg2 == 0) {
            for (int j = startVertexLastToCheckLevel; j <= lastVertexLastToCheckLevel; j++) {
                // Make impossible to be of maxDeg
                removeElement(g->indexMaxDeg, j);
                makeVertexIllegal(g, j);
            }
        }
    }

    /* Uncomment underneath prints to obtain more information */

//    fprintf(stderr,"Adjacency list Moore tree: \n");
//    printVerticesAtLevel(g);
//    printAdjacencyList(g);
//    printIndexMinDeg(g);
//    printIndexMaxDeg(g);
}


// Function to call to perform first phase
void generateGraphsFromMooreTree(int numVertices, struct options *options, struct counters *counters, int* config, int t) {
    struct graph g = makeEmptyGraph(numVertices, options->minDeg, options->maxDeg, t);
    createMooreTreeGivenConfig(&g, options, config);
    int startVertex = getVertexMinLegalChoiceMinDeg(&g);

    backtrack(&g, options, counters, startVertex, options->minDeg,
              false, 0, config);
    deleteGraph(&g);
    freeSplayTrees();
}



/************************************************
 * Second phase functions
 ***********************************************/

void doSecondPhase(struct graph *g, int currentVertex, int currentDeg, int searchDepth,
                   struct options *options, struct counters *counters, int comb[]) {

    backtrack(g, options, counters, currentVertex, currentDeg, allVerticesAtLeastMinDeg(g),
              searchDepth, comb);

    deleteGraph(g);
    freeSplayTrees();
}


// Read audit file and sets forbidden edges
void readAuditFileAndDoSecondPhase(int numVertices, int res, struct options *options, struct counters *counters,
                                   int comb[]) {
    // Calculate the length of the fileName
    int fileNameLength = snprintf(NULL, 0, AUDIT_FILE_STR, res) + 1;
    // Dynamically allocate memory for fileName
    char *fileName = (char *)malloc(fileNameLength * sizeof(char));
    if (fileName == NULL) {
        printf("Memory allocation failed.\n");
        exit(-1);
    }
    // Using sprintf to format the string
    sprintf(fileName, AUDIT_FILE_STR, res);

    // Retrying code based on Chat-GPT
    // Why retrying? When running the code for larger graphs, it might be that you write more than a million
    // of audit files. For this quantity of files, the chance of something going wrong here becomes not negligible
    // anymore, especially on complex hardware systems, like supercomputers.
    FILE *file;
    int attempt = 0;

    while (attempt < MAX_RETRIES) {
        file = fopen(fileName, "r");
        if (file != NULL) {
            break;
        } else {
            fprintf(stderr, "Error opening audit%d.txt for reading. Attempt %d/%d\n",
                    res, attempt+1, MAX_RETRIES);
            attempt++;
            sleep(RETRY_DELAY_SECONDS); // Wait before retrying
        }
    }

    if (file == NULL) {
        // If the loop completes without opening the file, print error and return NULL
        fprintf(stderr, "Failed to open audit%d.txt file after %d attempts.\n",
                res, MAX_RETRIES);
        exit(-1);
    }


    // Keep on reading audit of sub-searches and computing them until file is 'finished'
    while(true) {
        if (feof(file)) {
            fprintf(stderr, "End of file reached\n%d searches done\n", nOfSubSearches);
            break;
        }

        int t = (options->girth % 2) ? options->girth/2: options->girth/2 - 1;
        struct graph g = makeEmptyGraph(numVertices, options->minDeg, options->maxDeg, t);
        int currentVertex, currentDeg, searchDepth;

        // Read the adjacency list from the file
        char adjListLine[MAXN]; // MAXN is a rough estimate, but should usually be enough
        for (int i = 0; i < g.numberOfVertices; i++) {
            if (fgets(adjListLine, sizeof(adjListLine), file) != NULL) {
                // Count the number of integers in the line
                int count = 0;
                for (int j = 0; adjListLine[j] != '\n'; j++) {
                    if (adjListLine[j] == ' ') {
                        count++;
                    }
                }
                if (count > 0) {
                    // Parse integers from the line
                    char *token = strtok(adjListLine, " ");
                    while (token != NULL) {
                        addEdge((&g),i,atoi(token));
                        token = strtok(NULL, " ");
                    }
                }
            } else {
                fprintf(stderr, "Error reading file.\n");
                fclose(file);
                exit(-1);
            }
        }

        fscanf(file, "%d %d %d %d\n", &g.numEdgesAfterTree, &currentVertex, &currentDeg, &searchDepth);


        // read indexMaxDeg
        g.indexMaxDeg = EMPTY;
        // Based on output of ChatGPT
        char buffer[7*MAXEDGES];
        if (fgets(buffer, sizeof(buffer), file) != NULL) {
            // Parse integers from the line
            char *token = strtok(buffer, " ");
            while (token != NULL) {
                add(g.indexMaxDeg, atoi(token));
                token = strtok(NULL, " ");
            }
        } else {
            fprintf(stderr, "Error reading file.\n");
            fclose(file);
            exit(-1);
        }

        // Based on output of ChatGPT
        // read legal choices
        for (int i = 0; i < g.numberOfVertices; i++) {
            g.legalChoices[i] = EMPTY;
        }
        int v1, v2;
        while (fscanf(file, "%d %d", &v1, &v2) == 2) {
            add(g.legalChoices[v1], v2);
        }

        // Set distMatrix
        for (int i = 0; i < g.numberOfVertices; i++) {
            find_dist(g.nautyGraph, MAXM, g.numberOfVertices, i, g.distMatrix[i]);
        }


        // Read delimiter
        char delim[LEN_DELIMITER];
        // This does not check that the delimiter is correct!
        fscanf(file, "%s\n", delim);

        doSecondPhase(&g, currentVertex, currentDeg, searchDepth, options, counters, comb);

        if (nOfSubSearches % 500 == 0 && nOfSubSearches != 0)
            fprintf(stderr, "Sub-search %d done\n", nOfSubSearches);

        nOfSubSearches++;
    }

    // Close the file
    fclose(file);
    free(fileName);
}


/************************************************
 * Execution information printing
 ***********************************************/


void printExecutionInformation(int numVertices, double timeSpent, struct options *options, struct counters *counters) {
    if(options->firstPhase && nOfSubSearches != 0)
        fprintf(stderr, "Wrote %d sub-searches to %d files\n",
                nOfSubSearches, options->mod);

    /* Uncomment these prints to obtain more information on program execution */

    fprintf(stderr,"Number of graphs visited with %d vertices: %llu\n",
            numVertices, counters->nOfGraphsVisited);
//    fprintf(stderr,"Number of unique graphs visited with %d vertices: %llu\n",
//            numVertices, counters->nOfGraphsVisited-counters->nOfPrunedAlreadyVisited);
//    fprintf(stderr,"Number of graphs visited per level:\n");
//    for (int i = 0; i < 2*MAXN; i++) {
//        fprintf(stderr,"%llu ",
//                counters->nOfGraphsVisitedPerLevel[i]);
//    }
//    fprintf(stderr,"\n");
    fprintf(stderr, "Number of ({%d,%d};%d)-graphs found with %d vertices: %llu\n", options->minDeg,
            options->maxDeg, options->girth, numVertices, counters->nOfBiRegularGraphsFound);
//    fprintf(stderr,"Number of graphs pruned because already visited: %llu\n",
//            counters->nOfPrunedAlreadyVisited);
//    fprintf(stderr,"Number of graphs pruned when adding edge: %llu\n",
//            counters->nOfPrunedAddEdge);
//    fprintf(stderr,"Number of graphs pruned when adding edge per level:\n");
//
//    for (int i = 0; i < 2*MAXN; i++) {
//        fprintf(stderr,"%llu ",
//                counters->nOfPrunedAddEdgePerLevel[i]);
//    }
//    fprintf(stderr,"\n");
//
//    if (options->firstPhase) {
//        fprintf(stderr,"Number of graphs pruned with pruning rule: %llu\n",
//                counters->nOfPrunedPruningRule);
//        fprintf(stderr,"Number of graphs pruned with pruning rule per level:\n");
//        for (int i = 0; i < 2*MAXN; i++) {
//            fprintf(stderr,"%llu ",
//                    counters->nOfPrunedPruningRulePerLevel[i]);
//        }
//        fprintf(stderr,"\n");
//        fprintf(stderr,"Number of times pruning rule made graph illegal: %llu\n",
//                counters->nOfPrunedPruningRuleMadeIllegal);
//    }

    fprintf(stderr,"Time spent: %f s\n",
            timeSpent);

//    fprintf(stderr,"Time spent in createCanonicalForm: %f s\n",
//            timeCreateCanonical);
//    if (options->firstPhase)
//        fprintf(stderr,"Time spent in prunable: %f s\n",
//                timeSpentPrunable);

    fprintf(stderr, "Search ended successfully\n");
}



/************************************************
 * Main function
 ***********************************************/
 
 
// Helper function
bool getBooleanOutOfPhase(int phase) {
    if (phase == 1)
    	return true;
    if (phase == 2)
    	return false;
    fprintf(stderr, "Incorrect phase argument: %d\nShould be either 1 or 2", phase);
    exit(-1);
    return true;
}

int main(int argc, char ** argv) {
    if (argc != 11) {
        // Print the usage information and exit with an error code
        printf(USAGE);
        return 1;
    }

    int numVertices = atoi(argv[1]);

    struct options options = {
            .minDeg = atoi(argv[2]),
            .maxDeg = atoi(argv[3]),
            .girth = atoi(argv[4]),
            .level0 = atoi(argv[5]),
            .level1 = atoi(argv[6]),
            .level2 = atoi(argv[7]),
            .minDistCheckEnabled = true,
            .l2PruningFirst = true
    };

    int phase = atoi(argv[8]);
    options.firstPhase = getBooleanOutOfPhase(phase);
    
    options.res = atoi(argv[9]);
    options.mod = atoi(argv[10]);


    options.minDistMaxDeg =  minDistDegM(numVertices, options.minDeg, options.maxDeg, options.girth / 2,
                                         options.girth);

    // TODO
    fprintf(stderr, "\nOptions:\nn r m g l0 l1 l2: %d %d %d %d %d %d %d\n",
            numVertices, options.minDeg, options.maxDeg, options.girth, 
            options.level0, options.level1, options.level2);

    if(options.firstPhase)
        fprintf(stderr, "First Phase\n");
    else
        fprintf(stderr, "Second Phase: subbranches res %d\n", options.res);

//    fprintf(stderr,"min dist m = %d\n", options.minDistMaxDeg);

    // Initialize all counters to 0.
    struct counters counters = {};

    // From GenK2
    clock_t start = clock();

    int t = (options.girth % 2) ? options.girth/2: options.girth/2 - 1;
    if(options.firstPhase) {
        // For odd girth: 2t+1=g, for even girth: 2(t+1) = g
        configArray mooreTrees = getTreeCombs(numVertices, options.minDeg, options.maxDeg, options.girth, t);

        /* Uncomment two underneath lines to obtain more information on the possible Moore trees */
//        fprintf(stderr,"Moore trees\n");
//        printConfigArray(mooreTrees, t);

        int* combWithMaxes = getCombMaxEachLevel(mooreTrees, t);
        generateGraphsFromMooreTree(numVertices, &options, &counters, combWithMaxes, t);
        freeConfigArray(mooreTrees);
        free(combWithMaxes);
    } else {
        configArray mooreTrees = getTreeCombs(numVertices, options.minDeg, options.maxDeg, options.girth, t);
        int* combWithMaxes = getCombMaxEachLevel(mooreTrees, t);
        readAuditFileAndDoSecondPhase(numVertices, options.res, &options, &counters, combWithMaxes);
        freeConfigArray(mooreTrees);
        free(combWithMaxes);
    }


    // From GenK2
    clock_t end = clock();
    double timeSpent = (double)(end - start) / CLOCKS_PER_SEC;


    // Print information on program execution
    printExecutionInformation(numVertices, timeSpent, &options, &counters);

    return 0;
}
