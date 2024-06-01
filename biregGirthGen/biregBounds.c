#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * This file contains code related to computing
 *  - the minimum distande of degree m
 *  - the possible Moore trees
 *  of an ({r,m}:g)-graph.
 */


#define EVEN(n) (n % 2 == 0)
#define ODD(n) (n % 2 == 1)


/************************************************
 * Minimum distance vertices degree m
 ***********************************************/

//  Helper functions

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


// Returns the minimal distance vertices of degree m need to be
// according to 'Small bi-regular graphs of even girth'
// (https://doi.org/10.1016/j.disc.2015.10.009)
// For odd girth: 2t+1=g, for even girth: 2t = g
int minDistDegM(int numVertices, int r, int m, int t, int girth) {
    if (girth == 3 || girth == 4)
        return 1;

    int currNumVertices;
    int minDist = -1;
    for (int d = t; d > 0; d--) {
        if (d == t && girth % 2 == 0 && t % 2 == 1) {
            int s = t / 2;
            double temp1 = m * powerSum(r - 1, t - 2);
            double temp2 = (double) (m * ipow(r - 1, t - 1) + m) / (double) r;
            double temp3 =  ((double) m / (double) r - 1.0);
            double temp4 = 1 + m * powerSum(r - 1, s - 1);
            double temp = (temp1 + temp2) / (1.0 + temp3 / temp4);
            currNumVertices = ceil(temp);
        } else if (d < t && d != 1 && (d != 3 || girth % 2 == 1)) {
            if (girth % 2 == 0)
                currNumVertices = 1 + m * powerSum(r - 1, t - 2) +
                                (m - r) * powerSum(r - 1, t - 1 - d) + ipow(r - 1, t - 1);
            else
                currNumVertices = 1 + m * powerSum(r - 1, t - 1) +
                                  (m - r) * powerSum(r - 1, t - 1 - d);
        } else if (d == 1) {
            if (girth % 2 == 0)
                currNumVertices = ceil(2 + 2 * (m - 1) * (ipow(r - 1, t - 1) - 1) / (r - 2));
            else
                currNumVertices = ceil(2 + 2 * (m - 1) * (ipow(r - 1, t - 1) - 1) / (r - 2)
                                        + (m-1) * ipow(r-1,t-1));
        } else if (d == 3 && girth % 2 == 0) {
            double temp1 = (double)((ipow(r - 1, t - 1) - 1)) / (double)((r - 2));
            double temp2 = (double)(ipow(r - 1, t - 1));
            double temp3 = (double)(((m - r) * (r - 2) * ipow(r - 1, t - 2))) / (double)(r);
            double temp = 1.0 + (double)(m) * temp1 + temp2 + temp3;
            temp = ceil(temp);
            currNumVertices = (int) temp;
        } else {
            continue;
        }

        if (currNumVertices <= numVertices) {
            minDist = d;
        } else {
            if (minDist == -1)
                minDist = d+1;
        }

    }
    return minDist;
}


/************************************************
 * Computing possible configurations (=placements,
 * combinations) of vertices of degree m in the
 * Moore tree
 ***********************************************/


typedef struct {
    int** configs;
    size_t size;
} configArray;

void printConfig(int* config, int size) {
    for (int i = 0; i < size; i++) {
        fprintf(stderr,"%d\t", config[i]);
    }
    fprintf(stderr,"\n");
}

void printConfigArray(configArray configArray, int t) {
    for (int i = 0; i < configArray.size; i++) {
        printConfig(configArray.configs[i], t);
    }
}


int getNumVerticesLevel0(int g) {
    return (ODD(g) ? 1 : 2);
}

int getNumVerticesLevel1(int r, int m, int g, int config[]) {
    if(ODD(g))
        return m;
    return m + (config[0] == 1 ? r : m) - 2;
}

int getNumVerticesAtLevel(int r, int m, int g, int level, int* config) {
    if(level == 0)
        return getNumVerticesLevel0(g);
    int numVLevel1 = getNumVerticesLevel1(r, m, g, config);
    if(level == 1)
        return numVLevel1;

    int sum = numVLevel1 * ipow(r-1, level-1);
    for(int i = 1; i <= level-1; i++) {
        sum += (m-r) * config[i] * ipow(r-1, level-1-i);
    }
    return sum;
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


int* getCombMaxEachLevel(configArray configArray, int configSize) {
    int *combWithMaxes = (int *) malloc(configSize*sizeof(int));
    for (int j = 0; j < configSize; j++)
        combWithMaxes[j] = 0;

    for (int i = 0; i < configArray.size; i++)
        for (int j = 0; j < configSize; j++)
            if(configArray.configs[i][j] > combWithMaxes[j])
                combWithMaxes[j] = configArray.configs[i][j];

    return combWithMaxes;
}


// For odd girth: 2t+1=g, for even girth: 2(t+1) = g
configArray getTreeCombs(int n, int r, int m, int g, int t) {
    configArray configArray = {};
    int config[t];
    config[0] = 1;
    for (int i = 1; i < t; i++) {
        config[i] = 0;
    }
    int numVerticesComb, numVerticesAtLevel;
    int stop = 0;
    int currLevel = t-1;

    if(g == 3 || g == 4) {
        configArray.size++;
        configArray.configs = (int **)realloc(configArray.configs, configArray.size * sizeof(int *));
        int *newConfig = (int *)malloc(t * sizeof(int));
        for (int i = 0; i <= t-1; i++) {
            newConfig[i] = config[i];
        }
        configArray.configs[configArray.size-1] = newConfig;

    } else {
        while(!stop) {
            numVerticesComb = getNumVerticesComb(r, m, g, t, config);
            numVerticesAtLevel = getNumVerticesAtLevel(r, m, g, currLevel, config);

            if(config[currLevel] > numVerticesAtLevel || numVerticesComb > n) {
                do {
                    config[currLevel] = 0;
                    currLevel--;
                    if(currLevel == -1) {
                        stop = 1;
                        break;
                    }
                    config[currLevel]++;
                    numVerticesAtLevel = getNumVerticesAtLevel(r, m, g, currLevel, config);
                } while(config[currLevel] > numVerticesAtLevel);
                currLevel = t-1;

            } else {
                // Copy config to valid configs
                configArray.size++;
                configArray.configs = (int **)realloc(configArray.configs, configArray.size * sizeof(int *));
                int *newConfig = (int *)malloc(t * sizeof(int));
                for (int i = 0; i <= t-1; i++) {
                    newConfig[i] = config[i];
                }
                configArray.configs[configArray.size-1] = newConfig;
                config[currLevel]++;
            }
        }
    }

    return configArray;
}


void freeConfigArray(configArray configArray) {
    for (int i = 0; i < configArray.size; i++) {
        free(configArray.configs[i]);
    }
    free(configArray.configs);
}


