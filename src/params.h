#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include "matrix.h"

void Welcome(int help);

int handleParams(int argc, char* argv[], int *npop, int *n1, int *n2, int *discardEdge, int *burnin, int *nSNP, double *m1, double *m2, double *minf, int *Ne, int *nSamples, int *ngen, int *downSample, int *calculatePsi, int *Loss, int *nextSample, int **Times, int *rangeExpansion, int *timeStart, int *timeEnd, int *popStart, int *popEnd, int *nSteps, int *lengthSteps, int *admixture, int *ke, char **outputName);

void getCoordinates(int **coords, int *Times, int npop, int n1, int n2, int discardEdge, int nSamples);
