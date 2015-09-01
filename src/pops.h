#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>


void InitPops(double **Pops, double **OutsidePop, double **AF, int npop, int nSNP);

void evolve1D(double *Pops, double *OutsidePop, double *AF, int npop, int nSNP, double m1, double minf, int Ne, int Loss);

void evolve2D(double *Pops, double *OutsidePop, double *AF, int npop, int n1, int n2, int nSNP, double m1, double m2, double minf, int Ne, int Loss);

int burninPop(double *Pops, double *OutsidePop, double *AF, int burnin, int npop, int n1, int n2, int m1, int m2, int minf, int nSNP, int Ne);

