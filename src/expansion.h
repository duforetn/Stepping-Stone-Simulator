#include "random.h"
#include "pops.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void initExpansion1D(double *Pops, int nSNP, int npop, int popStart, int popEnd);

void initExpansion2D(double *Pops, int nSNP, int n1, int n2, int popStart, int popEnd);

void expand1D(double *Pops, double *OutsidePop, double *Propagule, double *AF, int npop, int nSNP, double m1, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture);

void expand2D(double *Pops, double *OutsidePop, double *AF, int n1, int n2, int nSNP, double m1, double m2, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture);

void expandNewpop1D(double *Pops, double *AF, double *Propagule, int pop, int nSNP, int Ne, int ke, int admixture);

void expandNewpop2D(double *Pops, double *AF, int pop, int nSNP, int n1, int n2, int Ne, int ke, int admixture);

