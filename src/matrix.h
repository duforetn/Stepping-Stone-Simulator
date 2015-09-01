#include "../src_Lapack/lpk.h"
#include <string.h>
#include "random.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int findMax(int *matrix, int nrow, int ncol);

int findNext(int *matrix, int nrow, int ncol, int current);

int sampleData(double *Data, int *sampledData, double *Pops, int npop, int n1, int n2, int discardEdge, int nSNP, int nSamples, int sample, int downSample);

void writeData(double *Data, int *sampledData, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP, int downSample, char *outputName);

void writeCoordinates(int *coords, int npop, int n1, int n2, int discardEdge, int nSamples, char *outputName);

void prodMatrix(double *A, double *B, double *res, int nrowA, int ncolA, int nrowB, int ncolB);

void tr(double *Matrix, int nrow, int ncol);

void tAA(double *A, double *tAA, int nrow, int ncol);

void Covariance(double *A, double *Sigma, int nrow, int ncol, int transpose);

void colMeans(double *Matrix, int nrow, int ncol, double *Means);

void getCor(double **Cor, double *Data, double *OutsidePop, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP);

int which(int *sample1, int g1, int *sample2, int g2, int nSNP);

double Psi_ij(int *sample1, int *sample2, int nSNP);

void getPsi(double **Psi, int *sampledData, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP);

void writeCor(double *Cor, double *Psi, int *coords, int npop, int n1, int n2, int discardEdge, int nSamples, int calculatePsi, char *outputName);
