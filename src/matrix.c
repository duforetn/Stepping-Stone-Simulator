#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src_Lapack/lpk.h"
#include "matrix.h"

#define NA 9

int findMax(int *matrix, int nrow, int ncol){

	int i, imax = 0;
	for (i=0; i<ncol*nrow; i++){
		if(matrix[imax] < matrix[i]){
			imax = i;
		}
	}
	return imax;
}

int findNext(int *matrix, int nrow, int ncol, int current){

	int i, inext = findMax(matrix, nrow, ncol);
	for (i=0; i<nrow*ncol; i++){
		if ((matrix[i] < matrix[inext]) && (matrix[i] > current)) inext = i;
	}
	return inext;
}

int sampleData(double *Data, int *sampledData, double *Pops, int npop, int n1, int n2, int discardEdge, int nSNP, int nSamples, int sample, int downSample){
	
	int l, p, g, p1, p2;
	int npopSampledGeneration = (n1 - 2*discardEdge)*(n2 - 2*discardEdge);
	/* 
	contiguously writing all populations for one time sample.
	*/
	if (n2 == 1){
		npopSampledGeneration = (n1 - 2*discardEdge);
		for (p=discardEdge; p<npop - discardEdge; p++){
			for (l=0; l<nSNP; l++){
				Data[(p - discardEdge)*nSNP + l + sample*(npop - 2*discardEdge)*nSNP] = Pops[p*nSNP + l];
				if (downSample){
					g = 0;
					if(rand_double(0.0, 1.0) < Data[(p - discardEdge)*nSNP + l + sample*(npop - 2*discardEdge)*nSNP]) g++;
					if(rand_double(0.0, 1.0) < Data[(p - discardEdge)*nSNP + l + sample*(npop - 2*discardEdge)*nSNP]) g++;
					sampledData[(p - discardEdge)*nSNP + l + sample*(npop - 2*discardEdge)*nSNP] = g;
				}
			}
		}
	} else {
                for (p2=discardEdge; p2<n2 - discardEdge; p2++){
			for (p1 = discardEdge; p1<n1 - discardEdge; p1++){
	                        for (l=0; l<nSNP; l++){
                        	        Data[((p2 - discardEdge)*(n1 - 2*discardEdge) + p1 - discardEdge)*nSNP + l + sample*npopSampledGeneration*nSNP] = Pops[(p2*n1 + p1)*nSNP + l];
                	                if (downSample){
        	                                g = 0;
	                                        if(rand_double(0.0, 1.0) < Data[((p2 - discardEdge)*n1 + p1)*nSNP + l + sample*npopSampledGeneration*nSNP]) g++;
                                        	if(rand_double(0.0, 1.0) < Data[((p2 - discardEdge)*n1 + p1)*nSNP + l + sample*npopSampledGeneration*nSNP]) g++;
                                	        sampledData[((p2 - discardEdge)*n1 + p1)*nSNP + l + sample*npopSampledGeneration*nSNP] = g;
                        	        }
                 	       }
			}
                }
	}
	return sample + 1;	
}

void writeData(double *Data, int *sampledData, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP, int downSample, char *outputName){

	FILE *Output, *Output2;
	int p, l, s;
        char *dataName = calloc(256, sizeof(char));
	char *sampleName = calloc(256, sizeof(char));
	int npopSampled = nSamples*(n1 - 2*discardEdge)*(n2 - 2*discardEdge);
	if (n2 == 1) npopSampled = nSamples*(n1 - 2*discardEdge);

        strcpy(dataName, outputName);
        strcat(dataName, ".dat");

	if((Output = fopen(dataName, "w")) == NULL) printf("ERROR, unable to open %s\n", dataName);

	strcpy(sampleName, outputName);
        strcat(sampleName, ".samples");

        if((Output2 = fopen(sampleName, "w")) == NULL) printf("ERROR, unable to open %s\n", sampleName);

	for (p=0; p<npopSampled; p++){
		for (l=0; l<nSNP; l++){
			fprintf(Output, "%g ", Data[p*nSNP + l]);
			if (downSample) fprintf(Output2, "%i ", sampledData[p*nSNP + l]);
		}
		if (downSample) fprintf(Output2, "\n");
		fprintf(Output, "\n");
	}
	fclose(Output);
	fclose(Output2);
}

void writeCoordinates(int *coords, int npop, int n1, int n2, int discardEdge, int nSamples, char *outputName){

        FILE *Output;
        int p1, p2, s;
        char *coordsName = calloc(256, sizeof(char));

        strcpy(coordsName, outputName);
        strcat(coordsName, ".coords");

        if((Output = fopen(coordsName, "w")) == NULL) printf("ERROR, unable to open %s\n", coordsName);

	if (n2 > 1){
        	int npopSampled = nSamples*(n1 - 2*discardEdge)*(n2 - 2*discardEdge);
	        for (s=0; s<npopSampled; s++){
        	        fprintf(Output, "%i %i %i\n", coords[3*s], coords[1 + 3*s], coords[s*3 + 2]);
	        }
	} else {
		for (s=0; s<nSamples; s++){
	                for (p1=discardEdge; p1<n1 - discardEdge; p1++){
                              fprintf(Output, "%i 1 %i\n", coords[s*3*(n1 - 2*discardEdge) + (p1 - discardEdge)*3],  coords[s*3*(n1 - 2*discardEdge) + (p1 - discardEdge)*3 + 2]);
                        }
                }
	}

	fclose(Output);
}

void prodMatrix(double *A, double *B, double *res, int nrowA, int ncolA, int nrowB, int ncolB){

                long int nrtA = nrowA, nctB = ncolB, nctA = ncolA, nrtB = nrowB;
                long int ldc = nctB, lda = nctA, ldb = nctB;
                double alpha = 1, beta = 0;
                int r = dgemmlig_("T", "T", &nrtA, &nctB, &nctA, &alpha, A, &lda, B, &ldb, &beta, res, &ldc);

}

void tr(double *Matrix, int nrow, int ncol){

        int i, j;
        double *tmp = calloc(nrow*ncol, sizeof(double));
        for (i=0; i<ncol*nrow; i++) tmp[i] = Matrix[i];
        for (i=0; i<nrow; i++){
                for (j=0; j<ncol; j++){
                        Matrix[j*nrow+ i] = tmp[i*ncol + j];
                }
        }
        free(tmp);
}

void tAA(double *A, double *tAA, int nrow, int ncol){

		int i;
                double *tA = malloc(sizeof(double)*nrow*ncol);
                for (i=0; i<nrow*ncol; i++) tA[i] = A[i];
                long int nrtA = ncol, nctB = ncol, nctA = nrow, nrtB = nrow;
                long int ldc = nrtA, lda = nrtA, ldb = nrtB;
                double alpha = 1, beta = 0;
                tr(A, nrow, ncol);
                int r = dgemm_("N", "N", &nrtA, &nctB, &nrtB, &alpha, tA, &lda, A, &ldb, &beta, tAA, &ldc);
                tr(A, ncol, nrow);
                free(tA);
	
}

void Covariance(double *A, double *Sigma, int nrow, int ncol, int transpose){

        int i;
        if (transpose){
                tAA(A, Sigma, nrow, ncol);
        } else {
                double *tA = calloc(nrow*ncol, sizeof(double));
		for (i=0; i<nrow*ncol; i++) tA[i] = A[i];
                tr(tA, nrow, ncol);
                prodMatrix(A, tA, Sigma, nrow, ncol, ncol, nrow);
                free(tA);
        }

}

void colMeans(double *Matrix, int nrow, int ncol, double *Means){

        int row, col, na;
        for (col=0; col<ncol; col++){
                Means[col] = 0;
                na = 0;
                for (row=0; row<nrow; row++){
                        if (Matrix[row*ncol + col] != NA){
                                Means[col]+=Matrix[row*ncol + col];
                        } else {
                                na++;
                        }
                }
                Means[col] = Means[col]/(nrow - na);
        }

}

void getCor(double **Cor, double *Data, double *OutsidePop, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP){

	int s, i, j, nrow;
	double *varS;
	double *Means = malloc(sizeof(double)*nSNP);
	if (n2 > 1){
		nrow = (n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples;
		*Cor = malloc(sizeof(double)*(n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples*(n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples);
	} else {
                nrow = (n1 - 2*discardEdge)*nSamples;
                *Cor = malloc(sizeof(double)*nrow*nrow);
	}

	for (s=0; s<nSamples; s++) {
		colMeans(Data + (nrow/nSamples)*s*nSNP, nrow/nSamples, nSNP, Means);
		for (i=0; i<nrow/nSamples; i++){
			for (j=0; j<nSNP; j++){
				Data[i*nSNP + j] -= (Means[j] - OutsidePop[j]);
				Data[s*nSNP*(nrow/nSamples) + i*nSNP + j] -= OutsidePop[j];
			}
		}
	}
	Covariance(Data, *Cor, nrow, nSNP, 0);

	varS = malloc(sizeof(double)*nrow);
	for (i=0; i<nrow; i++) varS[i] = *(*(Cor) + i*nrow + i);

	for (i=0; i<nrow; i++){
		for (j=0; j<nrow; j++){
			*(*(Cor) + i*nrow + j) /= sqrt(varS[i]*varS[j]);
		}
	}
	free(varS);

}

int which(int *sample1, int g1, int *sample2, int g2, int nSNP){

	int i, c = 0;
	for (i=0; i<nSNP; i++){
		if ((*(sample1 + i) == g1) && (*(sample2 + i) == g2)) c++;
	}
	
	return c;

}

double Psi_ij(int *sample1, int *sample2, int nSNP){

	int f12 = 0, f21 = 0, f11 = 0;
	double psi;

	f11 = which(sample1, 1, sample2, 1, nSNP);
        f12 = which(sample1, 1, sample2, 2, nSNP);
        f21 = which(sample1, 2, sample2, 1, nSNP);

	psi = (double) (f21 - f12)/(f12 + f11 + f21);
	if (isnan(psi)) return 0;
	return psi;

}

void getPsi(double **Psi, int *sampledData, int npop, int n1, int n2, int discardEdge, int nSamples, int nSNP){

        int i, j, nrow;
	double p;
        if (n2 > 1){
                nrow = (n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples;
                *Psi = calloc((size_t) (n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples*(n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples, sizeof(double));
        } else {
                nrow = (n1 - 2*discardEdge)*nSamples;
                *Psi = calloc((size_t) nrow*nrow, sizeof(double));
        }

	for (i=0; i<(nrow - 1); i++){
		for (j=(i + 1); j<nrow; j++){
			p = Psi_ij(sampledData + i*nSNP, sampledData + j*nSNP, nSNP);
			*(*(Psi) + i*nrow + j) = p;
                        *(*(Psi) + j*nrow + i) = -p;
		}
	}

}

void writeCor(double *Cor, double *Psi, int *coords, int npop, int n1, int n2, int discardEdge, int nSamples, int calculatePsi, char *outputName){

        FILE *Output;
        int i, j, nrow;
        char *corName = calloc(256, sizeof(char));

        strcpy(corName, outputName);
        strcat(corName, ".r");


        if((Output = fopen(corName, "w")) == NULL) printf("ERROR, unable to open %s\n", corName);

        if (n2 > 1){
                nrow = (n1 - 2*discardEdge)*(n2 - 2*discardEdge)*nSamples;
        } else {
                nrow = (n1 - 2*discardEdge)*nSamples;
        }

	fprintf(Output, "n1 n2 Time n1 n2 Time r");
	if (calculatePsi) fprintf(Output, " Psi");
	fprintf(Output, "\n");
	for (i=0; i<nrow - 1; i++){
		for (j=i + 1; j<nrow; j++){
			fprintf(Output, "%i %i %i %i %i %i %g", coords[i*3], coords[i*3 + 1], coords[i*3 + 2], coords[j*3], coords[j*3 + 1], coords[j*3 + 2], Cor[i*nrow + j]);
		        if (calculatePsi) fprintf(Output, " %g", Psi[i*nrow + j]);
			fprintf(Output, "\n");
		}
	}

	fclose(Output);
}

