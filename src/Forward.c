#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "pops.h"
#include "random.h"
#include "matrix.h"

int main(int argc, char* argv[]){

	int Ne = 1000;
	double m1=.1, minf = .00004, m2 = 0;
	int n1 = 50, n2 = 1, npop, nSNP = 100, discardEdge = 0, Loss = 0, npopSampled;
	int i, j, burnin = 10, gen, ngen, nSamples = 1, nextSample, sample = 0, downSample = 0, calculatePsi = 0;
	int *Times, *coords, *sampledData;
	int rangeExpansion = 0, timeStart, admixture = 0, timeEnd, popStart, popEnd, nSteps, lengthSteps, ke;
	double *OutsidePop;
	char *outputName = malloc(sizeof(char)*128);
	double *Pops, *AF, *Data, *Cor, *Psi;

	if (argc == 1){
		Welcome(1);
		return 0;
	}

	if(handleParams(argc, argv, &npop, &n1, &n2, &discardEdge, &burnin, &nSNP, &m1, &m2, &minf, &Ne, &nSamples, &ngen, &downSample, &calculatePsi, &Loss, &nextSample, &Times, &rangeExpansion, &timeStart, &timeEnd, &popStart, &popEnd, &nSteps, &lengthSteps, &admixture, &ke, &outputName)) return 0;
	npopSampled = nSamples*(n1 - 2*discardEdge)*(n2 - 2*discardEdge);
	if (n2 == 1){
		Data = malloc(sizeof(double)*npopSampled*nSNP);
		if (downSample) sampledData = malloc(sizeof(int)*npopSampled*nSNP);
	} else {
		Data = malloc(sizeof(double)*npopSampled*nSNP);
		if (downSample) sampledData = malloc(sizeof(int)*npopSampled*nSNP);
	}

	InitPops(&Pops, &OutsidePop, &AF, npop, nSNP);
//for (i=0; i<n1 - 2*discardEdge; i++) printf("%g ", Pops[i*nSNP]);
//printf("\n");
//	if (burninPop(Pops, OutsidePop, AF, burnin, npop, n1, n2, m1, m2, minf, nSNP, Ne) != burnin) return 0;
	getCoordinates(&coords, Times, npop, n1, n2, discardEdge, nSamples);

        for (gen=1 - burnin; gen<ngen + 1; gen++){
		if (rangeExpansion && (gen > timeStart - 1) && (gen < timeEnd + 1)){
			if (n2 == 1) expand1D(Pops, OutsidePop, AF, npop, nSNP, m1, minf, Ne, Loss, gen, timeStart, timeEnd, popStart, popEnd, nSteps, lengthSteps, ke, admixture);
			if (n2 > 1) expand2D(Pops, OutsidePop, AF, n1, n2, nSNP, m1, m2, minf, Ne, Loss, gen, timeStart, timeEnd, popStart, popEnd, nSteps, lengthSteps, ke, admixture);
		} else {
			if (n2 == 1) evolve1D(Pops, OutsidePop, AF, npop, nSNP, m1, minf, Ne, Loss);
			if (n2 > 1) evolve2D(Pops, OutsidePop, AF, npop, n1, n2, nSNP, m1, m2, minf, Ne, Loss);
		}
		if (gen == nextSample){
                        printf("Sampling... %i\n", nextSample);
			nextSample = Times[findNext(Times, 1, nSamples, nextSample)];
			sample = sampleData(Data, sampledData, Pops, npop, n1, n2, discardEdge, nSNP, nSamples, sample, downSample);
		}
	}
	writeData(Data, sampledData, npop, n1, n2, discardEdge, nSamples, nSNP, downSample, outputName);
	writeCoordinates(coords, npop, n1, n2, discardEdge, nSamples, outputName);
	getCor(&Cor, Data, OutsidePop, npop, n1, n2, discardEdge, nSamples, nSNP);
	if (calculatePsi) getPsi(&Psi, sampledData, npop, n1, n2, discardEdge, nSamples, nSNP);
	writeCor(Cor, Psi, coords, npop, n1, n2, discardEdge, nSamples, calculatePsi, outputName);

	free(coords);
	free(Cor);
	free(AF);
	free(Pops);
	if (downSample) free(sampledData);
	if (calculatePsi) free(Psi);
	free(OutsidePop);
	free(Data);
	return 0;

}
