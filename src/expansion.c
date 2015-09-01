#include "random.h"
#include "pops.h"
#include "expansion.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void initExpansion1D(double *Pops, int nSNP, int npop, int popStart, int popEnd){

	int p, l;
	for (p=popStart; p<popEnd; p++){
		for (l=0; l<nSNP; l++){
			Pops[p*nSNP + l] = 0;
		}
	}

}

void initExpansion2D(double *Pops, int nSNP, int n1, int n2, int popStart, int popEnd){

        int p1, p2, l;
	int s1 = popStart%n1, s2 = (int) popStart/n1;
	int e1 = popEnd%n1, e2 = (int) popEnd/n1;
	for (p2=s2; p2<e2 + 1; p2++){
	        for (p1=s1; p1<e1 + 1; p1++){
	                for (l=0; l<nSNP; l++){
        	                Pops[(p1 + n1*p2)*nSNP + l] = .0;
                	}
		}
        }
/*for (p2=0; p2<n2; p2++){
	for (p1=0; p1<n1; p1++){
		for (l=0; l<nSNP; l++){
			printf("%g ", Pops[(p1 + n1*p2)*nSNP + l]);
		}
		printf("\n");
	}
}
*/
}


/* This function treats the case where an expansion is currently happening, with isolation between populations involved in the process.*/
void expand1D(double *Pops, double *OutsidePop, double *AF, int npop, int nSNP, double m1, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture){

	int npopColonized = (int) (gen - timeStart - 1)/lengthSteps + 1;
	int p, l;
	if ((gen == timeStart) && (!admixture)) initExpansion1D(Pops, nSNP, npop, popStart, popEnd);
        if (gen>timeStart){
		for (p=popStart - 1; p<popStart - 1 + npopColonized; p++){
			for (l=0; l<nSNP; l++){
				double af = Pops[p*nSNP + l];
				if ((af > 0) && (af < 1)) af = (double) rand_normal(af, af*(1 - af)/(2*Ne));
				if (af > 1) af = 1;
				if (af < 0) af = 0;
				Pops[p*nSNP + l] = af;
			}
		}
		if (((gen - timeStart)%lengthSteps == 0) && (npopColonized < (popEnd - popStart + 1))) expandNewpop1D(Pops, AF, popStart + npopColonized - 1, nSNP, Ne, ke, admixture);
	}

	/* Treat the populations outside the expansion */
	if (popStart > 1) evolve1D(Pops, OutsidePop, AF, popStart - 1, nSNP, m1, minf, Ne, Loss);
	if (popEnd < npop) evolve1D(Pops + nSNP*popEnd, OutsidePop, AF, npop - popEnd, nSNP, m1, minf, Ne, Loss);

}

void expand2D(double *Pops, double *OutsidePop, double *AF, int n1, int n2, int nSNP, double m1, double m2, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture){

/*        int npopColonized = (int) (gen - timeStart - 1)/lengthSteps + 1;
        int p, l;
        int s1 = popStart%n1, s2 = (int) popStart/n1;
        int e1 = popEnd%n1, e2 = (int) popEnd/n1;
p
        if ((gen == timeStart) && (!admixture)) initExpansion2D(Pops, nSNP, n1, n2, popStart, popEnd);
*/	
}

/* pop -1 is the last colonized deme, the new deme is pop*/
void expandNewpop1D(double *Pops, double *AF, int pop, int nSNP, int Ne, int ke, int admixture){

	int l, i, count;

	printf("Expand new pop: %i\n", pop + 1);
	for (l=0; l<nSNP; l++){
		count = 0;
		for (i=0; i<ke; i++) if (drand() < Pops[nSNP*(pop - 1) + l]) count++;
		if (admixture){
			for (i=ke; i<Ne; i++) if (drand() < Pops[nSNP*(pop) + l]) count ++;
		}
                Pops[nSNP*pop + l] = (double) count/ke;
		if (admixture) Pops[nSNP*pop + l] = (double) count/Ne;
	}

}
