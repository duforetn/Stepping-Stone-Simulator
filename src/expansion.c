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

/* Initialize all the demes with 0 allele frequencies for coordinates between popstart and popends */
void initExpansion2D(double *Pops, int nSNP, int n1, int n2, int popStart, int popEnd){

        int p1, p2, l;
	for (p1=0; p1<n1; p1++){
	        for (p2=popStart; p2<popEnd;  p2++){
	                for (l=0; l<nSNP; l++){
        	                Pops[(p1 + n1*p2)*nSNP + l] = .0;
                	}
		}
        }
}


/* This function treats the case where an expansion is currently happening, with isolation between populations involved in the process.*/
void expand1D(double *Pops, double *OutsidePop, double *Propagule, double *AF, int npop, int nSNP, double m1, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture){

        int npopColonized = (int) (gen - timeStart - 1)/lengthSteps + 1;
        int p, l;
        if ((gen == timeStart) && (!admixture)) initExpansion1D(Pops, nSNP, npop, popStart, popEnd);
        if (admixture && (gen == timeStart)) for (l=0; l<nSNP; l++) Propagule[l] = Pops[popStart*nSNP + l];
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
                if (((gen - timeStart)%lengthSteps == 0) && (npopColonized < (popEnd - popStart + 1))) expandNewpop1D(Pops, AF, Propagule, popStart + npopColonized - 1, nSNP, Ne, ke, admixture);
        }

        /* Treat the populations outside the expansion */
        if (popStart > 1) evolve1D(Pops, OutsidePop, AF, popStart - 1, nSNP, m1, minf, Ne, Loss);
        if (popEnd < npop) evolve1D(Pops + nSNP*popEnd, OutsidePop, AF, npop - popEnd, nSNP, m1, minf, Ne, Loss);
}

/* Perform the IBD and expansion during the range expansion time */
void expand2D(double *Pops, double *OutsidePop, double *AF, int n1, int n2, int nSNP, double m1, double m2, double minf, int Ne, int Loss, int gen, int timeStart, int timeEnd, int popStart, int popEnd, int nSteps, int lengthSteps, int ke, int admixture){

        int npopColonized = (int) (gen - timeStart - 1)/lengthSteps + 1;
        int p1, p2, l;
        if ((gen == timeStart) && (!admixture)) initExpansion2D(Pops, nSNP, n1, n2, popStart, popEnd);
        if (gen>timeStart){
		for (p1=0; p1<n1; p1++){
                for (p2=popStart - 1; p2<popStart - 1 + npopColonized; p2++){
                        for (l=0; l<nSNP; l++){
                                double af = Pops[(p1 + p2*n1)*nSNP + l];
                                if ((af > 0) && (af < 1)) af = (double) rand_normal(af, af*(1 - af)/(2*Ne));
                                if (af > 1) af = 1;
                                if (af < 0) af = 0;
                                Pops[(p1 + p2*n1)*nSNP + l] = af;
                        }
                }
                if (((gen - timeStart)%lengthSteps == 0) && (npopColonized < (popEnd - popStart + 1))) expandNewpop2D(Pops, AF, popStart + npopColonized - 1, nSNP, n1, n2, Ne, ke, admixture);
		}
        }

        /* Treat the populations outside the expansion */
        if (popStart > 1) evolve2D(Pops, OutsidePop, AF, (popStart - 1)*n1, n1, popStart - 1, nSNP, m1, m2, minf, Ne, Loss);
        if (popEnd < n2) evolve2D(Pops + nSNP*n1*popEnd, OutsidePop, AF, n1*(n2 - popEnd), n1, n2 - popEnd, nSNP, m1, m2, minf, Ne, Loss);


}

/* pop -1 is the last colonized deme, the new deme is pop*/
void expandNewpop1D(double *Pops, double *AF, double *Propagule, int pop, int nSNP, int Ne, int ke, int admixture){

        int l, i, count;
printf("Expand new pop: %i\n", pop + 1);
        for (l=0; l<nSNP; l++){
                count = 0;
                if (admixture){
                        for (i=0; i<ke; i++) if (drand() < Propagule[l]) count++;
                        Propagule[l] = (double) count/ke;
                        for (i=ke; i<Ne; i++) if (drand() < Pops[nSNP*(pop) + l]) count ++;
                        Pops[nSNP*pop + l] = (double) count/Ne;
                } else {
                        for (i=0; i<ke; i++) if (drand() < Pops[nSNP*(pop - 1) + l]) count++;
                        Pops[nSNP*pop + l] = (double) count/ke;
                }
        }

}

/* pop isthe coordinates on the second dimension of the newly colonized pop.*/
void expandNewpop2D(double *Pops, double *AF, int pop, int nSNP, int n1, int n2, int Ne, int ke, int admixture){

        int p1, l, i, count;
printf("Expand new pops in row: %i\n", pop + 1);
	for (p1=0; p1<n1; p1++){
        for (l=0; l<nSNP; l++){
                count = 0;
                for (i=0; i<ke; i++) if (drand() < Pops[nSNP*((pop - 1)*n1 + p1) + l]) count++;
                if (admixture){
                        for (i=ke; i<Ne; i++) if (drand() < Pops[nSNP*(pop*n1 + p1) + l]) count ++;
                }
                Pops[nSNP*(pop*n1 + p1) + l] = (double) count/ke;
                if (admixture) Pops[nSNP*(pop*n1 + p1) + l] = (double) count/Ne;
        }
	}

}

