#include "random.h"
#include "pops.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void InitPops(double **Pops, double **OutsidePop, double **AF, int npop, int nSNP){

	int i, j;
	*OutsidePop = malloc(sizeof(double)*nSNP);
	*Pops = malloc(sizeof(double)*nSNP*npop);
	*AF = malloc(sizeof(double)*npop);

	init_random();

	for (i=0; i<nSNP; i++){
		*(*(OutsidePop) + i) = (double) rand_beta(1, 9);
		//*(*(OutsidePop) + i) = (double) 100/1000;
	}

	for (j=0; j<npop; j++){
		for (i=0; i<nSNP; i++){
        		*(*(Pops) + j*nSNP + i) = *(*(OutsidePop) + i);
//			*(*(Pops) + j*nSNP + i) = rand_double(0, 1);
		}
	//	printf("%g\n", *(*(Pops) + j*nSNP));
        }

}

void evolve1D(double *Pops, double *OutsidePop, double *AF, int npop, int nSNP, double m1, double minf, int Ne, int Loss){

	int l, p, ret;
	double af, r_af;
        for (l=0; l<nSNP; l++){
		for (p=0; p<npop; p++){
			if (Loss && (Pops[p*nSNP + l] == .0 || Pops[p*nSNP + l] == 1.0)){
				ret = 0;
				if (AF[p] == 0){
					if(rand_double(0, 1) > minf*OutsidePop[l]) ret++;
					if (p > 0) if(rand_double(0, 1) > m1*Pops[p - 1]) ret++;
					if (p < npop) if(rand_double(0, 1) > m1*Pops[p + 1]) ret++;
					AF[p] = ret/Ne;
				} else {
                                        if(rand_double(0, 1) > minf*(1 - OutsidePop[l])) ret++;
                                        if (p > 0) if(rand_double(0, 1) > m1*(1 - Pops[p - 1])) ret++;
                                        if (p < npop) if(rand_double(0, 1) > m1*(1 - Pops[p + 1])) ret++;
					AF[p] = (Ne - ret)/Ne;
				}
			} else {
				if (p == 0){
					AF[p] = (1 - m1 - minf)*Pops[p*nSNP + l] + m1*(Pops[(p + 1)*nSNP + l]) + minf*OutsidePop[l];
				} else if (p == (npop - 1)){
        	        	        AF[p] = (1 - m1 - minf)*Pops[p*nSNP + l] + m1*(Pops[(p - 1)*nSNP + l]) + minf*OutsidePop[l];
				} else {
					AF[p] = (1 - m1 - minf)*Pops[p*nSNP + l] + m1*.5*(Pops[(p - 1)*nSNP + l] + Pops[(p + 1)*nSNP + l]) + minf*OutsidePop[l];
				}
			}
		}
		for (p=0; p<npop; p++){
			if(!Loss || ((AF[p] > .0) && (AF[p] < 1.0))){
				af = AF[p];
				r_af = (double) rand_normal(af, af*(1 - af)/(2*Ne));
				if (r_af > 1) r_af = 1;
				if (r_af < 0) r_af = 0;
				Pops[p*nSNP + l] = r_af;
			}
		}
	}

}

void evolve2D(double *Pops, double *OutsidePop, double *AF, int npop, int n1, int n2, int nSNP, double m1, double m2, double minf, int Ne, int Loss){

        int l, p1, p2, p;
        double af, r_af;
        for (l=0; l<nSNP; l++){
                for (p2=0; p2<n2; p2++){
			if (p2 == 0){

                               for(p1=0; p1<n1; p1++){
                                        if (p1 == 0){
//                                                AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 + 1)*nSNP] + l) + m2*(Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
						AF[0] = (1 - m1 - m2 - minf)*Pops[l] + m1*(Pops[(nSNP + l)]) + m2*(Pops[n1*nSNP + l]) + minf*OutsidePop[l];
                                        } else if (p1 == n1 - 1){
                                               AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l]) + m2*(Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
                                        } else {
                                               AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + .5*m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l] + Pops[(p1 + p2*n1 + 1)*nSNP + l]) + m2*(Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
                                        }
                                }


			} else if (p2 == n2 - 1) {

                               for(p1=0; p1<n1; p1++){
                                        if (p1 == 0){
                                                AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 + 1)*nSNP + l]) + m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
                                        } else if (p1 == n1 - 1){
                                               AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l]) + m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
                                        } else {
                                               AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + .5*m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l] + Pops[(p1 + p2*n1 + 1)*nSNP + l]) + m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
                                        }
                                }

			} else {

				for(p1=0; p1<n1; p1++){
					if (p1 == 0){
						AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 + 1)*nSNP + l]) + .5*m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l] + Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
					} else if (p1 == n1 - 1){
                                	       AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l]) + .5*m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l] + Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
					} else {
                	                       AF[p1 + p2*n1] = (1 - m1 - m2 - minf)*Pops[(p1 + p2*n1)*nSNP + l] + .5*m1*(Pops[(p1 + p2*n1 - 1)*nSNP + l] + Pops[(p1 + p2*n1 + 1)*nSNP + l]) + .5*m2*(Pops[(p1 + (p2 - 1)*n1)*nSNP + l] + Pops[(p1 + (p2 + 1)*n1)*nSNP + l]) + minf*OutsidePop[l];
					}
				}
			}
		}
                for (p=0; p<npop; p++){
                        if(!Loss || ((AF[p] > .0) && (AF[p] < 1.0))){
                                af = AF[p];
                                r_af = (double) rand_normal(af, af*(1 - af)/(2*Ne));
//printf("%g:%g ", af, r_af);
                                if (r_af > 1) r_af = 1;
                                if (r_af < 0) r_af = 0;
                                Pops[p*nSNP + l] = r_af;
                        }
                }
        }
}

int burninPop(double *Pops, double *OutsidePop, double *AF, int burnin, int npop, int n1, int n2, int m1, int m2, int minf, int nSNP, int Ne){

	int i;
	int Loss = 0;
        printf("Burnin...\n");
	for (i=0; i<burnin; i++){
                evolve1D(Pops, OutsidePop, AF, npop, nSNP, m1, minf, Ne, Loss);
	}
	return i;
}

