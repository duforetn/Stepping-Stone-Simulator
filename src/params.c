#include "params.h"
#include "matrix.h"

void Welcome(int help){
        printf("\t\t/******************************\\\n");
        printf("\t\t|*** Forward-Stepping-stone ***|\n");
        printf("\t\t\\******************************/\n\n");
        printf("\n\n");

        if (help) printf("Command line Help:\n\t-o output files\n\t-m m1 m2 minf -n npop1 npop2 Ne\n\t-l number of loci\n\t-d number of demes to discard on each side\n\t-t nTimeSamples t1 ... tn\n\t-E 1/0 (admixture or not) popStart popEnd timeStart timeEnd ke\n");
}


int handleParams(int argc, char* argv[], int *npop, int *n1, int *n2, int *discardEdge, int *burnin, int *nSNP, double *m1, double *m2, double *minf, int *Ne, int *nSamples, int *ngen, int *downSample, int *calculatePsi, int *Loss, int *nextSample, int **Times, int *rangeExpansion, int *timeStart, int *timeEnd, int *popStart, int *popEnd, int *nSteps, int *lengthSteps, int *admixture, int *ke, char **outputName){

        int i, ii, tmp;
        for (i=0; i<argc; i++){
                if (argv[i][0] == '-'){
                        switch(argv[i][1]){
                                case 'n':
					*n1 = atoi(argv[i + 1]);
					*n2 = atoi(argv[i + 2]);
					*Ne = atoi(argv[i + 3]);
					*npop = (*n1) * (*n2);
					printf("%i pops of effective size %i\n", (*n1) * (*n2), *Ne);
				break;
				case 'm':
					*m1 = (double) atof(argv[i + 1]);
					*m2 = (double) atof(argv[i + 2]);
                                        *minf = (double) atof(argv[i + 3]);
				break;
				case 'l':
					*nSNP = atoi(argv[i + 1]);
					printf("%i Snps in the data.\n", *nSNP);
				break;
				case 'd':
					*discardEdge = atoi(argv[i + 1]);
				break;
                                case 'b':
                                        *burnin = atoi(argv[i + 1]);
                                break;
				case 't': 
					*nSamples = atoi(argv[i + 1]);
					*Times = (int *) malloc((*nSamples)*sizeof(int));
					for (ii=0; ii<*nSamples; ii++) *(*(Times) + ii) = atoi(argv[i + 2 + ii]);
					printf("Sampling at generation: ");
					for (ii=0; ii<*nSamples; ii++) printf("%i ", *(*(Times) + ii));
					printf("\n");
				break;
				case 'o':
					*outputName = strdup(argv[i + 1]);
				break;
				case 's':
					*downSample = atoi(argv[i + 1]);
				break;
				case 'p':
					*downSample = atoi(argv[i + 1]);
					*calculatePsi = atoi(argv[i + 1]);
				break;
				case 'L':
					*Loss = atoi(argv[i + 1]);
				break;
				case 'E':
					*admixture = atoi(argv[i + 1]);
					*popStart = atoi(argv[i + 2]);
					*popEnd = atoi(argv[i + 3]);
					*timeStart = atoi(argv[i + 4]);
                                        *timeEnd = atoi(argv[i + 5]);
					*ke = atoi(argv[i + 6]);
					*nSteps = (*popEnd) - (*popStart);
printf("nsteps %i\n", *nSteps);
					*lengthSteps = (int) (*timeEnd - *timeStart)/(*nSteps);
					if (*nSteps == ((int) (*timeEnd - *timeStart) + 1)) *lengthSteps = 1;
					if ((*lengthSteps > 0) && (*ke < *Ne) && (*nSteps > 0)){ *rangeExpansion = 1; printf("Range Expansion between population %i and %i\n\t-Starting at time %i\n\t-Ending at time %i\n\t-Founder effect %i/%i\n\t-%i steps of length %i\n", *popStart, *popEnd, *timeStart, *timeEnd, *ke, *Ne, *nSteps, *lengthSteps);}
				break;
			}
		}
	}

        *ngen = *(*(Times) + findMax(*Times, 1, *nSamples));
        *nextSample = *(*(Times) + findNext(*Times, 1, *nSamples, 0));
	if ((*discardEdge)*2 >= *npop){
		printf("incompatible population parameters\n");
		return 1;
	}

        if(!strcmp(*outputName, "")) *outputName = "output";
        printf("results in files %s*\n", *outputName);

	return 0;
}

void getCoordinates(int **coords, int *Times, int npop, int n1, int n2, int discardEdge, int nSamples){

	int p1, p2, s, nextSample;

	if (n2 == 1) *coords = malloc(sizeof(int)*3*(n1 - 2*discardEdge)*nSamples);
	else *coords = malloc(sizeof(int)*3*nSamples*((npop - 2*discardEdge)*n2 - 2*(n1 - 2*discardEdge)*discardEdge));

	nextSample = Times[findNext(Times, 1, nSamples, 0)];
	if (n2 > 1){
		for (s=0; s<nSamples; s++){
			for (p1=discardEdge; p1<n1 - discardEdge; p1++){
				for (p2=discardEdge; p2<n2 - discardEdge; p2++){
					*(*(coords) + 3*s*(n1 - 2*discardEdge)*(n2 - 2*discardEdge) + (p2 - discardEdge)*(n1 - 2*discardEdge)*3 + (p1 - discardEdge)*3) = p1 + 1;
					*(*(coords) + 3*s*(n1 - 2*discardEdge)*(n2 - 2*discardEdge) + (p2 - discardEdge)*(n1 - 2*discardEdge)*3 + (p1 - discardEdge)*3 + 1) = p2 + 1;
					*(*(coords) + 3*s*(n1 - 2*discardEdge)*(n2 - 2*discardEdge) + (p2 - discardEdge)*(n1 - 2*discardEdge)*3 + (p1 - discardEdge)*3 + 2) = nextSample;;
				}
			}
        	        nextSample = Times[findNext(Times, 1, nSamples, nextSample)];
		}
	} else {
		for (s=0; s<nSamples; s++){
	                for (p1=discardEdge; p1<n1 - discardEdge; p1++){
                                *(*(coords) + 3*s*(n1 - 2*discardEdge) + (p1 - discardEdge)*3) = p1 + 1;
                                *(*(coords) + 3*s*(n1 - 2*discardEdge) + (p1 - discardEdge)*3 + 1) = 1;
                                *(*(coords) + 3*s*(n1 - 2*discardEdge) + (p1 - discardEdge)*3 + 2) = nextSample;
                	}
	                nextSample = Times[findNext(Times, 1, nSamples, nextSample)];
	        }
	}

}
