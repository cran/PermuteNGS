////////////////////////////////////////////////////////////////////
//compile : gcc distribution.c -odis
//execute : ./dis Human 
//input file : Human.count, Human_IPI_ID_chr_location_version_3.sort 
//output file : Human.distributon 
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define randomize() srand((unsigned)time(NULL));
#define random(num) (rand()%(num));

float z_test(int x, int y, int xsum, int ysum, int tail);
SEXP PermuteDEG(SEXP dataA, SEXP dataB, SEXP gene_count, SEXP samplingNumber, SEXP tail);

SEXP PermuteDEG(SEXP dataA, SEXP dataB, SEXP gene_count, SEXP samplingNumber, SEXP tail){
	int 	i,k,j;
	int	totdataA = 0,totdataB = 0;
	int	**BOOT_LIST;
	float	**Z_BOOT;
	int	*sumAB;
	int	*WD;
	int	WD_COUNT = 0;
//	float	*Z_MINE;
	int	gene_count2, samplingNumber2, tail2;

	int *PdataA, *PdataB;
	PROTECT( dataA= AS_INTEGER(dataA));
	PROTECT( dataB= AS_INTEGER(dataB));
	PdataA = INTEGER_POINTER(dataA);
	PdataB = INTEGER_POINTER(dataB);
	gene_count2 = INTEGER_VALUE(gene_count);
	samplingNumber2 = INTEGER_VALUE(samplingNumber);

	tail2 = INTEGER_VALUE(tail);

	SEXP PVALUE, Z_MINE;
	PROTECT( PVALUE = allocVector( REALSXP, gene_count2));
	PROTECT( Z_MINE = allocVector( REALSXP, gene_count2));
	double *Pval, *pZ_MINE;
	Pval = REAL(PVALUE);
	pZ_MINE  = REAL(Z_MINE);


/*
	SEXP BOOT_LIST;
	PROTECT( BOOT_LIST = allocMatrix( REALSXP, gene_count2,samplingNumber2));
	//int * BOOT_A;
	//BOOT_A = COMPLEX_POINTER(BOOT_LIST);	

*/
	/*
	for(i=0;i<20;i++){
		printf("%d\n",PdataA[i]);
		printf("%d\n",PdataB[i]);
	}
	*/

	BOOT_LIST = (int **)malloc(gene_count2 * sizeof(int *));
	Z_BOOT = (float **)malloc(gene_count2 * sizeof(float *));
	sumAB = (int *)malloc(gene_count2 * sizeof(int));
//	Z_MINE = (float *)malloc(gene_count2 * sizeof(float));
//	PVALUE = (float *)malloc(gene_count2 * sizeof(float));

	for(i=0;i<gene_count2;i++){
		sumAB[i] = PdataA[i]+PdataB[i];
		totdataA += PdataA[i];
		totdataB += PdataB[i];
	}

	WD = (int *)malloc((totdataA+totdataB+1000) * sizeof(int));
	for(i=0;i<gene_count2;i++){
		BOOT_LIST[i] = (int *)calloc(samplingNumber2, sizeof(int));
		for(k=0;k<samplingNumber2;k++){
			BOOT_LIST[i][k] = 0;
		}

		for(k=0;k<PdataA[i];k++){
			WD[WD_COUNT] = i;
			WD_COUNT++;
		}
		for(k=0;k<PdataB[i];k++){
			WD[WD_COUNT] = i;
			WD_COUNT++;
		}
		Z_BOOT[i] = (float *)malloc(samplingNumber2 * sizeof(float));
		for(k=0;k<samplingNumber2;k++){
			Z_BOOT[i][k] = 0.0;
		}
		pZ_MINE[i] = 0.0;
		Pval[i] = 0.0;
	}
	
	/////////////////////////////////////////////////////////////////////////
	//for( SAMPLING_COUNT_LIST(N) count ){
	//	for( 10000 ){
	//		N sampling in LIST_WHOLE(variable) => compute score
	//	}
	//}
	/////////////////////////////////////////////////////////////////////////
	int 	temp;
	int	*SELECT;
	int	cuindex;
	int	SAMPLING_CONT;

	SAMPLING_CONT = totdataA;

	randomize();

	//float Z = z_test(1,2,100,100);
	//printf("%f\n",Z);
	//exit(1);

	for(i=0;i<samplingNumber2;i++){
		if((i+1)%10 == 0 ){
			printf("random sampling without replace: %d\n",i+1);
		}
		SELECT = (int *)malloc(WD_COUNT *  sizeof(int));
		for(j=0;j<WD_COUNT;j++){
			SELECT[j] = 0;
		}
		for(j=0;j<SAMPLING_CONT;j++){
			temp = random(WD_COUNT);
			if( SELECT[temp] == 0 ){
				cuindex = WD[temp];
				BOOT_LIST[cuindex][i] += 1;
				SELECT[temp] = 1;
			} else {
				j--;
			}
		}
		free(SELECT);
	}


	int 	X,Y;
	int  	XSUM = totdataA,YSUM = totdataB;
	float 	Z;
	int	FP;

	for(i=0;i<gene_count2;i++){
		//mine
		X = PdataA[i];
		Y = PdataB[i];
		Z = z_test(X,Y,XSUM,YSUM, tail2);
		pZ_MINE[i] = Z;

		//boot 
		FP = 0;
		for(j=0;j<samplingNumber2;j++){
			X = BOOT_LIST[i][j];
			Y = sumAB[i] - X;
			Z = z_test(X,Y,XSUM,YSUM, tail2);
			Z_BOOT[i][j] = Z;
			if( pZ_MINE[i] < 0 ){
				if(Z_BOOT[i][j] <= pZ_MINE[i]) FP++;
			} else if( pZ_MINE[i] > 0 ) {
				if(Z_BOOT[i][j] >= pZ_MINE[i]) FP++;
			} else {
				continue;
			}


		}
		Pval[i] = FP/(samplingNumber2*1.0);
	}

/*
	SEXP list, list_name;
	char *names[1] = {"Pvalue"};
	//char *names[3] = {"Pvalue", "originalZ", "sampledA"};
	PROTECT(list_name = allocVector( STRSXP,1));

	for (i =0; i <1; i++){
		SET_STRING_ELT( list_name, i, mkChar(names[i])); 
	}

	PROTECT(list = allocVector( VECSXP,1));
	SET_VECTOR_ELT( list, 0, PVALUE);
//	SET_VECTOR_ELT( list, 1, Z_MINE);
//	SET_VECTOR_ELT( list, 2, BOOT_LIST);
	setAttrib(list, R_NamesSymbol, list_name);
*/
	free(sumAB);
	free(WD);
	for(i=0;i<gene_count2;i++){
		free(BOOT_LIST[i]);
		free(Z_BOOT[i]);
	}


	free(BOOT_LIST);
	free(Z_BOOT);
//	free(Z_MINE);
//	free(PVALUE);
	UNPROTECT(4);	
	return PVALUE;
}

float z_test(int x, int y, int xsum, int ysum, int tail2 ){
	if( x == 0 && y == 0 )
		return 1.0;

	float  newxsum = xsum*1.0;
	float  newysum = ysum*1.0;

	float  P1 = x/newxsum;
	float  P2 = y/newysum;
	float  z = 0;
	if ( tail2 == 1 ){
		z = (P1-P2);
	}else {
		if (P1 > P2)
			z = (P1-P2);
		else 
			z = (P2-P1);
	}
	return z;
}
