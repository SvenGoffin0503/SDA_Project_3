#include "dynamicTimeWarping.h"
#include "Signal.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>


/* ------------------------------------------------------------------------- *
 * This function computes a local cost measure (by evaluating the average
 * absolute distance) between two elements of two different signals.
 *
 * PARAMETERS
 *
 * elSignalA		Element of the signal A
 * elSignalB		Element of the signal B
 * n_coef			Number of coefficient used in the MFCC transformation (it
 *					corresponds to the length of elSignalA and elSignalB)
 * ------------------------------------------------------------------------- */
static int cost (Signal* signA, Signal* signB, size_t tStepA, size_t tStepB){
	
	double absDist = 0;
	size_t n_coef = signA->n_coef;
	
	for(size_t i=0; i < signA->n_coef; i++)
		absDist += fabs(signA->mfcc[i][tStepA] - signB->mfcc[i][tStepB]);
	
	absDist = (double) absDist/n_coef;
	return absDist;
}

static double min(double a, double b, double c, double d, double e){
	double tmp = DBL_MAX;
	if(a < tmp)
		tmp = a;
	if(b < tmp)
		tmp = b;
	if(c < tmp)
		tmp = c;
	if(d < tmp)
		tmp = d;
	if(e < tmp)
		tmp = e;
	
	return tmp;
}

static double** accuCostMat (Signal* signA, Signal* signB){
	
	size_t height = signA->size+3;
	size_t width = signB->size+3;
	double ACMat[height][width];
	
	
	// Fills the first three columns and lines with infinite value
	for(size_t j=0; j < width; j++){
		ACMat[0][j] = DBL_MAX;
		ACMat[1][j] = DBL_MAX;
		ACMat[2][j] = DBL_MAX;
	}
	
	for(size_t i=3; i < height; i++){
		ACMat[i][0] = DBL_MAX;
		ACMat[i][1] = DBL_MAX;
		ACMat[i][2] = DBL_MAX;
	}
	
	ACMat[3][3] = cost(signA, signB, 3, 3);
	
	for(size_t j=4; j < width; j++)
		ACMat[3][j] = ACMat[3][j-1] + cost(signA, signB, 3, j);
	
	for(size_t i=4; i < height; i++)
		ACMat[i][3] = ACMat[i-1][3] + cost(signA, signB, i, 3);
		
	for(size_t i=4; i < height; i++){
		for(size_t j=4; j < width; j++){
			double a = ACMat[i-1][j-1] + cost(signA, signB, i, j);
			double b = ACMat[i-2][j-1] + cost(signA, signB, i-1, j)
						+ cost(signA, signB, i, j);
			double c = ACMat[i-1][j-2] + cost(signA, signB, i, j-1)
						+ cost(signA, signB, i, j);
			double d = ACMat[i-3][j-1] + cost(signA, signB, i-1, j)
						+ cost(signA, signB, i-2, j) + cost(signA, signB, i, j);
			double e = ACMat[i-1][j-3] + cost(signA, signB, i, j-1)
						+ cost(signA, signB, i, j-2) + cost(signA, signB, i, j);
			ACMa min(a, b, c, d, e);
		}
	}
}
