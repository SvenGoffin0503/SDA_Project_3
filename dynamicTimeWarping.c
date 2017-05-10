#include <stddef.h>
#include <float.h>
#include <math.h>
#include "dynamicTimeWarping.h"
#include "Signal.h"


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

/* ------------------------------------------------------------------------- *
 * This function returns the minimum between the variables a and b of type 
 * size_t.
 * ------------------------------------------------------------------------- */
static size_t min(size_t a, size_t b){
	
	return ( a < b ? a : b );
}

/* ------------------------------------------------------------------------- *
 * This function returns the maximum between the variables a and b of type
 * size_t.
 * ------------------------------------------------------------------------- */
static size_t max(size_t a, size_t b){
	
	return ( a > b ? a : b );
}

/* ------------------------------------------------------------------------- *
 * This function returns the minimum between three variables of type double :
 * a, b and c.
 * ------------------------------------------------------------------------- */
static double minVec(double a, double b, double c){
	
	double tmp = (a < b ? a : b);
	return (c < tmp ? c : tmp);
}

/** ------------------------------------------------------------------------ *
 * This function computes the dynamic time warping between two signals 
 * subject to a locality constraint.
 * ------------------------------------------------------------------------- */
double dtw(Signal* signal1, Signal* signal2, size_t locality){
	
	size_t height = signal1->size;
	size_t width = signal2->size;
	double ACMat[height][width];
	
	// Locality constraint prevents the dtw from computing the score
	if(locality < abs((int)(height - width)))
		return DBL_MAX;
	
	for(size_t i = 0; i < height; i++)
		for(size_t j = 0; j < width ; j++)
			ACMat[i][j] = DBL_MAX;
	
	ACMat[0][0] = 0.0;
	
	// Computation of the dtw score
	for(size_t i = 1; i < height; i++){
		for (size_t j = max(1, i-locality); j < min(width, i + locality); j++){
			
			ACMat[i][j] = cost(signal1, signal2,i, j)
						+ minVec(ACMat[i-1][j], ACMat[i][j-1], ACMat[i-1][j-1]);
		}
	}
	
	return ACMat[height][width];
}
