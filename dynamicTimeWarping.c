#include <float.h>
#include <math.h>
#include <stdint.h>
#include "dynamicTimeWarping.h"


/* ------------------------------------------------------------------------- *
 * This function computes a local cost measure (by evaluating the average
 * absolute distance) between two elements of two different signals.
 *
 * PARAMETERS
 *
 * signal1			First signal
 * signal2			Second signal
 * tStep1			The time step of signal1 which must be compared to compute
 *					the local cost measure
 * tStep2			The time step of signal2 which must be compared to compute
 *					the local cost measure
 * ------------------------------------------------------------------------- */
static double cost(Signal* signal1, Signal* signal2, size_t tStep1,
				   size_t tStep2){
	
	double absDist = 0;
	size_t n_coef = signal1->n_coef;
	
	for(size_t i=0; i < signal1->n_coef; i++)
		absDist += fabs(signal1->mfcc[i][tStep1] - signal2->mfcc[i][tStep2]);
	
	absDist = (double) absDist/n_coef;
	return absDist;
}

/* ------------------------------------------------------------------------- *
 * This function returns the minimum between three variables of type double :
 * a, b and c.
 * ------------------------------------------------------------------------- */
static double minVec(double a, double b, double c){
	
	double tmp = (a < b ? a : b);
	return (c < tmp ? c : tmp);
}

/* ------------------------------------------------------------------------- *
 * This function computes the dynamic time warping between two signals
 * subject to a locality constraint.
 * ------------------------------------------------------------------------- */
double dtw(Signal* signal1, Signal* signal2, size_t locality){
	
	size_t height;
	size_t width;
	Signal* shortSignal;
	Signal* longSignal;
	
	if(signal1->size < signal2->size){
		height = signal1->size;
		width = signal2->size;
		shortSignal = signal1;
		longSignal = signal2;
	}
	else{
		height = signal2->size;
		width = signal1->size;
		shortSignal = signal2;
		longSignal = signal1;
	}
	
	// Locality constraint prevents the dtw from computing the score
	if(locality < width - height)
		return DBL_MAX;
	
	// Initialization of the Accumulated Cost Matrix
	double ACMat[height + 1][width + 1];
	
	for(size_t i = 1; i <= height; i++)
		ACMat[i][0] = DBL_MAX;
	
	for(size_t j = 1; j <= width ; j++)
		ACMat[0][j] = DBL_MAX;
	
	ACMat[0][0] = cost(shortSignal, longSignal, 0, 0);
	
	// Computation of the dtw score
	for(size_t i = 1; i <= height; i++){
		// init = max(1, i - locality)
		size_t init = (locality >= i)? 1 : i - locality;
		
		// cond = min(width, i + locality)
		size_t cond = (locality != SIZE_MAX &&
					   width > i + locality)? i + locality : width;
		
		ACMat[i][init - 1] = DBL_MAX;
		ACMat[i][cond + 1] = DBL_MAX;
		
		for (size_t j = init; j <= cond; j++){
			
			ACMat[i][j] = cost(shortSignal, longSignal, i - 1, j - 1)
			+ minVec(ACMat[i-1][j], ACMat[i][j-1], ACMat[i-1][j-1]);
		}
	}
	
	return ACMat[height][width];
}
