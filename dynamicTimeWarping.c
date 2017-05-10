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

static size_t min(size_t a, size_t b){
    size_t tmp = DBL_MAX;
    if(a < b)
        tmp = a;
    else
        tmp = b;
    return tmp;
}

static size_t max(size_t a, size_t b){
    size_t tmp = DBL_MAX;
    if(a > b)
        tmp = a;
    else
        tmp = b;
    return tmp;
}


static double minVec(double a, double b, double c){
    double tmp = DBL_MAX;
    if(a < tmp)
        tmp = a;
    if(b < tmp)
        tmp = b;
    if(c < tmp)
        tmp = c;
    
    return tmp;
}

double dtw(Signal* signal1, Signal* signal2, size_t locality){
    size_t height = signal1->size;
    size_t width = signal2->size;
    double ACMat[height][width];
    
    if(locality > abs(height - width))
        return DBL_MAX;
    
    for(size_t i = 0; i < height; i++)
        for(size_t j = 0; j < width ; j++)
            ACMat[i][j] = DBL_MAX;
    
    ACMat[0][0] = 0;
    
    for(size_t i = 0; i < height; i++){
        for (size_t j = max(1, i-locality); j < min(width, i + locality); j++) {
            ACMat[i][j] = cost(signal1, signal2,i, j) + minVec(ACMat[i-1][j], ACMat[i][j-1], ACMat[i-1][j-1]);
        }
    }
    
    return ACMat[height][width];
}
