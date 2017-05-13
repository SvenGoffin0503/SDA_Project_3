#include "splitSequence.h"
#include "predictDigit.h"
#include <float.h>
#include <stdio.h>

typedef struct {
	double totScore;
	int digit;
	size_t splitInd;
}splittedSeq;



static Signal* cuttingSignal(Signal* signal, size_t indStart, size_t indEnd){
	
	if(indStart > indEnd)
		return NULL;
	
	Signal* cutSignal = (Signal*) malloc(sizeof(Signal));
	if(!cutSignal)
		return NULL;
	
	cutSignal->n_coef = signal->n_coef;
	cutSignal->size = indEnd - indStart + 1;
	
	// allocate matrix
	cutSignal->mfcc = malloc(sizeof(double*) * cutSignal->n_coef);
	if (!cutSignal->mfcc) {
		free(signal);
		return NULL;
	}

	for (size_t i = 0; i < cutSignal->n_coef; ++i) {
		cutSignal->mfcc[i] = malloc(sizeof(double) * cutSignal->size);
		if (!signal->mfcc[i]) {
			// deallocate previously allocated arrays
			for(size_t j = i - 1; j < i; --j) {
				free(signal->mfcc[j]);
			}
			free(signal->mfcc);
			free(signal);
			return NULL;
		}
	}
	
	// fill matrix
	for (size_t i = 0; i < cutSignal->n_coef; ++i) {
		for (size_t j = 0; j < cutSignal->size; ++j) {
			cutSignal->mfcc[i][j] = signal->mfcc[i][indStart + j];
		}
	}

	return cutSignal;
}


/* ------------------------------------------------------------------------- *
 * Given a signal containing a sequence of digits, this function finds the best
 * split of the signal isolating each digit in the sequence.
 *
 * Returns the optimal digit sequence, its score, its corresponding split
 * indexes as a structure.
 * ------------------------------------------------------------------------- */
DigitSequence bestSplit(Signal* signal, Database* database,
						size_t locality, size_t lMin, size_t lMax){
	
	DigitSequence digitSeq = {0, 0, NULL, NULL};
	splittedSeq splitSeq[signal->size];
	
	for(size_t i = 0; i < signal->size; i++){
		splitSeq[i].totScore = DBL_MAX;
		splitSeq[i].splitInd = 0;
		splitSeq[i].digit = -1;
				printf("%d\n", i);
	}

	for(size_t i = lMin - 1; i < signal->size; i++){

		DigitScore curDigit;
		
		if(i < lMax){
			Signal* cutSignal = cuttingSignal(signal, 0, i);
			if(!cutSignal){
				return digitSeq;
			}

			curDigit = predictDigit(cutSignal, database, locality);
			splitSeq[i].totScore = curDigit.score;
			splitSeq[i].digit = curDigit.digit;
			splitSeq[i].splitInd = 0;
			freeSignal(cutSignal);
		}
		
		if(i >= (2 * lMin) - 1){
			for(size_t j = lMin - 1; j < lMax; j++){
				if(i - j == 0)
					break;
				Signal* cutSignal = cuttingSignal(signal, i - j, i);
				if(!cutSignal){
					return digitSeq;
				}
				curDigit = predictDigit(cutSignal, database, locality);
				curDigit.score += splitSeq[i - j - 1].totScore;
				freeSignal(cutSignal);
				
				if(curDigit.score < splitSeq[i].totScore){
					splitSeq[i].totScore = curDigit.score;
					splitSeq[i].digit = curDigit.digit;
					splitSeq[i].splitInd = i - j;
				}
			}
		}
	}
	
	size_t k = 0;
	size_t i = signal->size - 1;
	
	while(i < signal->size){
		i = splitSeq[i].splitInd - 1;
		k++;
	}
	
	digitSeq.digits = malloc(sizeof(int) * k);
	if(!digitSeq.digits)
		return digitSeq;
	
	digitSeq.splits = malloc(sizeof(size_t) * k);
	if(!digitSeq.splits){
		free(digitSeq.digits);
		return digitSeq;
	}
	
	digitSeq.nDigits = k;
	i = signal->size - 1;
	
	while(i < signal->size){
		k--;
		digitSeq.digits[k] = splitSeq[i].digit;
		digitSeq.splits[k] = splitSeq[i].splitInd;
		i = splitSeq[i].splitInd - 1;
	}
	digitSeq.score = splitSeq[signal->size - 1].totScore;
	
	return digitSeq;
}
