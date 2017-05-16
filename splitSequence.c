#include <float.h>
#include "splitSequence.h"
#include "predictDigit.h"


/* ------------------------------------------------------------------------- *
 * Structure used to keep usefull information about a digit sequence.
 *
 * FIELDS
 * totScore			The lowest DTW score of the sequence
 * digit			The digit represented by the sequence
 * splitInd         The index of the beginning of the current sequence in the
 *					original sequence
 * ------------------------------------------------------------------------- */
typedef struct {
	double totScore;
	int digit;
	size_t splitInd;
}splittedSeq;


/* ------------------------------------------------------------------------- *
 * This function modies the variable signal in order to isolate one portion 
 * of the signal. It returns false if a problem occurs, true otherwise.
 *
 * PARAMETERS
 *
 * signal			The original signal
 * indStart			The index of the beginning of the signal portion that will 
 *					be isolated
 * indEnd			The index of the end of the signal portion that will be
 *					isolated
 * ------------------------------------------------------------------------- */
static bool cuttingSignal(Signal* signal, size_t indStart, size_t indEnd){
	
	if(indStart > indEnd || indStart >= signal->size || indEnd >= signal->size)
		return false;
	
	for (size_t i = 0; i < signal->n_coef; ++i) {
		signal->mfcc[i] = signal->mfcc[i] + indStart;
	}

	signal->size = indEnd - indStart + 1;
	return true;
}


/* ------------------------------------------------------------------------- *
 * This function reiitializes the original signal which has been modified by
 * the function cuttingSignal.
 *
 * PARAMETERS
 *
 * cutSignal			A valid pointer to the modified signal
 * originalsize			The size of the original signal
 * originalInd			A valid pointer to an array containing the start 
 *						indices of the original array
 * ------------------------------------------------------------------------- */
static void reinitSignal(Signal* cutSignal, size_t originalSize,
						 double** originalInd){
	
	cutSignal->size = originalSize;
	for(size_t i = 0; i < cutSignal->n_coef; i++)
		cutSignal->mfcc[i] = originalInd[i];
}


/* ------------------------------------------------------------------------- *
 * This function rebuilds the unknown sequence from results computed in the 
 * function bestSplit.
 *
 * PARAMETERS
 *
 * digitSeq			A valid pointer to the DigitSequence structure which will
 *					contain all information about the unknown sequence
 * splitSeq			A valid pointer to a splittedSeq structure which contains
 *					all the results computed by the function bestSplit
 * signal			A valid pointer to the unknown sequence
 * ------------------------------------------------------------------------- */
static void retrieveSeq(DigitSequence* digitSeq, splittedSeq* splitSeq,
								 Signal* signal){
	
	size_t k = 0;
	size_t i = signal->size - 1;
	
	// Computes the length of the unknown sequence
	while(i < signal->size){
		i = splitSeq[i].splitInd - 1;
		k++;
	}
	
	// Allocates necessary memory
	digitSeq->digits = malloc(sizeof(int) * k);
	if(!digitSeq->digits)
		return;
	
	digitSeq->splits = malloc(sizeof(size_t) * k);
	if(!digitSeq->splits){
		free(digitSeq->digits);
		return;
	}
	
	// Reconstitution of the unknown sequence
	digitSeq->nDigits = k;
	i = signal->size - 1;
	
	while(i < signal->size){
		k--;
		digitSeq->digits[k] = splitSeq[i].digit;
		digitSeq->splits[k] = splitSeq[i].splitInd;
		i = splitSeq[i].splitInd - 1;
	}
	digitSeq->score = splitSeq[signal->size - 1].totScore;
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
	size_t originalSize = signal->size;
	double* originalInd[signal->n_coef];
	
	for(size_t i = 0; i < signal->n_coef; i++)
		originalInd[i] = signal->mfcc[i];
	
	for(size_t i = lMin - 1; i < signal->size; i++){
		DigitScore curDigit;
		splitSeq[i].totScore = DBL_MAX;
		
		// The current signal portion can't be divided into several sequences
		if(i < lMax){
			if(!cuttingSignal(signal, 0, i))
				return digitSeq;
			
			curDigit = predictDigit(signal, database, locality);
			splitSeq[i].totScore = curDigit.score;
			splitSeq[i].digit = curDigit.digit;
			splitSeq[i].splitInd = 0;
			reinitSignal(signal, originalSize, originalInd);
		}
		
		// The current signal portion can be divided into several sequences
		if(i >= (2 * lMin) - 1){
			for(size_t j = lMin - 1; j < lMax; j++){
				if(i - j < lMin)
					break;
				
				if(!cuttingSignal(signal, i - j, i))
					return digitSeq;
		
				curDigit = predictDigit(signal, database, locality);
				curDigit.score += splitSeq[i - j - 1].totScore;
				reinitSignal(signal, originalSize, originalInd);
				
				// Maintains the computed split with the lowest score
				if(curDigit.score < splitSeq[i].totScore){
					splitSeq[i].totScore = curDigit.score;
					splitSeq[i].digit = curDigit.digit;
					splitSeq[i].splitInd = i - j;
				}
			}
		}
	}
	
	// Rebuilds the unknown sequence
	retrieveSeq(&digitSeq, splitSeq, signal);
	return digitSeq;
}
