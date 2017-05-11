#include <float.h>
#include "predictDigit.h"
#include "dynamicTimeWarping.h"


/* ------------------------------------------------------------------------- *
 * This function predicts the digit represented by the signal thanks to the 
 * given database. It returns a DigitScore structure containing the predicted
 * digit and the dtw score.
 * ------------------------------------------------------------------------- */
DigitScore predictDigit(Signal* signal, Database* database, size_t locality) {

	DigitScore digitScore = {DBL_MAX, 0};
	
	for(size_t i = 0; i < 10; i++){
		LinkedList* curLList = database->samples[i];
		LLNode* curNode = curLList->head;
		
		while(curNode != NULL){
			double score = dtw(signal, (Signal*)curNode->value, locality);
			
			if(score < digitScore.score){
				digitScore.score = score;
				digitScore.digit = (int)i;
			}
			curNode = curNode->next;
		}
	}
	
	return digitScore;		
}
