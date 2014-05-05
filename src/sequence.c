#include <stdlib.h>
#include "sequence.h"

void freeSeq( seq_t *S){
	free( S->S);
	if( S->RS) free( S->RS);
	free( S->name);
}

