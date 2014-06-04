/**
 * @file
 * @brief Sequence utilities
 *
 * This file contains utility functions for working with DNA sequences.
 */
#include <stdlib.h>
#include <string.h>
#include "sequence.h"

/**
 * @brief Frees the memory of a given sequence.
 * @param S - The sequence to free.
 */
void freeSeq( seq_t *S){
	free( S->S);
	free( S->RS);
	free( S->name);
	S->S = S->RS = S->name = NULL;
}


/**
 * Compute the reverse complement.
 * @param str The master string.
 * @param len The length of the master string
 * @return The reverse complement. The caller has to free it!
 */
char *revcomp( const char *str, size_t len){
	char *rev = (char*) malloc( len + 1);
	if( !str || !rev) return NULL;
	
	char *r = rev;
	const char *s = str + len-1;
	
	rev[len] = '\0';
	
	char c, d;
	
	while( len --> 0 ){
		c = *s--;
		
		switch( c){
			case 'A': d = 'T'; break;
			case 'T': d = 'A'; break;
			case 'G': d = 'C'; break;
			case 'C': d = 'G'; break;
			default: continue;
		}
		
		*r++ = d;
	}
	
	return rev;
}

/**
 * This function concatenates the reverse complement to a given master string. A
 * `#` sign is used as a separator.
 * @param s The master string.
 * @param len Its length.
 * @return The newly concatenated string.
 */
char *catcomp( char *s , size_t len){
	if( !s) return NULL;
	
	char *rev = revcomp( s, len);
	
	rev = (char*) realloc( rev, 2 * len + 2);
	if( !rev) return NULL;
	
	rev[len] = '#';
	
	memcpy( rev+len+1, s, len+1);
	
	return rev;
}

double calc_gc( seq_t *S){
	size_t GC = 0;
	char *p = S->S;
	
	for(; *p; p++){
		if( *p == 'G' || *p == 'C'){
			GC++;
		}
	}
	
	return S->gc = (double)GC/S->len;
} 

