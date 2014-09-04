/**
 * @file
 * @brief Sequence utilities
 *
 * This file contains utility functions for working with DNA sequences.
 */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include "sequence.h"
#include "global.h"

void normalize( seq_t *S);

using namespace std;

/**
 * @brief Convert an array of multiple sequences into a single sequence.
 * 
 * This function joins all sequences containd in an array into one
 * long sequence. The sequences are seperated by a `!` character. The
 * caller has to free the initial array.
 *
 * @returns A new sequence represention the union of the array.
 */
seq_t dsa_join( dsa_t *dsa){
	seq_t joined = { NULL, NULL, 0, 0, NULL, 0.0};
	if( dsa->size() == 0){
		return joined;
	}
	
	// Compute the total length
	size_t total = 0;
	for( auto it = dsa->begin(); it != dsa->end(); ++it){
		total += it->len + 1;
	}
	
	// A single malloc for the whole new sequence
	char *ptr = (char *) malloc( total);
	if( ptr == NULL ){
		return joined;
	}
	char *next = ptr;
	
	// Copy all old sequences and add a `!` in between
	auto it = dsa->begin();
	memcpy( next, it->S, it->len);
	next += it->len;
	for( it++; it != dsa->end(); it++){
		*next++ = '!';
		memcpy( next, it->S, it->len);
		next += it->len;
	}
	
	// Dont forget the null byte.
	*next = '\0';
	
	joined.S = ptr;
	joined.len = total -1; // subtract the null byte
	
	return joined;
}

/**
 * @brief Frees the memory of a given sequence.
 * @param S - The sequence to free.
 */
void free_seq( seq_t *S){
	free( S->S);
	free( S->RS);
	free( S->name);
	S->S = S->RS = S->name = NULL;
	S->len = S->RSlen = 0;
	S->gc = 0.0;
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
	char local_non_acgt = 0;
	
	while( len --> 0 ){
		c = *s--;
		
		switch( c){
			case 'A': d = 'T'; break;
			case 'T': d = 'A'; break;
			case 'G': d = 'C'; break;
			case 'C': d = 'G'; break;
			case '!': d = ';'; break; // rosebud
			default:
				local_non_acgt = 1; 
				continue;
		}
		
		*r++ = d;
	}
	
	if( local_non_acgt ){
		#pragma omp atomic
		FLAGS |= F_NON_ACGT;
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

/**
 * @brief Calculates the GC content of a sequence.
 *
 * This function computes the relative amount of G and C in the total sequence.
 */
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

/**
 * @brief Initializes a sequence.
 *
 * This function prepares a sequence for further processing.
 */
void init_seq( seq_t *S){
	normalize( S);
	
	if( !S->len){
		S->len = strlen(S->S);
	}
	calc_gc(S);
	
	S->RS = catcomp(S->S, S->len);
	S->RSlen = 2 * S->len + 1;
}

/**
 * @brief Restricts a sequence characters set to ACGT.
 *
 * This function strips a sequence of non ACGT characters and converts acgt to
 * the upper case equivalent. A flag is set if a non-canonical character was encountered.
 */
void normalize( seq_t *S){
	char *p, *q;
	char local_non_acgt = 0;
	for( p= q= S->S; *p; p++){
		switch( *p){
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case '!':
				*q++ = *p;
				break;
			case 'a':
			case 'c':
			case 'g':
			case 't':
				*q++ = toupper( (unsigned char)*p);
				break;
			default:
				local_non_acgt = 1;
				break;
				
		}
	}
	*q = '\0';
	if ( local_non_acgt ){
		#pragma omp atomic
		FLAGS |= F_NON_ACGT;
	}
}

