/**
 * @file
 * @brief Functions and structures for DNA sequences
 *
 */
#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <err.h>
#include <errno.h>
#include <stdlib.h>

/**
 * @brief A structure for sequences.
 *
 * This structure is used to represent a DNA sequence of some kind.
 */
typedef struct seq_s {
	/** This is the DNAs forward strand as a string. */
	char *S;
	/** This member contains first the reverse strand and then the
		forward strand. */
	char *RS;
	/** The length of the forward strand. */
	size_t len;
	/** Corresponds to strlen(RS) */
	size_t RSlen;
	/** A name for this sequence */
	char *name;
	/**
	 * @brief GC-Content
	 *
	 * The relative amount of G or C in the DNA.
	 */
	double gc;
} seq_t;

void seq_free( seq_t *S);
void seq_subject_init( seq_t *S);
void seq_subject_free( seq_t *S);
int seq_init( seq_t *S, const char *seq, const char *name);

/**
 * A dynamically growing structure for sequences.
 */

typedef struct dsa_s {
	seq_t* data;
	size_t capacity, size;
} dsa_t;



/** Create a new dynamic array for sequences. */
inline int dsa_init(dsa_t *A){
	A->data = (seq_t*) malloc(sizeof(seq_t) * 2);
	if(!A->data){
		return 1;
	}
	//!FIXME check for null!
	A->capacity = 2;
	A->size = 0;
	return 0;
}

/** Add a sequence to an array. */
inline void dsa_push( dsa_t *A, seq_t S){
	if( A->size < A->capacity){
		A->data[A->size++] = S;
	} else {
		seq_t* ptr = (seq_t*) realloc(A->data, sizeof(seq_t) * A->capacity * 2);
		if(ptr == NULL){
			errx(errno, "out of memory?");
		}

		A->capacity *= 2;
		A->data = ptr;
		A->data[A->size++] = S;
	}
}

/** Frees the array and all sequences stored within. */
inline void dsa_free( dsa_t *A){
	size_t i;
	for( i=0; i< A->size; i++){
		seq_free(&A->data[i]);
	}

	free(A->data);
	A->data = NULL;
	A->capacity = A->size = 0;
}

/** Returns the number of sequences stored within an array. */
inline size_t dsa_size( dsa_t *A){
	return A->size;
}

/** Get the raw C array. */
inline seq_t* dsa_data( dsa_t *A){
	return A->data;
}

seq_t dsa_join( dsa_t *dsa);

#endif
