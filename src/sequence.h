/**
 * @file
 * @brief Functions and structures for DNA sequences
 *
 */
#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <vector>

/**
 * @brief A structure for sequences.
 *
 * This structure is used to represent a DNA sequence of some kind.
 */
typedef struct {
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

void free_seq( seq_t *S);
void init_seq( seq_t *S);

/**
 * A dynamicly growing structure for sequences.
 */
typedef std::vector<seq_t> dsa_t;

/** Create a new dynamic array for sequences. */
inline dsa_t *init_dsa(){
	return new std::vector<seq_t>();
}

/** Add a sequence to an array. */
inline void push_dsa( dsa_t *dsa, seq_t S){
	dsa->push_back( S);
}

/** Frees the array and all sequences stored within. */
inline void free_dsa( dsa_t *dsa){
	size_t i = 0;
	for(i=0; i < dsa->size(); i++){
		free_seq( &(dsa->at(i)));
	}
	delete dsa;
}

/** Returns the number of sequences stored within an array. */
inline size_t size_dsa( dsa_t *dsa){
	return dsa->size();
}

/** Get the raw C array. */
inline seq_t* data_dsa( dsa_t *dsa){
	return dsa->data();
}

seq_t join_dsa( dsa_t *dsa);

#endif

