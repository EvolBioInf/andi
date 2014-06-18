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

dsa_t *init_dsa();
void   push_dsa( dsa_t *dsa, seq_t S);
void   free_dsa( dsa_t *dsa);
size_t size_dsa( dsa_t *dsa);
seq_t* data_dsa( dsa_t *dsa);
seq_t join_dsa( dsa_t *dsa);

#endif

