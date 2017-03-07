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
	/** The length of the forward strand. */
	size_t len;
	/** A name for this sequence */
	char *name;
} seq_t;

/**
 * @brief This structure enhances the usual sequence with its reverse
 * complement.
 */
typedef struct seq_subject {
	/** This member contains first the reverse strand and then the
		forward strand. */
	char *RS;
	/** Corresponds to strlen(RS) */
	size_t RSlen;
	/**
	 * @brief GC-Content
	 *
	 * The relative amount of G or C in the DNA.
	 */
	double gc;
} seq_subject;

void seq_free(seq_t *S);
int seq_subject_init(seq_subject *S, const seq_t *);
void seq_subject_free(seq_subject *S);
int seq_init(seq_t *S, const char *seq, const char *name);

/**
 * @brief A dynamically growing structure for sequences.
 */
typedef struct dsa_s {
	seq_t *data;
	size_t capacity, size;
} dsa_t;

int dsa_init(dsa_t *A);
void dsa_push(dsa_t *A, seq_t S);
void dsa_free(dsa_t *A);
size_t dsa_size(dsa_t *A);
seq_t *dsa_data(dsa_t *A);

seq_t dsa_join(dsa_t *dsa);

#endif
