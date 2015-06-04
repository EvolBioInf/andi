/**
 * @file
 * @brief This header contains function declarations for io procedures.
 */
#ifndef _IO_H_
#define _IO_H_

#include <err.h>
#include <errno.h>
#include "sequence.h"
#include <stdio.h>

/** @brief The data structure can be used to store output data resulting 
 * from the computation of distance.
 */
typedef struct data_s {
	/** The distance */
	double distance;
	/** The coverage */
	double coverage;
} data_t;


/**
 * This is a neat hack for dealing with matrices.
 */
#define D( X, Y) (D[ (X)*n + (Y) ])
#define M( X, Y) (M[ (X)*n + (Y) ])

void read_fasta( FILE *in, dsa_t *dsa);
void read_fasta_join( FILE *in, dsa_t *dsa, const char *name);

void print_distances( const double *D, const seq_t *sequences, size_t n);
void print_coverages( const data_t *D, size_t n);

#endif // _IO_H_
