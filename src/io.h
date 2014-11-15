/**
 * @file
 * @brief This header contains function declarations for io procedures.
 */
#ifndef _IO_H_
#define _IO_H_

#include <err.h>
#include <errno.h>
#include "sequence.h"

typedef struct data_s {
	double distance;
	double coverage;
} data_t;


/**
 * This is a neat hack for dealing with matrices.
 */
#define D( X, Y) (D[ (X)*n + (Y) ])
#define M( X, Y) (M[ (X)*n + (Y) ])

void readFile( FILE *in, dsa_t *dsa);
void joinedRead( FILE *in, dsa_t *dsa, const char *name);

void printDistMatrix( double *D, seq_t *sequences, size_t n);
void printCovMatrix( data_t *D, size_t n);

#endif // _IO_H_
