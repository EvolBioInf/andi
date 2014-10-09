/**
 * @file
 * @brief This header contains function declaratins for io procedures. It also 
 * defines the WARN and FAIL macro.
 */
#ifndef _IO_H_
#define _IO_H_

#include "sequence.h"

typedef struct data_s {
	double distance;
	double coverage;
} data_t;


/**
 * This is a neat hack for dealing with matrices.
 */
#define D( X, Y) (D[ (X)*n + (Y) ])

void readFile( FILE *in, dsa_t *dsa);
void joinedRead( FILE *in, dsa_t *dsa, char *name);

void printDistMatrix( data_t *D, seq_t *sequences, size_t n);
void printCovMatrix( data_t *D, size_t n);

/** @brief Exit early on a critical error.
 *
 * @returns never.
 */
#define FAIL(X) do { \
	const char str[] = { "ERROR: " X "\n"}; \
	fprintf(stderr, str); \
	exit(EXIT_FAILURE); \
} while (0)

/** @brief Print a warning.
 *
 */
#define WARN(X, ...) do { \
	_Pragma("omp critical") \
	{ \
		const char str[] = { "WARNING: " X "\n"}; \
		fprintf(stderr, str, __VA_ARGS__); \
	} \
} while (0)

#endif // _IO_H_
