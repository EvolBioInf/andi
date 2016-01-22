/**
 * @file
 * @brief This header contains function declarations for io procedures.
 */
#ifndef _IO_H_
#define _IO_H_

#include <err.h>
#include <errno.h>
#include <stdio.h>
#include "sequence.h"
#include "model.h"

/**
 * This is a neat hack for dealing with matrices.
 */
#define D(X, Y) (D[(X)*n + (Y)])
#define M(X, Y) (M[(X)*n + (Y)])

void read_fasta(const char *, dsa_t *dsa);
void read_fasta_join(const char *, dsa_t *dsa);

void print_distances(const struct model *, const seq_t *, size_t, int);
void print_coverages(const struct model *, size_t);

#endif // _IO_H_
