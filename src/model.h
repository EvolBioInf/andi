/** @file
 * @brief This header contains all structures and prototypes for creating a
 * mutation matrix and estimating distances trough an evolutionary model
 * thereof.
 */
#pragma once

#include <stdlib.h>

/**
 * This enum contains all possible mutations. The total number
 * of different possible mutations is MUTCOUNTS.
 */
enum {
	AtoA,
	AtoC,
	AtoG,
	AtoT,
	CtoA,
	CtoC,
	CtoG,
	CtoT,
	GtoA,
	GtoC,
	GtoG,
	GtoT,
	TtoA,
	TtoC,
	TtoG,
	TtoT,
	MUTCOUNTS
};

/** @brief The mutation matrix.
 *
 * We need to keep track of the different types of mutations between two
 * sequences. For this the following matrix is filled.
 *
 *  To   A  C  G  T
 * From
 *  A  (            )
 *  C  (            )
 *  G  (            )
 *  T  (            )
 *
 * The cells are absolute counts. Together with seq_len (the query length),
 * we can deduce the substitution rate and coverage.
 *
 * As libdivsufsort is 32 bit the sequence length is limited to (INT_MAX-1)/2.
 * We can thus use the same limit for the counts.
 */
typedef struct model {
	/** The absolute counts of mutation types. */
	unsigned int counts[MUTCOUNTS];
	/** The query length. */
	unsigned int seq_len;
} model;

void model_count_equal(model *, const char *, size_t);
void model_count(model *, const char *, const char *, size_t);
model model_average(const model *, const model *);
double model_coverage(const model *);
double estimate_RAW(const model *);
double estimate_JC(const model *);
double estimate_KIMURA(const model *);
double estimate_LOGDET(const model *);
model model_bootstrap(model);
