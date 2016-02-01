/** @file
 * @brief This header contains all structures and prototypes for creating a
 * mutation matrix and estimating distances trough an evolutionary model
 * thereof.
 */
#pragma once

#include <stdlib.h>

/**
 * This enum contains all possible mutations. Some of them are considered
 * symmetric (e.g. AtoT and TtoA) and thus have the same value. The total number
 * of different possible mutations is MUTCOUNTS.
 */
enum {
	AtoA,
	AtoC,
	AtoG,
	AtoT,
	CtoC,
	CtoG,
	CtoT,
	GtoG,
	GtoT,
	TtoT,
	MUTCOUNTS,
	CtoA = AtoC,
	GtoA = AtoG,
	GtoC = CtoG,
	TtoA = AtoT,
	TtoC = CtoT,
	TtoG = GtoT
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
 */
typedef struct model {
	/** The absolute counts of mutation types. */
	size_t counts[MUTCOUNTS];
	/** The query length. */
	size_t seq_len;
} model;

void model_count_equal(model *, const char *, size_t);
void model_count(model *, const char *, const char *, size_t);
model model_average(const model *, const model *);
double model_coverage(const model *);
double estimate_RAW(const model *);
double estimate_JC(const model *);
double estimate_KIMURA(const model *);
model model_bootstrap(const model);
