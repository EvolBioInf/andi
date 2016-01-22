/** @file
 * @brief This file contains all functions for the mutation matrix and the
 * estimation of evolutionary distances thereof.
 */

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>
#include "global.h"
#include "model.h"

/**
 * @brief Sum some mutation count specified by `summands`. Intended to be used
 * through the `model_sum` macro.
 *
 * @param MM - The mutation matrix.
 * @param summands - The mutations to add.
 * @returns The sum of mutations.
 */
static size_t model_sum_types(const model *MM, const int summands[]) {
	size_t total = 0;
	for (int i = 0; summands[i] != MUTCOUNTS; ++i) {
		total += MM->counts[summands[i]];
	}
	return total;
}

#define model_sum(MM, ...)                                                     \
	model_sum_types((MM), (int[]){__VA_ARGS__, MUTCOUNTS})

/**
 * @brief Average two mutation matrices.
 *
 * @param MM - One matrix
 * @param NN - Second matrix
 * @returns The average (sum) of two mutation matrices.
 */
model model_average(const model *MM, const model *NN) {
	model ret = *MM;
	for (int i = 0; i != MUTCOUNTS; ++i) {
		ret.counts[i] += NN->counts[i];
	}
	ret.seq_len += NN->seq_len;
	return ret;
}

/**
 * @brief Compute the total number of nucleotides in the pairwise alignment.
 *
 * @param MM - The mutation matrix.
 * @returns The length of the alignment.
 */
size_t model_total(const model *MM) {
	size_t total = 0;
	for (size_t i = 0; i < MUTCOUNTS; ++i) {
		total += MM->counts[i];
	}
	return total;
}

/**
 * @brief Compute the coverage of an alignment.
 *
 * @param MM - The mutation matrix.
 * @returns The relative coverage
 */
double model_coverage(const model *MM) {
	size_t covered = model_total(MM);
	size_t actual = MM->seq_len;

	return (double)covered / (double)actual;
}

/**
 * @brief Estimate the uncorrected distance of a pairwise alignment.
 *
 * @param MM - The mutation matrix.
 * @returns The uncorrected substitution rate.
 */
double estimate_RAW(const model *MM) {
	size_t nucl = model_total(MM);
	size_t SNPs = model_sum(MM, AtoC, AtoG, AtoT, CtoG, CtoT, GtoT);

	// Insignificant results. All abort the fail train.
	if (nucl <= 3) {
		return NAN;
	}

	return (double)SNPs / (double)nucl;
}

/**
 * @brief Compute the Jukes-Cantor distance.
 *
 * @param MM - The mutation matrix.
 * @returns The corrected JC distance.
 */
double estimate_JC(const model *MM) {
	double dist = estimate_RAW(MM);
	dist = -0.75 * log(1.0 - (4.0 / 3.0) * dist); // jukes cantor

	// fix negative zero
	return dist <= 0.0 ? 0.0 : dist;
}

/* @brief Bootstrap a mutation matrix.
 *
 * The classical bootstrapping process, as described by Felsenstein, resamples
 * all nucleotides of a MSA. As andi only computes a pairwise alignment, this
 * process boils down to a simple multinomial distribution. We just have to
 * resample the elements of the mutation matrix. (Paper in review.)
 *
 * @param MM - The original mutation matrix.
 * @returns A bootstrapped mutation matrix.
 */
model model_bootstrap(const model MM) {
	model datum = MM;

	size_t nucl = model_total(&MM);
	double p[MUTCOUNTS];
	for (size_t i = 0; i < MUTCOUNTS; ++i) {
		p[i] = MM.counts[i] / (double)nucl;
	}

	unsigned int n[MUTCOUNTS];

	gsl_ran_multinomial(RNG, MUTCOUNTS, nucl, p, n);

	for (size_t i = 0; i < MUTCOUNTS; ++i) {
		datum.counts[i] = n[i];
	}

	return datum;
}

/**
 * @brief Count the substitutions and add them to the mutation matrix.
 *
 * @param MM - The mutation matrix.
 * @param S - The subject
 * @param Q - The query
 * @param len - The length of the alignment
 */
void model_count(model *MM, const char *S, const char *Q, size_t len) {
	for (; len--; S++, Q++) {
		char s = *S;
		char q = *Q;

		// Skip special characters.
		if (s == ';' || s == '!' || s == '#' || q == ';' || q == '!') {
			continue;
		}

		size_t index = 0;

		if (s == q) {
			switch (s) {
				case 'A': index = AtoA; break;
				case 'C': index = CtoC; break;
				case 'G': index = GtoG; break;
				case 'T': index = TtoT; break;
			}

			MM->counts[index]++;
			continue;
		}

		if (s > q) {
			char c = s;
			s = q;
			q = c;
		}

		switch (s) {
			case 'A':
				switch (q) {
					case 'C': index = AtoC; break;
					case 'G': index = AtoG; break;
					case 'T': index = AtoT; break;
				}
				break;
			case 'C':
				switch (q) {
					case 'G': index = CtoG; break;
					case 'T': index = CtoT; break;
				}
				break;
			case 'G': index = GtoT; break;
		}

		MM->counts[index]++;
	}
}
