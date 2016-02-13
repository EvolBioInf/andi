/**
 * @file
 * @brief This header contains the declarations for functions in esa.c.
 *
 */
#ifndef _ESA_H_
#define _ESA_H_

#include <sys/types.h>
#include "sequence.h"
#include "config.h"

#ifdef HAVE_LIBDIVSUFSORT
#include <divsufsort.h>
#else

#include "../opt/psufsort/interface.h"

typedef int saidx_t;

#endif

/**
 * @brief Represents LCP-Intervals.
 *
 * This struct is used to represent LCP-intervals. The member `i` should
 * coincide with the lower bound whereas `j` is the upper bound. Both bounds
 * are inclusive. So if `i == j` the interval contains exactly one element,
 * namely `i`. To represent an empty interval please use `i == j == -1`.
 * Other variants, such as `i == j == -2` can be used to indicate an error.
 * The common prefix length is denoted by l and should always be non-negative.
 * Variables of this type are often called `ij`.
 */
typedef struct {
	/** @brief The common prefix length */
	saidx_t l;
	/** @brief lower bound */
	saidx_t i;
	/** @brief upper bound */
	saidx_t j;
	/** The new middle. */
	saidx_t m;
} lcp_inter_t;

/**
 * @brief The ESA type.
 *
 * This structure holds arrays and objects associated with an enhanced
 * suffix array (ESA).
 */
typedef struct esa_s {
	/** The base string from which the ESA was generated. */
	const char *S;
	/** The actual suffix array with indexes into S. */
	saidx_t *SA;
	/** The LCP holds the number of letters up to which a suffix `S[SA[i]]`
		equals `S[SA[i-1]]`. Hence the name longest common prefix. For `i = 0`
		and `i = len` the LCP value is -1. */
	saidx_t *LCP;
	/** The length of the string S. */
	saidx_t len;
	/** A cache for lcp-intervals */
	lcp_inter_t *cache;
	/** The FVC array holds the character after the LCP. */
	char *FVC;
	/** This is the child array. */
	saidx_t *CLD;
} esa_s;

lcp_inter_t get_match_cached(const esa_s *, const char *query, size_t qlen);
lcp_inter_t get_match(const esa_s *, const char *query, size_t qlen);
int esa_init(esa_s *, const seq_t *S);
void esa_free(esa_s *);

#ifdef DEBUG

char code2char(ssize_t code);

#endif // DEBUG

#endif
