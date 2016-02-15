/**
 * @file
 * @brief ESA functions
 *
 * This file contains various functions that operate on an enhanced suffix
 * array. The basic algorithms originate from the book of Ohlebusch
 * "Bioinformatics Algorithms" (2013). Most of these were heavily modified
 * for improved performance. One example is the lcp-cache.
 *
 * The ESA structure defined in esa.h contains a `cache` field. This cache is
 * used to quickly look up lcp-intervals. Consider the queries "AAGT" and
 * "AACG". In both cases the interval for "AA" has to be looked up in the
 * ESA. If we simply store the interval for "AA" in the cache, once and use it
 * for each query we are significantly faster (up to 7 times).
 */
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "esa.h"
#include "global.h"

static void esa_init_cache_dfs(esa_s *, char *str, size_t pos, lcp_inter_t in);
static void esa_init_cache_fill(esa_s *, char *str, size_t pos, lcp_inter_t in);

static lcp_inter_t get_interval(const esa_s *, lcp_inter_t ij, char a);
lcp_inter_t get_match(const esa_s *, const char *query, size_t qlen);
static lcp_inter_t get_match_from(const esa_s *, const char *query, size_t qlen,
								  saidx_t k, lcp_inter_t ij);

static int esa_init_SA(esa_s *);
static int esa_init_LCP(esa_s *);
static int esa_init_CLD(esa_s *);

/** @brief The prefix length up to which LCP-intervals are cached. */
const size_t CACHE_LENGTH = 10;

/** @brief Map a code to the character. */
char code2char(ssize_t code) {
	switch (code & 0x3) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
	}
	return '\0';
}

/** @brief Map a character to a two bit code. */
ssize_t char2code(const char c) {
	ssize_t result = -1;
	switch (c) {
		case 'A': result = 0; break;
		case 'C': result = 1; break;
		case 'G': result = 2; break;
		case 'T': result = 3; break;
	}
	return result;
}

#define R(CLD, i) ((CLD)[(i)])
#define L(CLD, i) ((CLD)[(i)-1])

/** @brief Fills the LCP-Interval cache.
 *
 * Traversing the virtual suffix tree, created by SA, LCP and CLD is rather
 * slow. Hence we create a cache, holding the LCP-interval for a prefix of a
 * certain length ::CACHE_LENGTH. This function it the entry point for the
 * cache filling routine.
 *
 * @param self - The ESA.
 * @returns 0 iff successful
 */
int esa_init_cache(esa_s *self) {
	lcp_inter_t *cache = malloc((1 << (2 * CACHE_LENGTH)) * sizeof(*cache));
	CHECK_MALLOC(cache);

	self->cache = cache;

	char str[CACHE_LENGTH + 1];
	str[CACHE_LENGTH] = '\0';

	saidx_t m = L(self->CLD, self->len);
	lcp_inter_t ij = {.i = 0, .j = self->len - 1, .m = m, .l = self->LCP[m]};

	esa_init_cache_dfs(self, str, 0, ij);

	return 0;
}

/** @brief Fills the cache — one char at a time.
 *
 * This function is a depth first search on the virtual suffix tree and fills
 * the cache. Or rather it calls it self until some value to cache is
 * calculated. This function is a recursive version of get_inteval but with more
 * edge cases.
 *
 * @param C - The ESA.
 * @param str - The current prefix.
 * @param pos - The length of the prefix.
 * @param in - The LCP-interval of prefix[0..pos-1].
 */
void esa_init_cache_dfs(esa_s *C, char *str, size_t pos, const lcp_inter_t in) {
	// we are not yet done, but the current strings do not exist in the subject.
	if (pos < CACHE_LENGTH && in.i == -1 && in.j == -1) {
		esa_init_cache_fill(C, str, pos, in);
		return;
	}

	// we are past the caching length
	if (pos >= CACHE_LENGTH) {
		esa_init_cache_fill(C, str, pos, in);
		return;
	}

	lcp_inter_t ij;

	// iterate over all nucleotides
	for (int code = 0; code < 4; ++code) {
		str[pos] = code2char(code);
		ij = get_interval(C, in, str[pos]);

		// fail early
		if (ij.i == -1 && ij.j == -1) {
			esa_init_cache_fill(C, str, pos + 1, ij);
			continue;
		}

		if (ij.l <= (ssize_t)(pos + 1)) {
			// Continue one level deeper
			// This is the usual case
			esa_init_cache_dfs(C, str, pos + 1, ij);
			continue;
		}

		// The LCP-interval is deeper than expected
		// Check if it still fits into the cache
		if ((size_t)ij.l >= CACHE_LENGTH) {
			// If the lcp-interval exceeds the cache depth, stop here and fill
			esa_init_cache_fill(C, str, pos + 1, in);
			continue;
		}

		/* At this point the prefix `str` of length `pos` has been found.
		 * However, the call to `getInterval` above found an interval with
		 * an LCP value bigger than `pos`. This means that not all elongations
		 * (more precise: just one) of `str` appear in the subject. Thus fill
		 * all values with the matched result to far and continue only with
		 * the one special substring.
		 */
		esa_init_cache_fill(C, str, pos + 1, in);

		char non_acgt = 0;

		// fast forward
		size_t k = pos + 1;
		for (; k < (size_t)ij.l; k++) {
			// In some very edgy edge cases the lcp-interval `ij`
			// contains a `;` or another non-acgt character. Since we
			// cannot cache those, break.
			char c = C->S[C->SA[ij.i] + k];
			if (char2code(c) < 0) {
				non_acgt = 1;
				break;
			}

			str[k] = c;
		}

		if (non_acgt) {
			esa_init_cache_fill(C, str, k, ij);
		} else {
			esa_init_cache_dfs(C, str, k, ij);
		}
	}
}

/** @brief Fills the cache with a given value.
 *
 * Given a prefix and a value this function fills the cache beyond this point
 * the value.
 *
 * @param C - The ESA.
 * @param str - The current prefix.
 * @param pos - The length of the prefix.
 * @param in - The LCP-interval of prefix[0..pos-1].
 */
void esa_init_cache_fill(esa_s *C, char *str, size_t pos, lcp_inter_t in) {
	if (pos < CACHE_LENGTH) {
		for (int code = 0; code < 4; ++code) {
			str[pos] = code2char(code);
			esa_init_cache_fill(C, str, pos + 1, in);
		}
	} else {
		ssize_t code = 0;
		for (size_t i = 0; i < CACHE_LENGTH; ++i) {
			code <<= 2;
			code |= char2code(str[i]);
		}

		C->cache[code] = in;
	}
}

/**
 * @brief Initializes the FVC (first variant character) array.
 *
 * The FVC is of my own invention and simply defined as
 * `FVC[i] = S[SA[i]+LCP[i]]`. This expression is constantly used in
 * get_interval. By precomputing the result, we have less memory
 * accesses, less cache misses, and thus improved runtimes of up to 15%
 * faster matching. This comes at a negligible cost of increased memory.
 *
 * @param self - The ESA
 * @returns 0 iff successful
 */
int esa_init_FVC(esa_s *self) {
	size_t len = self->len;

	char *FVC = self->FVC = malloc(len);
	CHECK_MALLOC(FVC);

	const char *S = self->S;
	const int *SA = self->SA;
	const int *LCP = self->LCP;

	FVC[0] = '\0';
	for (size_t i = len; i; i--, FVC++, SA++, LCP++) {
		*FVC = S[*SA + *LCP];
	}

	return 0;
}

/** @brief Initializes an ESA.
 *
 * This function initializes an ESA with respect to the provided sequence.
 * @param C - The ESA to initialize.
 * @param S - The sequence
 * @returns 0 iff successful
 */
int esa_init(esa_s *C, const seq_t *S) {
	if (!C || !S || !S->S) return 1;

	*C = (esa_s){.S = S->RS, .len = S->RSlen};

	int result;

	result = esa_init_SA(C);
	if (result) return result;

	result = esa_init_LCP(C);
	if (result) return result;

	result = esa_init_CLD(C);
	if (result) return result;

	result = esa_init_FVC(C);
	if (result) return result;

	result = esa_init_cache(C);
	if (result) return result;

	return 0;
}

/** @brief Free the private data of an ESA. */
void esa_free(esa_s *self) {
	free(self->SA);
	free(self->LCP);
	free(self->CLD);
	free(self->cache);
	free(self->FVC);
	*self = (esa_s){};
}

/**
 * Computes the SA given a string S. To do so it uses libdivsufsort.
 * @param C The enhanced suffix array to use. Reads C->S, fills C->SA.
 * @returns 0 iff successful
 */
int esa_init_SA(esa_s *C) {
	// assert c.S
	if (!C || !C->S) {
		return 1;
	}

	C->SA = malloc(C->len * sizeof(*C->SA));
	CHECK_MALLOC(C->SA);

	saidx_t result = 1;

#ifdef HAVE_LIBDIVSUFSORT
	result = divsufsort((const unsigned char *)C->S, C->SA, C->len);
#else
	result = c_psufsort(C->S, C->SA);
#endif

	return result;
}

/** @brief Initializes the CLD (child) array.
 *
 * See Ohlebusch.
 *
 * @param C - The ESA
 */
int esa_init_CLD(esa_s *C) {
	if (!C || !C->LCP) {
		return 1;
	}
	saidx_t *CLD = C->CLD = malloc((C->len + 1) * sizeof(*CLD));
	CHECK_MALLOC(CLD);

	saidx_t *LCP = C->LCP;

	typedef struct pair_s { saidx_t idx, lcp; } pair_t;

	pair_t *stack = malloc((C->len + 1) * sizeof(*stack));
	CHECK_MALLOC(stack);
	pair_t *top = stack; // points at the topmost filled element
	pair_t last;

	R(CLD, 0) = C->len + 1;

	top->idx = 0;
	top->lcp = -1;

	// iterate over all elements
	for (size_t k = 1; k < (size_t)(C->len + 1); k++) {
		while (LCP[k] < top->lcp) {
			// top->lcp is a leaf
			last = *top--;

			// link all elements of same lcp value in a chain
			while (top->lcp == last.lcp) {
				R(CLD, top->idx) = last.idx;
				last = *top--;
			}

			// store the l-index of last
			if (LCP[k] < top->lcp) {
				R(CLD, top->idx) = last.idx;
			} else {
				L(CLD, k) = last.idx;
			}
		}

		// continue one level deeper
		top++;
		top->idx = k;
		top->lcp = LCP[k];
	}

	free(stack);
	return 0;
}

/**
 * This function computed the LCP array, given the suffix array. Thereto it uses
 * a special `phi` array, which makes it slightly faster than the original
 * linear-time algorithm by Kasai et al.
 *
 * @param C The enhanced suffix array to compute the LCP from.
 * @returns 0 iff successful
 */
int esa_init_LCP(esa_s *C) {
	const char *S = C->S;
	saidx_t *SA = C->SA;
	saidx_t len = C->len;

	// Trivial safety checks
	if (!C || !S || !SA || len == 0) {
		return 1;
	}

	// Allocate new memory
	// The LCP array is one element longer than S.
	saidx_t *LCP = C->LCP = malloc((len + 1) * sizeof(*LCP));
	CHECK_MALLOC(LCP);

	LCP[0] = -1;
	LCP[len] = -1;

	// Allocate temporary arrays
	saidx_t *PHI = malloc(len * sizeof(*PHI));
	saidx_t *PLCP = PHI;
	CHECK_MALLOC(PHI);

	PHI[SA[0]] = -1;
	saidx_t k;
	ssize_t i;

	for (i = 1; i < len; i++) {
		PHI[SA[i]] = SA[i - 1];
	}

	ssize_t l = 0;
	for (i = 0; i < len; i++) {
		k = PHI[i];
		if (k != -1) {
			while (S[k + l] == S[i + l]) {
				l++;
			}
			PLCP[i] = l;
			l--;
			if (l < 0) l = 0;
		} else {
			PLCP[i] = -1;
		}
	}

	// unpermutate the LCP array
	for (i = 1; i < len; i++) {
		LCP[i] = PLCP[SA[i]];
	}

	free(PHI);
	return 0;
}

/** @brief For the lcp-interval of string `w` compute the interval for `wa`
 *
 * Say, we already know the LCP-interval ij for a string `w`. Now we want to
 * check if `wa` may also be found in the ESA and thus in the subject. So we
 * look for the sub interval of `ij` in which all strings feature an `a` as
 * the next character. If such a sub interval is found, its boundaries are
 * returned.
 *
 * @param self - The ESA.
 * @param ij - The lcp-interval for `w`.
 * @param a - The next character.
 * @returns The lcp-interval one level deeper.
 */
static lcp_inter_t get_interval(const esa_s *self, lcp_inter_t ij, char a) {
	saidx_t i = ij.i;
	saidx_t j = ij.j;

	const saidx_t *SA = self->SA;
	const saidx_t *LCP = self->LCP;
	const char *S = self->S;
	const saidx_t *CLD = self->CLD;
	const char *FVC = self->FVC;
	// check for singleton or empty interval
	if (i == j) {
		if (S[SA[i] + ij.l] != a) {
			ij.i = ij.j = -1;
		}
		return ij;
	}

	int m = ij.m;
	int l = ij.l;

	char c = S[SA[i] + l];
	goto SoSueMe;

	do {
		c = FVC[i];

	SoSueMe:
		if (c == a) {
			/* found ! */
			saidx_t n = L(CLD, m);

			ij = (lcp_inter_t){.i = i, .j = m - 1, .m = n, .l = LCP[n]};

			return ij;
		}

		if (c > a) {
			break;
		}

		i = m;

		if (i == j) {
			break; // singleton interval, or `a` not found
		}

		m = R(CLD, m);
	} while (/*m != "bottom" && */ LCP[m] == l);

	// final sanity check
	if (i != ij.i ? FVC[i] == a : S[SA[i] + l] == a) {
		ij.i = i;
		ij.j = j;
		/* Also return the length of the LCP interval including `a` and
		 * possibly even more characters. Note: l + 1 <= LCP[m] */
		ij.l = LCP[m];
		ij.m = m;
	} else {
		ij.i = ij.j = -1;
	}

	return ij;
}

/** @brief Compute the longest match of a query with the subject.
 *
 * The *longest match* is the core concept of `andi`. Its simply defined as the
 * longest prefix of a query Q appearing anywhere in the subject S. Talking
 * about genetic sequences, a match is a homologous region, likely followed by a
 * SNP.
 *
 * This function returns the interval for where the longest match of the query
 * can be found in the ESA. Thereto it expects a starting interval for the
 * search.
 *
 * @param C - The enhanced suffix array for the subject.
 * @param query - The query sequence.
 * @param qlen - The length of the query. Should correspond to `strlen(query)`.
 * @param k - The starting index into the query.
 * @param ij - The LCP interval for the string `query[0..k]`.
 * @returns The LCP interval for the longest prefix.
 */
lcp_inter_t get_match_from(const esa_s *C, const char *query, size_t qlen,
						   saidx_t k, lcp_inter_t ij) {

	if (ij.i == -1 && ij.j == -1) {
		return ij;
	}

	// fail early on singleton intervals.
	if (ij.i == ij.j) {

		// try to extend the match. See line 513 below.
		saidx_t p = C->SA[ij.i];
		size_t k = ij.l;
		const char *S = (const char *)C->S;

		for (; k < qlen && S[p + k]; k++) {
			if (S[p + k] != query[k]) {
				ij.l = k;
				return ij;
			}
		}

		ij.l = k;
		return ij;
	}

	saidx_t l, i, j;

	lcp_inter_t res = ij;

	const saidx_t *SA = C->SA;
	const char *S = C->S;

	// Loop over the query until a mismatch is found
	do {
		// Get the subinterval for the next character.
		ij = get_interval(C, ij, query[k]);
		i = ij.i;
		j = ij.j;

		// If our match cannot be extended further, return.
		if (i == -1 && j == -1) {
			res.l = k;
			return res;
		}

		res.i = ij.i;
		res.j = ij.j;

		l = qlen;
		if (i < j && ij.l < l) {
			/* Instead of making another look up we can use the LCP interval
			 * calculated in get_interval */
			l = ij.l;
		}

		// By definition, the kth letter of the query was matched.
		k++;

		// Extend the match
		for (int p = SA[i]; k < l; k++) {
			if (S[p + k] != query[k]) {
				res.l = k;
				return res;
			}
		}
	} while (k < (ssize_t)qlen);

	res.l = qlen;
	return res;
}

/** @brief Get a match.
 *
 * Given an ESA and a string Q find the longest prefix of Q that matches
 * somewhere in C. This search is done entirely via jumping around in the ESA,
 * and thus is slow.
 *
 * @param C - The ESA.
 * @param query - The query string — duh.
 * @param qlen - The length of the query.
 * @returns the lcp interval of the match.
 */
lcp_inter_t get_match(const esa_s *C, const char *query, size_t qlen) {
	// sanity checks
	if (!C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->CLD) {
		return (lcp_inter_t){-1, -1, -1, -1};
	}

	saidx_t m = L(C->CLD, C->len);
	lcp_inter_t ij = {.i = 0, .j = C->len - 1, .m = m, .l = C->LCP[m]};

	return get_match_from(C, query, qlen, 0, ij);
}

/** @brief Compute the LCP interval of a query. For a certain prefix length of
 * the query its LCP interval is retrieved from a cache. Hence this is faster
 * than the naive `get_match`. If the cache fails to provide a proper value, we
 * fall back to the standard search.
 *
 * @param C - The enhanced suffix array for the subject.
 * @param query - The query sequence.
 * @param qlen - The length of the query. Should correspond to `strlen(query)`.
 * @returns The LCP interval for the longest prefix.
 */
lcp_inter_t get_match_cached(const esa_s *C, const char *query, size_t qlen) {
	if (qlen <= CACHE_LENGTH) return get_match(C, query, qlen);

	ssize_t offset = 0;
	for (size_t i = 0; i < CACHE_LENGTH && offset >= 0; i++) {
		offset <<= 2;
		offset |= char2code(query[i]);
	}

	if (offset < 0) {
		return get_match(C, query, qlen);
	}

	lcp_inter_t ij = C->cache[offset];

	if (ij.i == -1 && ij.j == -1) {
		return get_match(C, query, qlen);
	}

	return get_match_from(C, query, qlen, ij.l, ij);
}
