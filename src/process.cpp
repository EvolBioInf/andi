
/**
 * @file
 * @brief This file contains various distance methods.
 */
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "esa.h"
#include "global.h"
#include "process.h"
#include "sequence.h"
#include "io.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RMQ_n_1_improved.hpp>




double shuprop( size_t x, double g, size_t l);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor is allowed to be random.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject.
 * @returns The minimum length of an anchor.
 */
size_t minAnchorLength( double p, double g, size_t l){
	size_t x = 1;
	
	double prop = 0.0;
	while( prop < 1 - p){
		prop = shuprop( x, g/2, l);
		x++;
	}
	
	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We used to use gsl_sf_lnchoose(xx,kk) for this functionality.
 * After all, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL dependency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
 * @returns (n choose k)
 */
size_t binomial_coefficient( size_t n, size_t k){
	if( n <= 0 || k > n){
		return 0;
	}
	
	if( k == 0 || k == n ){
		return 1;
	}
	
	if( k > n-k ){
		k = n-k;
	}
	
	size_t res = 1;

	for( size_t i= 1; i <= k; i++){
		res *= n - k + i;
		res /= i;
	}
	
	return res;
}

/**
 * @brief Given `x` this function calculates the probability of a shustring 
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position. Then
 * this function computes P{X <= x} with respect to the given parameter set. See
 * Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shuprop( size_t x, double p, size_t l){
	double xx = (double)x;
	double ll = (double)l;
	size_t k;
	
	double s = 0.0;
	
	for(k=0; k<= x; k++){
		double kk = (double)k;
		double t = pow(p,kk) * pow(0.5 - p, xx - kk);
		
		s += pow(2,xx) * (t * pow(1-t,ll)) * (double)binomial_coefficient(x,k);
		if( s >= 1.0){
			s = 1.0;
			break;
		}
	}

	return s;
}

/**
 * @brief Divergence estimation using the anchor technique.
 *
 * The dist_anchor() function estimates the divergence between two
 * DNA sequences. The subject is given as an ESA, whereas the query
 * is a simple string. This function then looks for *anchors* -- long
 * substrings that exist in both sequences. Then it manually checks for
 * mutations between those anchors.
 * 
 * @return An estimate for the number of mutations within homologous regions.
 * @param C - The enhanced suffix array of the subject.
 * @param query - The actual query string.
 * @param query_length - The length of the query string. Needed for speed reasons.
 */
data_t dist_anchor( const esa_t *C, const char *query, size_t query_length, double gc){
	size_t snps = 0; // Total number of found SNPs
	size_t homo = 0; // Total number of homologous nucleotides.
	
	lcp_inter_t inter;
	
	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a pair.
	size_t last_was_right_anchor = 0;
	
	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;
	
	size_t num_right_anchors = 0;

#ifdef DEBUG
	size_t num_matches = 0;
	size_t num_anchors = 0;
	size_t num_anchors_in_rc = 0;
	size_t length_anchors = 0;
#endif

	size_t threshold = minAnchorLength( 1-sqrt(1-RANDOM_ANCHOR_PROP), gc, C->len);

	data_t retval = {0.0,0.0};

	// Iterate over the complete query.
	while( this_pos_Q < query_length){
		inter = get_match_cached( C, query + this_pos_Q, query_length - this_pos_Q);

#ifdef DEBUG
		num_matches++;
#endif

		this_length = inter.l <= 0 ? 0 : inter.l;
		
		if( inter.i == inter.j && this_length >= threshold)
		{
			// We have reached a new anchor.
			this_pos_S = C->SA[ inter.i];

#ifdef DEBUG
			num_anchors++;
			length_anchors += this_length;
			if( this_pos_S < (size_t)(C->len / 2)){
				num_anchors_in_rc++;
			}
#endif

			// Check if this can be a right anchor to the last one.
			if( this_pos_Q - last_pos_Q == this_pos_S - last_pos_S ){
				num_right_anchors++;
			
				// Count the SNPs in between.
				size_t i;
				for( i= last_length; i< this_pos_Q - last_pos_Q; i++){
					if( C->S[ last_pos_S + i] != query[ last_pos_Q + i] ){
						snps++;
					}
				}
				homo += this_pos_Q - last_pos_Q;
				last_was_right_anchor = 1;
			} else {
				if( last_was_right_anchor){
					// If the last was a right anchor, but with the current one, we 
					// cannot extend, then add its length.
					homo += last_length;
				}
				
				last_was_right_anchor = 0;
			}
			
			// Cache values for later
			last_pos_Q = this_pos_Q;
			last_pos_S = this_pos_S;
			last_length= this_length;
		}
		
		// Advance
		this_pos_Q += this_length + 1;
	}

#ifdef DEBUG
	if( FLAGS & F_EXTRA_VERBOSE ){
		const char str[] = {
			"- threshold: %ld\n"
			"- matches: %lu\n"
			"- anchors: %lu\n"
			"- in reverse complement: %lu\n"
			"- right anchors: %lu\n"
			"- avg length: %lf\n"
			"\n"
		};

		#pragma omp critical
		{
			fprintf(stderr, str, threshold, num_matches, num_anchors, num_anchors_in_rc, num_right_anchors, (double)length_anchors/ num_anchors );
		}
	}
#endif
	
	// Very special case: The sequences are identical
	if( last_length >= query_length ){
		retval.coverage = 1.0;
		return retval;
	}
	
	// We might miss a few nucleotides if the last anchor was also a right anchor.
	if( last_was_right_anchor ){
		homo += last_length;
	}
	
	// Nearly identical sequences
	if( homo == query_length){
		retval.distance = (double)snps/(double)homo;
		retval.coverage = 1.0;
		return retval;
	}

	// Insignificant results. All abort the fail train.
	if ( homo <= 3){
		retval.distance = log(-1.0);
		return retval;
	}
	
	// Abort if we have more homologous nucleotides than just nucleotides. This might
	// happen with sequences of different lengths.
	if( homo >= (size_t) C->len ){
		retval.distance = log(-1.0);
		retval.coverage = 1.0;
		return retval;
	}
	
	retval.distance = (double)snps/(double)homo;
	retval.coverage = (double)homo/(double)query_length;
	return retval;
}

/**
 * @brief Computes the distance matrix.
 *
 * The distMatrix() populates the D matrix with computed distances. It allocates D and
 * fills it with useful values, but the caller has to free it!
 * @return The distance matrix
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
#define FAST
#include "dist_hack.h"

/**
 * @brief Computes the distance matrix.
 *
 * The distMatrixLM() populates the D matrix with computed distances. It allocates D and
 * filles it with useful values, but the caller has to free it!
 * @return The distance matrix
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
#undef FAST
#include "dist_hack.h"

/**
 * @brief Calculates and prints the distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void calcDistMatrix( seq_t* sequences, int n){
	int i;

	// check the sequences
	for( i=0;i<n;i++){
		if( sequences[i].S == NULL){
			errx(1,"Missing sequence.");
		}
	}
	
	// Warn about non ACGT residues.
	if( FLAGS & F_NON_ACGT ){
		warnx("The input sequences contained characters other than acgtACGT. "
			"These were automatically stripped to ensure correct results.");
	}

	data_t *M = NULL;
	
	if( FLAGS & F_VERBOSE){
		errno = 0;
		M = (data_t*) malloc(n*n*sizeof(data_t));
		if( !M){
			warn("Couldn't allocate enough memory for verbose mode; Continuing without.");
			FLAGS &= ~F_VERBOSE;
		}
	}

	// compute the distances
	double *D = FLAGS & F_LOW_MEMORY ? distMatrixLM( sequences, n, M) : distMatrix( sequences, n, M);
	
	// print the results
	printDistMatrix( D, sequences, n);

	// print additional information.
	if( FLAGS & F_VERBOSE){
		printCovMatrix( M, n);
	}
	
	free(D);
	free(M);
}


