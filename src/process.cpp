/**
 * @file
 * @brief This file contains various distance methods.
 */
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "esa.h"
#include "global.h"
#include "process.h"
#include "sequence.h"

#include <RMQ_n_1_improved.hpp>

/**
 * This is a neat hack for dealing with matrices.
 */
#define D( X, Y) (D[ (X)*n + (Y) ])

double shuprop( size_t x, double g, size_t l);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The propability with which an anchor is allowed to be random.
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
 * Afterall, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL depenency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
 * @returns (n choose k)
 */
size_t binomial_coefficient( size_t n, size_t k){
	if( n <= 0 || k < 0 || k > n){
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
 * @brief Given `x` this function calculates the propability of a shustring 
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position. Then
 * this function computes P{X <= x} with respect to the given parameter set. See
 * Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The propability of a certain shustring length.
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
double dist_anchor( const esa_t *C, const char *query, size_t query_length, double gc){
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
	
	size_t threshhold = minAnchorLength( 1-sqrt(RANDOM_ANCHOR_PROP), gc, C->len);
	if( FLAGS & F_VERBOSE){
		fprintf(stderr, "threshhold: %ld\n", threshhold);
	}
	
	// Iterate over the complete query.
	while( this_pos_Q < query_length){
		inter = getLCPInterval( C, query + this_pos_Q, query_length - this_pos_Q);
		this_length = inter.l;
		if( this_length == 0) break;
		
		if( inter.i == inter.j && this_length >= threshhold)
		{
			// We have reached a new anchor.
			this_pos_S = C->SA[ inter.i];
			
			// Check if this can be a right anchor to the last one.
			if( this_pos_Q - last_pos_Q == this_pos_S - last_pos_S ){
				num_right_anchors++;
			
				// Count the SNPs in between.
				size_t i;
				for( i= 0; i< this_pos_Q - last_pos_Q; i++){
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
	
	// Very special case: The sequences are identical
	if( last_length >= query_length ){
		return 0.0;
	}
	
	// We might miss a few nucleotides if the last anchor was also a right anchor.
	if( last_was_right_anchor ){
		homo += last_length;
	}
	
	// Nearly identical sequences
	if( homo == query_length){
		return (double)snps/(double)homo;
	}
	
	if ( num_right_anchors <= 1 || snps <= 2 || homo <= 3){
		// Insignificant results. All abort the fail train.
		return log(-1.0);
	}
	
	// Abort if we have more homologous nucleotides than just nucleotides. This might
	// happen with sequences of different lengths.
	if( homo >= C->len ){
		return log(-1.0);
	}

	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "snps: %lu, homo: %lu\n", snps, homo);
		fprintf( stderr, "number of right anchors: %lu\n", num_right_anchors);
	}
	
	return (double)snps/(double)homo;
}

/**
 * @brief Computes the distance matrix.
 *
 * The distMatrix() populates the D matrix with computed distances. It allocates D and
 * filles it with useful values, but the caller has to free it!
 * @return The distance matrix
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
double *distMatrix( seq_t* sequences, int n){
	double *D = (double*)malloc( n * n * sizeof(double));
	assert(D);
	
	double d;
	
	int i;

	#pragma omp parallel for num_threads( THREADS)
	for(i=0;i<n;i++){
		esa_t E = {NULL,NULL,NULL,NULL,0,NULL};
		
		// initialize the enhanced suffix array
		E.S = (const char*) sequences[i].RS;
		E.len = sequences[i].RSlen;
		
		int result;

		result = compute_SA( &E);
		assert( result == 0); // zero errors
		result = compute_LCP_PHI( &E);
		assert( result == 0);
	
		//E.rmq_lcp = new RMQ_succinct(E.LCP, E.len);
		E.rmq_lcp = new RMQ_n_1_improved(E.LCP, E.len);
		
		// now compare every other sequence to i
		int j;
		for(j=0; j<n; j++){
			if( j == i) {
				D(i,j) = 0.0;
				continue;
			}
			
			// TODO: Provide a nicer progress indicator.
			if( FLAGS & F_VERBOSE ){
				#pragma omp critical
				{
					fprintf( stderr, "comparing %d and %d\n", i, j);
				}
			}

			size_t ql = sequences[j].len;
			
			d = dist_anchor( &E, sequences[j].S, ql, sequences[i].gc);
			
			if( !(FLAGS & F_RAW)){
				d = -0.75 * log(1.0- (4.0 / 3.0) * d ); // jukes cantor
			}
			// fix negative zero
			if( d <= 0.0 ){
				d = 0.0;
			}
			D(i,j) = d;
		}
		
		delete E.rmq_lcp;
		free( E.SA);
		free( E.ISA);
		free( E.LCP);
	}
	
	return D;
}

/**
 * @brief Prints the distance matrix.
 *
 * This function pretty prints the distance matrix. For small distances
 * scientific notation is used.
 * @param D - The distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void printDistMatrix( double *D, seq_t* sequences, size_t n){

	int use_scientific = 0;
	size_t i,j;
	
	for( i=0; i<n && !use_scientific; i++){
		for( j=0; j<n; j++){
			if( D(i,j) > 0 && D(i,j) < 0.001 ){
				use_scientific = 1;
				break;
			}
		}
	}
	
	printf("%lu\n", n);
	for( i=0;i<n;i++){
		// Print exactly nine characters of the name. Padd wih spaces if necessary.
		printf("%-9.9s", sequences[i].name);
		
		for( j=0;j<n;j++){
			if( use_scientific){
				printf(" %1.4e", (D(i,j)+D(j,i))/2 );
			} else {
				printf(" %1.4f", (D(i,j)+D(j,i))/2 );
			}
		}
		printf("\n");
	}
}

/**
 * @brief Calculates and prints the distance matrix
 * @param sequences - An array of pointers to the sequences.
 * @param n - The number of sequences.
 */
void calcDistMatrix( seq_t* sequences, int n){
	int i;

	// initialise the sequences
	#pragma omp parallel for num_threads( THREADS)
	for( i=0;i<n;i++){
		if( sequences[i].S == NULL){
			fprintf(stderr, "missing sequence %d\n", i);
			exit(1);
		}
		
		init_seq( &sequences[i]);
	}
	
	// Warn about non ACGT residues.
	if( FLAGS & F_NON_ACGT ){
		const char str[] = {
			"The input sequences contained characters other than acgtACGT. "
			"These were automatically stripped to ensure correct results.\n"
		};
		fprintf( stderr, "%s", str);
	}
	
	// compute the distances
	double *D = distMatrix( sequences, n);
	
	// print the results
	printDistMatrix( D, sequences, n);
	
	free(D);
}


