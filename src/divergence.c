/***** divergence.c ***********************************************************
 * Description: divergence calculation
 * The divergence computation is based on the mathematical model described in:
 * Haubold, B., Pfaffelhuber, P., Domazet-Loso, M., and Wiehe, T. 2009.
 * Estimating mutation distances from unaligned genomes, J. Comput. Biol.
 *
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Modified by: Mirjana Domazet-Loso, 02/12/2008
 * Modified by: Bernhard Haubold, August 24, 2013
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "divergence.h"
#include <gsl/gsl_sf_gamma.h>
#include "gsl/gsl_nan.h"

double expShulen(double d, double p, long l, double pS);
double pmax(long x, double p, long l, double pS);
double expShulenSimple( double d, long seqLen);

/** 
 * divergence
 * @returns the divergence (aka. mutation rate) which best explains the observed mean shulength
 * @param shulen  The mean observed length f shustrings.
 * @param seqLen  The length of the subject sequence.
 * @param gc  Relative gc-content of the query.
 * @param gcS  Relative gc-content of the subject.
 */
double divergence(double shulen, long seqLen, double gc, double gcS) {
	double p, q;
	double du, dl, dm, t, d, pS, es;
	//double errd = args -> E; /* relative error for shulen length */

	p = gc / 2.0;
	q = (1.0 - gc) / 2.0;
	pS = gcS / 2.0;
	du = 0;
	dl = 1.0 - (2 * p * p + 2 * q * q); // dl < 0.75
	t = 1.0e-3;//THRESHOLD;

	while ((dl - du) / 2.0 > t) {
		dm = (du + dl) / 2.0;
		
		//if (args -> f){
		//	es = expShulenSimple(dm, seqLen);
		//} else {
			es = expShulen(dm, p, seqLen, pS);
		//}
		if (shulen < es) {
			du = dm;
		} else {
			dl = dm;
		}
		/* test the relative error between du and dl; if it is smaller than some threshold, then break !!*/
		if (fabs(dl - du) / dl <= 1.0e-3) {
			break;
		}
	}
	d = (du + dl) / 2.0;
	return d;
}

double expShulenSimple( double d, long seqLen) {
	int x;
	float p, sl, cur, prev;
	float pow4;

	sl = 0;
	x = 1;
	p = 0;
	prev = 0;
	pow4 = 4;
	while (1) {
		cur = (1. - exp(-x * d)) * pow(1. - 1. / pow4, seqLen);
		p = (cur - prev) * x;
		if (p <= FLT_MIN && sl > 0){
			break;
		}
		sl += p;
		prev = cur;
		x++;
		pow4 *= 4;
	}
	return sl;
}

double expShulen(double d, double p, long l, double pS) {
	long i;

	double e = 0.0; /* expectation */
	double prob_i, delta;
	double factor;
	double probOld = 0;
	double last = 0.0;

	double t = 1.0 - d; // d = d' / l
	double p_t = t; /* pow(t, 1)*/

	//absolute error: double t1 = 1e-5; /* 1e-4 ok */
	double t1 = 1e-5; 
	//double t1 = args - > T;

	for (i = 1; i < l; i++) { //since for i = 0, the whole expression is 0
		factor = 1.0 - p_t; //factor = 1.0 - pow(1.0 - d, i);
		
		if( last < 1.0 ){
			last = pmax(i, p, l, pS);
			prob_i = factor * last;
		} else {
			prob_i = factor; /* prob_i = factor * s, where s = 1 */
		}
		delta = (prob_i - probOld) * i; /* delta should always be positive */
		e += delta; /* expectation of avg shulen l(Q, S)*/
		/* relative error - a little bit faster than the calculation with the absolute error */
		if (e >= 1 && delta / e <= t1) {
			break;
		}
		p_t *= t;
		probOld = prob_i;
	}
	return e;
}

double pmax(long x, double p, long l, double pS) {

	long k;
	double s = 0, x_choose_k, xx = (double)x, kk;
	double t, t1, m, delta;
	/* m_t value should be explored by simulation */
	const double m_t = pow(10, DBL_MIN_10_EXP); //args -> M; /* 10^(-307) */

	s = 0;
	for (k = 0; k <= x; k++) {
		kk = (double) k;
		x_choose_k = gsl_sf_lnchoose(x, k);
		m = pow(2.0, xx) * pow(p, kk) * pow(0.5 - p, xx-kk) * pow(1.0 - pow(pS, kk) * pow(0.5 - pS, xx-kk), (double) l);
		if (m == 0) {
			delta = 0;
		} else if (m >= m_t) {
			t = log(m);
			if (t == GSL_NEGINF) {
				delta = 0;
			} else {
				delta = exp(t + x_choose_k);
			}
		} else {
			t1 = log(1 + m); // for small values of m - to avoid overflow (-INF)
			delta = exp(t1 + x_choose_k) - exp(x_choose_k);
		}
		s += delta;
		if (s >= 1.0) {
			s = 1.0;
			break;
		}
	} /* end for */
	
	return s;
}


