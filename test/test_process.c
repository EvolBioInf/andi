#include "global.h"
#include "process.h"
#include <glib.h>
#include <math.h>

int FLAGS = 0;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
double ANCHOR_P_VALUE = 0.025;
gsl_rng *RNG = NULL;
int MODEL = M_JC;

double shuprop(size_t x, double g, size_t l);
size_t min_anchor_length(double p, double g, size_t l);

void test_shuprop() {
	int len = 100000;
	double gc = 0.5;
	double p_value = 0.025;

	size_t threshold = min_anchor_length(p_value, gc, len);

	g_assert_cmpfloat(1 - p_value, <, shuprop(threshold + 1, gc / 2, len));
	g_assert_cmpfloat(1 - p_value, <=, shuprop(threshold, gc / 2, len));
	g_assert_cmpfloat(1 - p_value, >, shuprop(threshold - 1, gc / 2, len));
}

int main(int argc, char *argv[]) {
	g_test_init(&argc, &argv, NULL);
	g_test_add_func("/process/shuprop", test_shuprop);

	return g_test_run();
}
