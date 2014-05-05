#include "RMQ_succinct.hpp"
#include "RMQ_1_n.hpp"
#include <iostream>
using namespace std;

/** Create a random array of length 5,000,000 and test the succinct data structure
 ** for RMQ by posing 1,000 random RMQs to the array and compare the result
 ** with the naive algorithm for RMQ.
 **/
int
main(int argc, char* argv[]) {
	DTidx n = 5000000;
 	srand((unsigned) time(NULL));
	cout << "constructing array of size "<< n << " and preprocessing for RMQ... \n";
	DT* a = new DT[n];
	for (DTidx i = 0; i < n; i++) a[i] = rand();

 	RMQ_1_n  RMQ(a, n);       // this is the naive algorithm
   	RMQ_succinct RMQ2(a, n);  // this is our succinct constant-time query algorithm
	cout << "... preprocessing done, now testing correctness (may take a while!)...\n";

	// test some combinations:
	DTidx i, j;
	for (DTidx x = 0; x < 1000; x++) {
		i = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		j = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		if (i > j) { DTidx tmp=i; i=j;j=tmp; }
		if (a[RMQ.query(i,j)] != a[RMQ2.query(i,j)]) {
			cout << "ERROR: " << i << "," << j << endl;
			cout << "ERROR: " << RMQ.query(i,j) << "," << RMQ2.query(i,j) << endl;
			cout << x << endl;
			exit(0);
		}
	}

	// demonstrate performance
	cout << "... testing done, now demonstrating performance. First, naive method..\n";
	for (DTidx x = 0; x < 1000; x++) {
		i = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		j = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		if (i > j) { DTidx tmp=i; i=j;j=tmp; }
		RMQ.query(i,j);
	}

	cout << "... now our method...\n";
	for (DTidx x = 0; x < 1000; x++) {
		i = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		j = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		if (i > j) { DTidx tmp=i; i=j;j=tmp; }
		RMQ2.query(i,j);
	}
	cout << "... done.\n";
}
