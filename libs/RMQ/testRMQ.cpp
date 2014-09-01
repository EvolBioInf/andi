#include "RMQ_n_1_improved.hpp"
#include <iostream>
using namespace std;

main(int argc, char* argv[]) {
	srand((unsigned) time(NULL));
	DTidx n = 3000000;
	DT* a = new DT[n];
	
	for (DTidx i = 0; i < n; i++) a[i] = 1+(int) (4.0*rand()/(RAND_MAX+1.0));
	cout << "preprocessing array of size 3,000,000...";
	RMQ_n_1_improved RMQ(a,n);
	cout << "done!\n";

	// perform some queries:
	cout << "performing 1,000,000 random queries...";
	DTidx i, j;
	for (DTidx x = 0; x < 1000000; x++) {
		i = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		j = (DTidx) (n * ((float) rand() / (RAND_MAX + 1.0)));
  		if (i > j) { DTidx tmp=i; i=j;j=tmp; } //swap
		RMQ.query(i,j);
	}
	cout << "done, good bye!\n";
}
