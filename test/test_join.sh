#!/bin/sh -f

./test/test_fasta -l 1000 -L 1000 -d 0.1 > p1.fasta
./test/test_fasta -l 1000 -L 1000 -d 0.1 > p2.fasta

head -qn 2 p1.fasta p2.fasta > S0.fasta
tail -qn 2 p1.fasta p2.fasta > S1.fasta

rm p1.fasta p2.fasta;

RES=$($srcdir/src/andi -r -j S0.fasta S1.fasta |
	tail -n 1 |
	awk '{print ($2 - 0.1)}' |
	awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) < 0.01}'
	)

if test $RES -ne 1; then
	echo "The last test computed a distance deviating more than one percent from its intended value."
	echo "See S1.fasta and S2.fasta for the used sequences."
	exit 1;
fi

rm S0.fasta S1.fasta
