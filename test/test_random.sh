#!/bin/sh -f

LENGTH=100000

for dist in 0.001 0.01 0.02 0.05 0.1 0.2 0.3
do
	for n in $(seq 10)
	do
		res=$($srcdir/test/test_fasta -l $LENGTH -d $dist |
			tee $srcdir/test/test_random.fasta |
			$srcdir/src/andi -r |
			tail -n 1 |
			awk -v dist=$dist '{print $2, dist}' |
			awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) < 0.01}')
		if test $res -ne 1; then
			echo "The last test computed a distance deviating more than one percent from its intended value."
			echo "See test_random.fasta for the used sequences."
			exit 1;
		fi
	done
done

rm $srcdir/test/test_random.fasta
