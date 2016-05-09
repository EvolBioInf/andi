#!/usr/bin/zsh

# Compute the number of failing comparisons for different distances.

DISTS=(0.1 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7)

LENGTH=100000

for dist in $DISTS; do
	echo "" > est_$dist.dist
	for (( i = 0; i < 1000; i++ )); do
		../test/test_fasta -l $LENGTH -d $dist > temp.fa
		../src/andi ./temp.fa > est.dist 2> /dev/null
		tail -n 1 est.dist >> est_$dist.dist
	done
	avg=$(cat est_$dist.dist | awk '"nan" !~ $2 {sum+=$2;c++}END{print sum/c}')
	sd=$(grep -v 'nan' est_$dist.dist | awk '{a[c++]=$2;aa+=$2}END{aa/=NR;for(c=0;c<NR;c++){t=a[c]-aa;sd+=t*t}print sqrt(sd/(NR-1))}')
	failed=$(cat est_$dist.dist | grep -c 'nan')
	echo $dist "\t" $avg"\t±" $sd "\t" $failed
done
