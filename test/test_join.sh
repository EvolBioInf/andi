#!/bin/sh -f

./src/andi --help > /dev/null || exit 1

SEED=${RANDOM_SEED:-0}
SEED2=0
SEED3=0
if test $SEED -ne 0; then
        SEED=$((SEED + 1))
        SEED2=$((SEED + 2))
        SEED3=$((SEED + 3))
fi

# Simple join test
./test/test_fasta -s $SEED -l 1000 -L 1000 -d 0.1 > p1_join.fasta
./test/test_fasta -s $SEED2 -l 1000 -L 1000 -d 0.1 > p2_join.fasta
./test/test_fasta -s $SEED3 -l 10000 -L 10000 -d 0.1 > p3_join.fasta

head -qn 2 p1_join.fasta p2_join.fasta p3_join.fasta > S0_join.fasta
tail -qn 2 p1_join.fasta p2_join.fasta p3_join.fasta > S1_join.fasta

rm p1_join.fasta p2_join.fasta p3_join.fasta;


RES=$(./src/andi -m RAW -t 1 -j S0_join.fasta S1_join.fasta |
	tail -n 1 |
	awk '{print ($2 - 0.1)}' |
	awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) < 0.03}'
	)

if test $RES -ne 1; then
	echo "The last test computed a distance deviating more than three percent from its intended value."
	echo "See S0_join.fasta and S1_join.fasta for the used sequences."
	exit 1;
fi

SEED=${RANDOM_SEED:-0}
SEED2=0
if test $SEED -ne 0; then
        SEED=$((SEED + 5))
        SEED2=$((SEED + 6))
fi

#unbalanced number of contigs
./test/test_fasta -s $SEED -l 1000 -L 1000 -d 0.1 > p2_join.fasta
./test/test_fasta -s $SEED2 -l 10000 -L 10000 -d 0.1 > p3_join.fasta

head -qn 2 p3_join.fasta > S0_join.fasta
tail -qn 2 p2_join.fasta p3_join.fasta > S1_join.fasta

rm p2_join.fasta p3_join.fasta;


RES=$(./src/andi -m RAW -t1 -j S0_join.fasta S1_join.fasta |
        tail -n 1 |
        awk '{print ($2 - 0.1)}' |
        awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) < 0.03}'
        )

if test $RES -ne 1; then
        echo "The last test computed a distance deviating more than three percent from its intended value."
        echo "See S0_join.fasta and S1_join.fasta for the used sequences."
        exit 1;
fi

SEED=${RANDOM_SEED:-0}
SEED2=0
SEED3=0
if test $SEED -ne 0; then
        SEED=$((SEED + 11))
        SEED2=$((SEED + 12))
        SEED3=$((SEED + 13))
fi

#unbalanced number of contigs 2
./test/test_fasta -s $SEED -l 1000 -L 1000 -d 0.1 > p1_join.fasta
./test/test_fasta -s $SEED2 -l 1000 -L 1000 -d 0.1 > p2_join.fasta
./test/test_fasta -s $SEED3 -l 10000 -L 10000 -d 0.1 > p3_join.fasta

head -qn 2 p1_join.fasta p3_join.fasta > S0_join.fasta
tail -qn 2 p1_join.fasta p2_join.fasta p3_join.fasta > S1_join.fasta

rm p1_join.fasta p2_join.fasta p3_join.fasta;


RES=$(./src/andi -mRAW -t 1 -j S0_join.fasta S1_join.fasta |
        tail -n 1 |
        awk '{print ($2 - 0.1)}' |
        awk 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($1-$2) < 0.03}'
        )

if test $RES -ne 1; then
        echo "The last test computed a distance deviating more than three percent from its intended value."
        echo "See S0_join.fasta and S1_join.fasta for the used sequences."
        exit 1;
fi


rm S0_join.fasta S1_join.fasta
