#!/bin/sh -f

# Test if andi exists, and can be executed
./src/andi --version || exit 1

# Test andi for more than just two sequences at a time
./test/test_fasta -l 100000 -d 0.01 -d 0.01 -d 0.01 -d 0.01 | ./src/andi > /dev/null


