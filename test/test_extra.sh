#!/bin/sh -f

# Test if andi exists, and can be executed
./src/andi --version > /dev/null || exit 1

# Test andi for more than just two sequences at a time
./test/test_fasta -l 100000 -d 0.01 -d 0.01 -d 0.01 -d 0.01 | ./src/andi > /dev/null || exit 1

# Test low-memory mode
./test/test_fasta -l 10000 > test_extra.fasta
./src/andi test_extra.fasta > extra.out
./src/andi test_extra.fasta --low-memory > extra_low_memory.out
diff extra.out extra_low_memory.out || exit 1

# Test file of filenames
./test/test_fasta -l 10000 > test_extra.fasta
echo "$PWD/test_extra.fasta" > fof.txt
./src/andi test_extra.fasta > extra.out
./src/andi --file-of-filenames fof.txt > fof.out
cat fof.txt | ./src/andi --file-of-filenames - > fof2.out
diff extra.out fof.out || exit 1
diff extra.out fof2.out || exit 1


rm -f test_extra.fasta extra.out extra_low_memory.out fof.out fof2.out fof.txt

