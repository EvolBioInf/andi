#/usr/bin/bash -f
# Simulate the anchor distance using vmatch.

# These are all sequences from the ECO29 set.
SEQS="AE005174.fasta AE005674.fasta AE014073.fasta AE014075.fasta AP009048.fasta AP009240.fasta BA000007.fasta CP000034.fasta CP000036.fasta CP000038.fasta CP000243.fasta CP000247.fasta CP000266.fasta CP000468.fasta CP000800.fasta CP000802.fasta CP000946.fasta CP000948.fasta CP000970.fasta CP001063.fasta CP001396.fasta CP001846.fasta CU928160.fasta CU928161.fasta CU928162.fasta CU928163.fasta CU928164.fasta FM180568.fasta U00096.fasta"

# Loop over all sequences
for S in $SEQS; do
	echo $S;
	Q=${SEQS/$S/}; # All sequences, except S

	# Create the index for S
	./mkvtree -db "$S" -dna -allout -pl

	# Recall, that an anchor is unique in S and of some minimum length. Hence
	# the parameters -mum cand and -l 16. The latter threshold was calculated
	# by andi.

	# Match all other sequences including their reverse against S
	./vmatch -q $Q -mum cand -l 16 -d -p -nodist -noscore -noevalue -noidentity $S > /dev/null

	# In theory we have to calculate the anchors distance here using the previously
	# computed matches. But since vmatch is already significantly slower than andi,
	# we skip this step.
done

