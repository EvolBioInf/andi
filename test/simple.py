#!/usr/bin/python3

# This scripts test the accuracy of andi with random inputs. For that
# it uses the small program test_random to generate pairs of sequences
# with a given distance. By default, test_random creates a new set of
# sequences each time it is called. Thus, this test has a small, but
# non-zero probability of failing. That is a problem with Debian's
# reproducible builds effort. So this script acts as a wrapper around
# this issue.
#
# Simply calling this script via
#     % ./test/test_random.sh
# checks a new test-case every time. But with the right parameter
#     % RANDOM_SEED=1729 ./test/test_random.sh
# one specific set of sequences is validated.

from __future__ import print_function
import sys
import os
import matrix
import subprocess
import re

def warnx(*objs):
	print("warnx:", *objs, file=sys.stderr)

SEED = 0
LENGTH = 100000


def tSimple(dist):
	"""
	Generates two sequences with distance $dist and pipe them
	through andi. Then checks if the estimated distance is reasonably
	close to $dist.
	"""
	global SEED

	cmd_genFasta = ('./test/test_fasta', '-l', str(LENGTH), '-d', str(dist), '-s', str(SEED))
	cmd_tee = ('tee', 'temp.fasta')
	cmd_andi = ('./src/andi', '-t1')

	if SEED != 0:
		SEED += 1

	try:
		sub_genFasta = subprocess.Popen(cmd_genFasta, stdout=subprocess.PIPE)
	except Exception as e:
		warnx("Calling ./test/test_fasta failed. Maybe it was not build, yet?")
		raise e

	sub_tee = subprocess.Popen(cmd_tee, stdin=sub_genFasta.stdout, stdout=subprocess.PIPE)

	try:
		res_andi = subprocess.check_output(cmd_andi, stdin=sub_tee.stdout).decode("utf-8")
		sub_tee.wait()
	except Exception as e:
		# pipeline failed
		warnx("Whoops, pipeline failed to run. Does andi execute successfully?")
		raise e

	try:
		ms = matrix.parse(res_andi)
	except Exception as e:
		warnx("That was not a distance matrix:\n" + res_andi)
		raise e

	if len(ms) != 1:
		warnx("Why did andi produce more than one distance matrix?")

	try:
		matrix.m2dApprox(ms[0], dist)
	except Exception as e:
		# read actual seed from temp.fasta
		that_seed = SEED
		if SEED == 0:
			file = open("temp.fasta").read()
			rSeed = re.compile(r'base_seed: (\d+)')
			that_seed = int(rSeed.search(file).groups()[0])

		cmd_genFasta = cmd_genFasta[:-1] + (str(that_seed),)
		warnx("The last distance diverges more than two percent from its intended value.\n" +
		 "Reproducible settings:\n" +
		 "  " + " ".join(cmd_genFasta) + " |\n" +
		 "  " + " ".join(cmd_tee) + " |\n" +
		 "  " + " ".join(cmd_andi))
		exit(1)
		


if __name__ == '__main__':
	# derp
	SEED = int(os.getenv('BASE_SEED', 0))
	for dist in [0.0, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]:
		for _ in range(10):
			tSimple(dist)

	#rm temp.fasta
