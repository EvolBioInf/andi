#!/usr/bin/python3

from __future__ import print_function
import sys
import os
import matrix
import subprocess
import re

def warnx(*objs):
	print("warnx:", *objs, file=sys.stderr)

SEED = 0
LENGTH = 10000


def tSimple(dist):
	"""
	Generates two sequences with distance $dist and pipe them
	through andi. Then checks if the estimated distance is reasonably
	close to $dist.
	"""

	cmd_genFasta = ('./test/test_fasta', '-l', str(LENGTH), '-d', str(dist), '-s', str(SEED))
	cmd_tee = ('tee', 'temp.fasta')
	cmd_andi = ('./src/andi', '-t1')

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

	try:
		for m in ms:
			matrix.m2dApprox(m, dist)
	except Exception as e:
		# read actual seed from temp.fasta
		that_seed = SEED
		if SEED == 0:
			file = open("temp.fasta").read()
			rSeed = re.compile(r'base_seed: (\d+)')
			that_seed = rSeed.search(file).groups()[0]

		cmd_genFasta = cmd_genFasta[:-1] + (that_seed,)
		warnx("The last distance diverges more than two percent from its intended value.\n" +
		 "Reproducible settings:\n" +
		 "  " + " ".join(cmd_genFasta) + " |\n" +
		 "  " + " ".join(cmd_tee) + " |\n" +
		 "  " + " ".join(cmd_andi))
		exit(1)
		


if __name__ == '__main__':
	# derp
	SEED = os.getenv('BASE_SEED', 0)
	tSimple(0.1)
	tSimple(0.2)
	tSimple(0.3)
	tSimple(0.4)
	tSimple(0.5)
