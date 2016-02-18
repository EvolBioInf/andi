#!/usr/bin/python

import sys
import re
import math

def map2d(fn, mat):
	return list(map(lambda row: list(map(fn,row)),mat))

DIVERGENCE = 0.02

def approx(f1, f2):
	"""Check absolute and relative divergence."""
	return abs(f1-f2) <= DIVERGENCE and abs(f1-f2) <= DIVERGENCE * f2

class Matrix(object):
	"""Distance Matrix"""
	def __init__(self, size, lines):
		rLine = re.compile(r'^(\S+)\s+((\S+\s*)*)$')
		self.size = len(lines)
		parts = [rLine.search(line).groups() for line in lines]
		self.names = [part[0] for part in parts]
		sdist = [part[1].split() for part in parts]
		self.dist = map2d(float, sdist)

	def __repr__(self):
		s = str(self.size) + "\n"
		formatted = map2d((lambda f: "%.4e" % f), self.dist)
		merge = [" ".join(row) for row in formatted]
		lines = ["%-9.9s %s" % (name, line) for (name,line) in zip(self.names, merge)]
		return s + "\n".join(lines)

	def proper(self):
		size = self.size
		dist = self.dist
		for i in range(size):
			assert len(dist[i]) == size
			for j in range(size):
				assert not math.isnan(dist[i][j])
				if i == j:
					assert dist[i][j] == 0
				if i < j:
					assert dist[i][j] == dist[j][i]

	def approx(self, other):
		assert self.size == other.size
		for i in range(size):
			for j in range(size):
				assert approx(self.size, other.size)

def m2dApprox(m, dist):
	"""Check if the matrix m is approximately [[0.0,dist],[dist,0.0]]"""
	m.proper()
	assert m.size == 2
	assert approx(m.dist[0][1], dist)
	assert approx(m.dist[1][0], dist)


def parse(s):
	lines = s.splitlines()
	ret = []
	i = 0

	while i < len(lines) and lines[i].isnumeric():
		size = int(lines[i])
		m = Matrix(size, lines[i+1:i+1+size])
		ret.append(m)
		i += size + 1

	assert i == len(lines)

	return ret



def main(args):
	for arg in args:
		if arg == "-":
			file = sys.stdin
		else:
			file = open(arg)

		matrices = parse(file.read())
		for m in matrices:
			m.proper()
			print(m)

		file.close()

# main
if __name__ == '__main__':
	main(sys.argv[1:])

