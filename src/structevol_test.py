#! python

import os, sys, math, unittest
import structevol



class test001(unittest.TestCase):
	def test_run(self):
		"""make alignment"""
		L = 20
		for n in range(5,20):
			a = structevol.makeAlignment(n,L, 0.1, 0)
			#print a
			self.assertTrue(len(a)==n)
			self.assertTrue(len(a[0])==L)

class test002(unittest.TestCase):
	def test_run(self):
		"""hamming distance"""
		L = 20
		a = structevol.makeAlignment(10,10, 0.1, 0)
		ha = structevol.ProfileAlignment(a)
		for m in ha.distances():
			print m

if __name__=="__main__":
	unittest.main(verbosity=2)
