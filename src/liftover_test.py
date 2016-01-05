#! python

import os, sys, math, unittest
# from [library] import [function]

class test001(unittest.TestCase):
	"""test something"""
	def test_run(self):
		self.assertTrue(True)

if __name__=="__main__":
	unittest.main(verbosity=2)
