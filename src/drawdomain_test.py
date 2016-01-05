#! python

import os, sys, math, unittest, colorsys
import drawdomain
import numpy as np

class test001(unittest.TestCase):
	"""color mapping"""
	def test_run(self):
		rgbstr='aabbcc'
		rgb = drawdomain.hexToRgb(rgbstr)
		self.assertTrue(rgb == (170, 187, 204))
		rgbtup = (12,50,100)
		rgbstr = drawdomain.rgbToHex(rgbtup)
		self.assertTrue(rgbstr == '0c3264')

class test002(unittest.TestCase):
	"""composition"""
	def test_run(self):
		seq = 'ACDEFACDEF'
		comp = drawdomain.sequenceCompositionNorm(seq)
		self.assertAlmostEqual((comp*comp).sum(), 1.0)

class test003(unittest.TestCase):
	"""colors"""
	def test_run(self):
		rgbtup = (0.1, 0.5, 0.8)
		hsv = colorsys.rgb_to_hsv(*rgbtup)
		print hsv
		rgb = colorsys.hsv_to_rgb(*hsv)
		print rgb
		vrgb = np.r_[rgb]
		print vrgb
		rgb8 = drawdomain.v01To8bit(vrgb)
		print rgb8
		hexc = drawdomain.rgbToHex(rgb8)
		print hexc
		rgb = drawdomain.hexToRgb(hexc)
		print rgb
		ycbcr = drawdomain.ycbcrToRgb(*drawdomain.v01To8bit(vrgb))
		print ycbcr
		rgb = drawdomain.rgbToYcbcr(*ycbcr)
		print rgb


if __name__=="__main__":
	unittest.main(verbosity=2)
