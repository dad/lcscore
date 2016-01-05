#! python

import os, sys, math, unittest
import slide, gelscore
# from [library] import [function]

class test001(unittest.TestCase):
	def test_run(self):
		"""composition comparison"""
		a = gelscore.Sequence('GYPMGGYPMGGYPMGGYPMGGYPMG')
		b = gelscore.Sequence('GYPMG')
		fac = slide.SequenceCompositionComparatorFactory()
		sca = fac.make(a,'sequence')
		scb = fac.make(b,'sequence')
		self.assertAlmostEqual(sca.compare(sca),1.0) # Equal to self
		self.assertAlmostEqual(sca.compare(scb),1.0)
		self.assertAlmostEqual(scb.compare(sca),1.0)

class test002(unittest.TestCase):
	def test_run(self):
		"""sliding"""
		a = gelscore.Sequence('GYPMGGYPMGGYPMGGYPMGGYPMG')
		b = gelscore.Sequence('GYPMG')
		fac = slide.SequenceCompositionComparatorFactory()
		scb = fac.make(b,'sequence')
		sw = slide.SequenceWindow(a, 5)
		res = sw.slide(scb)
		self.assertAlmostEqual(sum(res.scores),sw.maxNumberOfWindows())

class test003(unittest.TestCase):
	def test_run(self):
		"""sliding, no matching"""
		a = gelscore.Sequence('GYPMGGYPMGGYPMGGYPMGGYPMG')
		b = gelscore.Sequence('AAAAA')
		fac = slide.SequenceCompositionComparatorFactory()
		scb = fac.make(b,'sequence')
		sw = slide.SequenceWindow(a, 5)
		res = sw.slide(scb)
		self.assertAlmostEqual(sum(res.scores),0.0)

class test004(unittest.TestCase):
	def test_run(self):
		"""sliding, no matching"""
		a = gelscore.Sequence('QLAQQIQARNQMRYQQATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGPVYGVPPQGGFPRNANDNNQFYQ')
		b = gelscore.Sequence('MPQNGRA')
		fac = slide.SequenceCompositionComparatorFactory()
		scb = fac.make(b,'sequence')
		sw = slide.SequenceWindow(a, len(b))
		ress = sw.slide(scb)
		for (xi, res) in enumerate(ress.results()):
			self.assertTrue(a[res.position-1] == a[xi])
			self.assertTrue(res.score >= 0.0)
		#self.assertAlmostEqual(sum(res.scores),0.0)

class test005(unittest.TestCase):
	def test_run(self):
		"""searching"""
		a = gelscore.Sequence('ACDEFGHIKLAAAAAAAAAAAMNPQRSTVWY')
		sw = slide.SequenceWindow(a, 5)
		res = sw.search()
		self.assertTrue(str(res.currentSequence()) == 'AAAAA')

class test006(unittest.TestCase):
	def test_run(self):
		"""iterating"""
		a = gelscore.Sequence('ACDEFGHIKLAAAAAAAAAAAMNPQRSTVWY')
		#b = gelscore.Sequence('AAAAA')
		winsize = 5
		sw = slide.SequenceWindow(a, winsize)
		xi = 0
		while sw.isValid():
			seq = sw.currentSequence()
			#print seq, str(a[xi:(xi+winsize)])
			self.assertTrue(str(seq) == str(a[xi:(xi+winsize)]))
			sw.next()
			xi += 1

class test007(unittest.TestCase):
	def test_run(self):
		"""straight comparison"""
		a = gelscore.Sequence('AAAAA')
		b = gelscore.Sequence('AAAAA')
		c = gelscore.Sequence('LAAAA')
		fac = slide.SequenceCompositionComparatorFactory()
		ac = fac.make(a, 'sequence')
		bc = fac.make(b, 'sequence')
		cc = fac.make(c, 'sequence')
		self.assertTrue(ac.compare(bc)>ac.compare(cc))

class test008(unittest.TestCase):
	def test_run(self):
		"""window larger than sequence"""
		a = gelscore.Sequence('MLSLIFYLRFPSYIRG') # Actual protein sequence YJR151W-A
		sw = slide.SequenceWindow(a, 20)


if __name__=="__main__":
	unittest.main(verbosity=2)
