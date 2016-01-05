#! python

import os, sys, math, unittest
import gelscore
import numpy as np
import scipy as sp
# from [library] import [function]

def printRegions(res):
	print ""
	for r in res:
		print r.region

class test_contig_seq_intersection1(unittest.TestCase):
	"""Merging two ContiguousRegions"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 10)
		gs2 = gelscore.ContiguousRegion(gs, 5, 15)
		#print gs1.composition
		#print gs2.composition
		gs3 = gs1.intersection(gs2)
		#print gs3
		self.assertTrue(str(gs3) == 'GHIKL')

class test_contig_seq_intersection2(unittest.TestCase):
	"""Merging two ContiguousRegions"""
	def test_run(self):
		gs = gelscore.Sequence('MSEAQETHVEQLPESVVDAPVEEQHQEPPQAPDAPQEPQVPQESAPQESAPQEPPAPQEQNDVPPPSNAPIYEGEESHSVQDYQEAHQHHQPPEPQPYYPPPPPGEHMHGRPPMHHRQEGELSNTRLFVRPFPLDVQESELNEIFGPFGPMKEVKILNGFAFVEFEEAESAAKAIEEVHGKSFANQPLEVVYSKLPAKRYRITMKNLPEGCSWQDLKDLARENSLETTFSSVNTRDFDGTGALEFPSEEILVEALERLNNIEFRGSVITVERDDNPPPIRRSNRGGFRGRGGFRGGFRGGFRGGFSRGGFGGPRGGFGGPRGGYGGYSRGGYGGYSRGGYGGSRGGYDSPRGGYDSPRGGYSRGGYGGPRNDYGPPRGSYGGSRGGYDGPRGDYGPPRDAYRTRDAPRERSPTR')
		gs1 = gelscore.ContiguousRegion(gs, 74, 106)
		gs2 = gelscore.ContiguousRegion(gs, 88, 111)
		#print gs1.composition
		#print gs2.composition
		gs3 = gs1.intersection(gs2)
		#print gs3.start
		self.assertTrue(gs3.start==88)
		self.assertTrue(gs3.end==106)

class test_contig_seq_diff(unittest.TestCase):
	"""Difference between two ContiguousRegions"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 15)
		gs2 = gelscore.ContiguousRegion(gs, 5, 20)
		[left, right] = gs1.difference(gs2)
		self.assertTrue(str(left) == 'ACDEF')
		self.assertTrue(str(right) == 'STVWY')

class test_contig_seq_diff_no_overlap(unittest.TestCase):
	"""Difference between two ContiguousRegions that don't overlap"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 5)
		gs2 = gelscore.ContiguousRegion(gs, 15, 20)
		[left, right] = gs1.difference(gs2)
		self.assertTrue(str(left) == str(gs1))
		self.assertTrue(str(right) == str(gs2))

class test_contig_seq_diff_subset(unittest.TestCase):
	"""Difference between two ContiguousRegions where one is a subset of the other"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 20)
		gs2 = gelscore.ContiguousRegion(gs, 10, 15)
		gsl = gelscore.ContiguousRegion(gs, 0, 10)
		gsr = gelscore.ContiguousRegion(gs, 15, 20)
		[left, right] = gs1.difference(gs2)
		self.assertTrue(str(left) == str(gsl))
		self.assertTrue(str(right) == str(gsr))

class test_contig_seq_inter_no_overlap(unittest.TestCase):
	"""Intersection between two ContiguousRegions that don't overlap"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 5)
		gs2 = gelscore.ContiguousRegion(gs, 15, 20)
		inter = gs1.intersection(gs2)
		self.assertTrue(inter is None)

class test_contig_seq_merge(unittest.TestCase):
	"""Merging two ContiguousRegions"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs, 0, 5)
		gs2 = gelscore.ContiguousRegion(gs, 5, 10)
		gs3 = gs1.merge(gs2)

class test_contig_seq_trimright(unittest.TestCase):
	"""ContiguousRegion trim right"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		gs1 = gelscore.ContiguousRegion(gs)
		gs1.trimright(2)
		self.assertTrue(gs1[-1] == 'V')

class test_find_simple_regions_ncut(unittest.TestCase):
	"""find ncut regions in obvious cases"""
	def test_run(self):
		gs = gelscore.Sequence('AAAAAAAAAAAAGGGGGGGGGGG')
		dist = gelscore.SimilarityWeight(sim_add=1.0)
		max_regions = 4
		comp = gelscore.SequenceCompositionSimilarity()
		rf = gelscore.NormalizedCutRegionFinder(dist, min_region_size=6, score_threshold=0.9, max_regions=max_regions)
		res = rf.find(gs)
		self.assertTrue(len(res) < max_regions)
		self.assertTrue(len(res) == 2)

class test_find_fus_subregions_ncut(unittest.TestCase):
	"""find ncut regions in FUS"""
	def test_run(self):
		gs = gelscore.Sequence('RGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGN')
		dist = gelscore.SimilarityWeight(sim_add=1.0)
		rf = gelscore.NormalizedCutRegionFinder(dist, min_region_size=6, score_threshold=0.1)
		res = rf.find(gs)
		#printRegions(res)
		self.assertTrue(len(res) == 3)

class test_find_fus_regions_ncut(unittest.TestCase):
	"""find ncut regions in FUS"""
	def test_run(self):
		gs = gelscore.Sequence('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIESVADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFSGNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY')
		#gs = gelscore.Sequence('AAAAAAAAAAAAAYYGSGSGSGSGSAAAAAAAAA')
		dist = gelscore.SimilarityWeight(sim_add=1.0)
		rf = gelscore.NormalizedCutRegionFinder(dist, min_region_size=10, score_threshold=0.025)
		res = rf.find(gs)
		#printRegions(res)
		self.assertTrue(len(res) >= 7)

class test_entropy_quant_min(unittest.TestCase):
	"""entropy quant -- min entropy"""
	def test_run(self):
		gs = gelscore.Sequence('AAAAAAAAAAAAAAAAA')
		eq = gelscore.EntropyQuant()
		res = eq.quant(gs,20)
		self.assertAlmostEqual(res, 0.0)

class test_entropy_quant_max(unittest.TestCase):
	"""entropy quant -- max entropy"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		eq = gelscore.EntropyQuant()
		res = eq.quant(gs,20)
		self.assertAlmostEqual(res,1.0)

class test_iterator(unittest.TestCase):
	"""iterator over sequence"""
	def test_run(self):
		gs = gelscore.Sequence('ACDEFGHIKLMNPQRSTVWY')
		n = 0
		for (i,a) in enumerate(gs):
			n += 1
		self.assertTrue(n == len(gs))

class test_composition1(unittest.TestCase):
	"""sequence composition"""
	def test_run(self):
		gs1 = gelscore.Sequence('AAAA')
		gs2 = gelscore.Sequence('DDDD')
		comp1 = gelscore.SequenceComposition(gs1)
		comp2 = gelscore.SequenceComposition(gs2)
		comp3 = gelscore.SequenceComposition(gs1, weights=[1,0,0,0])
		self.assertAlmostEqual(comp1.dot(comp2),0.0)
		self.assertAlmostEqual(comp1.dot(comp3),1.0)
		self.assertAlmostEqual(np.linalg.norm(comp1.vector), 1)

class test_similar_score(unittest.TestCase):
	"""entropy quant -- max entropy"""
	def test_run(self):
		gs = gelscore.Sequence('GGGGGGGGGGAAAAAAAAAAPPPPPPPPPP')
		sim = gelscore.SimilarityWeight(sim_add=1.0)
		rf = gelscore.NormalizedCutRegionFinder(sim, min_region_size=5, score_threshold=0.8)
		D = sim.weights(gs)
		#print D
		#print D.sum()
		#print D[1,:].sum()
		#print D[:,1].sum()
		res = rf.find(gs)
		#print ""
		#for r in res:
		#	print r.score, r.norm_score, r.region
		self.assertTrue(len(res) == 3)

class test_find_regions2(unittest.TestCase):
	"""find regions in more realistic cases"""
	def test_run(self):
		# DAD: this is not a passing test right now, doesn't split into 3
		gs = gelscore.Sequence('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGACDEFGHIKLMNPQRSTVWYGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG')
		dist = gelscore.SimilarityWeight(sim_add=1.0)
		max_regions = 16
		rf = gelscore.NormalizedCutRegionFinder(dist, min_region_size=6, score_threshold=0.01, max_regions=max_regions)
		res = rf.find(gs)
		printRegions(res)
		self.assertTrue(len(res) == 3)

class test_find_region_homopolymer(unittest.TestCase):
	"""decline to separate homopolymer"""
	def test_run(self):
		gs = gelscore.Sequence('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG')
		#dist = gelscore.SimilarityDistanceWeight(sim_add=1.0, sd=10.0)
		dist = gelscore.SimilarityWeight(sim_add=1.0)
		max_regions = 16
		rf = gelscore.NormalizedCutRegionFinder(dist, min_region_size=6, score_threshold=0.3, max_regions=max_regions)
		res = rf.find(gs)
		#for r in res:
		#	print r.region.start, r.region.end, r.score, r.norm_score, r.region
		self.assertTrue(len(res) == 1)

class test_sim_score(unittest.TestCase):
	"""sim score"""
	def test_run(self):
		gs = gelscore.Sequence('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG')
		n = len(gs)
		#dist = gelscore.SimilarityDistanceWeight(sim_add=1.0, sd=10.0)
		sim = gelscore.SimilarityWeight(sim_add=1.0)
		W = sim.weights(gs)
		score = sim.score(W)
		self.assertAlmostEqual(score, 1.0)


if __name__=="__main__":
	unittest.main(verbosity=2)
