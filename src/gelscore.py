#! python

import sys, os, math, random, argparse
import util, translate
from numpy.lib.scimath import logn
import numpy as np
import scipy as sp
import ncut

class Sequence(object):
	# A sequence
	def __init__(self, seq, startindex=None, endindex=None):
		if startindex is None:
			self._seq = seq
		else:
			if endindex is None:
				self._seq = seq[startindex:]
				sefl._endindex = len(seq)
			else:
				assert startindex < endindex
				self._seq = seq[startindex:endindex]
	
	def count(self, x):
		return self._seq.count(x)
	
	def __len__(self):
		return len(self._seq)
	
	def __getitem__(self, it):
		return self._seq[it]
	
	def __str__(self):
		return str(self._seq)
	
	

class Comparator(object):
	def __init__(self, quant):
		pass

def entropy(num_list, base=2):
	assert base > 0
	ent = 0.0
	for p in num_list:
		if p>0.0:
			ent += -p*logn(base, p)
	return ent

class Quant(object):
	"""A class which takes a sequence and produces a numeric vector."""
	def __init__(self):
		self._aas = translate.AAs()
	
	def quant(self, seq):
		self._seq = seq
		return None
	
class CountQuant(Quant):
	"""A class which takes a sequence and produces a scalar number."""
	def __init__(self):
		super(CountQuant, self).__init__()
	
	def quant(self, seq):
		# count
		self._seq = seq
		counts = {}
		length = 0
		for aa in self._aas:
			counts[aa] = self._seq.count(aa)
		return [counts[aa] for aa in self._aas]
	
class AlphabetSizeQuant(CountQuant):
	"""A class which takes a sequence and produces a scalar number indicating the number of different amino acids used."""
	def __init__(self):
		super(AlphabetSizeQuant, self).__init__()
	
	def quant(self, seq):
		# count
		self._seq = seq
		counts = super(AlphabetSizeQuant).quant()
		return len([x for x in counts if x>0])
	

class EntropyQuant(Quant):
	"""A class which takes a sequence and produces a scalar number."""
	def __init__(self):
		super(EntropyQuant, self).__init__()
	
	def quant(self, seq, base=2.0):
		# count
		sorted_seq = sorted(seq)
		chars = sorted(list(set(sorted_seq)))
		cur_char = chars[0]
		char_counts = [0]*len(chars)
		char_ind = 0
		for i in xrange(len(sorted_seq)):
			aa = sorted_seq[i]
			if cur_char != aa:
				cur_char = aa
				char_ind += 1
			char_counts[char_ind] += 1
		length = float(len(seq))
		# Entropy
		score = entropy([v/length for v in char_counts], base=base)
		return score

def noop(x):
	"""No operation"""
	return x

class SimilarityWeight(object):
	def __init__(self, sim_add=1.0, max_dist=None):
		self._sim_add = float(sim_add)
		if max_dist is None:
			self._infer_max_dist = True
		else:
			self._infer_max_dist = False
			self._max_dist = max_dist
	
	def similarity(self, seq, i, j):
		if seq[i]==seq[j]:
			mult = 1.0
			if abs(j-i)>self._max_dist:
				mult = 0.0
			return 1.0*mult
		return 0.0
	
	def weights(self, seq): # return a weight matrix
		n = len(seq)
		D = np.zeros(shape=(n,n))
		sim = self.similarity
		simadd = self._sim_add
		# If we are guessing max distance, set it
		# to be the length of the sequence
		if self._infer_max_dist:
			self._max_dist = n

		for i in xrange(n):
			v = D[i,:]*0.0
			for j in xrange(0,n):
				v[j] += simadd*sim(seq, i, j)
			v[i] = 0.0
			D[i,:] = v
		return D
	
	def score(self, dist):
		n = dist.shape[0]
		res = 0.0
		# n*(n-1)/2 is the maximum score (all entries maximized)
		# in the upper triangle. 
		if n>1:
			res = np.triu(dist).sum()/(self._sim_add*n*(n-1)/2.0)
		return res

class DistanceSimilarityWeight(object):
	def __init__(self, sim_add=1.0, distance_falloff=1.0, max_distance=20, similarity_fxn=lambda x,y: x==y):
		self._sim_fxn = similarity_fxn
		self._max_dist = max_distance
		self._sim_add = sim_add
		self._dist_falloff = distance_falloff
	
	def similarity(self, x, y):
		if x==y:
			return 1.0
		return 0.0
	
	def weights(self, seq): # return a weight matrix
		n = len(seq)
		D = np.zeros(shape=(n,n))
		dfall = self._dist_falloff
		sim = self._sim_fxn
		maxdist = self._max_dist
		for i in xrange(n):
			v = D[i,:]*0.0
			for j in xrange(0,n):
				dist_mult = sp.exp(-(j-i)*(j-i)/10)
				v[j] += self._sim_add * self.similarity(seq[i],seq[j]) * dist_mult
			v[i] = 0.0
			D[i,:] = v
			D[:,i] = v
			#tot = v.sum()
			#if tot > 0.0:
			#	D[i,:] = v/tot # Ensures that each row sums to 1.0
		# Now make symmetric
		#Dsim = (D + D.transpose())/2.0
		#return Dsim
		return D

class NormalizedDistanceWeight(object):
	def __init__(self, sim_add=1.0, distance_falloff=1.0, max_distance=20, similarity_fxn=lambda x,y: x==y):
		self._sim_fxn = similarity_fxn
		self._max_dist = max_distance
		self._sim_add = sim_add
		self._dist_falloff = distance_falloff
	
	def weights(self, seq): # return a weight matrix
		n = len(seq)
		D = np.zeros(shape=(n,n))
		dfall = self._dist_falloff
		sim = self._sim_fxn
		maxdist = self._max_dist
		for i in xrange(0, n-1, 1):
			v = D[i,:]*0.0
			for j in xrange(i+1, min(i+1+maxdist,n), 1):
				v[j] = sp.exp(-dfall*(j-i))
				if sim(seq[i],seq[j]): # if similar, give bump to weights
					v[j] += self._sim_add
		v[i] = 0.0
		tot = v.sum()
		if tot > 0.0:
			D[i,:] = v/tot # Ensures that each row sums to 1.0
		else:
			D[i,:] = v
		# Now make symmetric
		Dsim = (D + D.transpose())/2.0
		return D

class GaussianDistanceWeight(object):
	def __init__(self, sim_add=1.0, variance=1.0, max_distance=20, similarity_fxn=lambda x,y: x==y):
		self._sim_fxn = similarity_fxn
		self._max_dist = int(max_distance)
		self._sim_add = float(sim_add)
		self._variance = float(variance)
	
	def weights(self, seq): # return a weight matrix
		n = len(seq)
		D = np.zeros(shape=(n,n))
		sim = self._sim_fxn
		maxdist = self._max_dist
		for i in xrange(0, n-1, 1):
			D[i,i] = 1.0
			for j in xrange(i+1, min(i+1+maxdist,n), 1):
				D[i,j] = sp.exp(-(j-i)*(j-i)/self._variance)
				if sim(seq[i],seq[j]): # if similar, give bump to weights
					D[i,j] += self._sim_add
				D[j,i] = D[i,j]
		D[n-1,n-1] = 1.0
		return D

class SimilarityComparator(object):
	def similarity(self, obj1, obj2):
		return None

class SequenceCompositionSimilarity(SimilarityComparator):
	def similarity(self, obj1, obj2):
		return obj1.composition.similarity(obj2.composition)

class SequenceComposition(object):
	def __init__(self, seq, weights=None, alphabet=translate.AAs()):
		self._seq = seq
		self._len = len(seq)
		if not weights is None:
			assert len(seq) == len(weights)
		else:
			weights = [1.0]*len(seq)
		d = dict([(aa,0) for aa in alphabet])
		for (i, aa) in enumerate(seq):
			try:
				d[aa] += weights[i]
			except KeyError:
				d[aa] = weights[i]
		self._alphabet = alphabet
		self._composition_vec = sp.r_[[float(d[a]) for a in alphabet]]
		self.normalize()
	
	def similarity(self, comp2):
		return self.dot(comp2)
	
	def normalize(self):
		"""Normalize vector to length 1."""
		if self._len > 0:
			norm = np.linalg.norm(self._composition_vec)
			if norm > 0.0:
				self._composition_vec = self._composition_vec / norm
	
	def dot(self, seqcomp):
		res = 0.0
		if self._len > 0:
			res = sp.dot(self._composition_vec, seqcomp._composition_vec)
		return res
	
	def __getitem__(self, it):
		return self._composition_vec[it]
	
	def __str__(self):
		return ",".join(["{}:{}".format(a,v) for (a,v) in zip(self._alphabet, self._composition_vec)])
	
	@property
	def vector(self):
		return self._composition_vec

class ZComposition(SequenceComposition):
	"""Composition quantified by Z scores from provided """
	def __init__(self, seq, distributions, alphabet=translate.AAs()):
		self._seq = seq
		self._alphabet = alphabet
		p = dict([(aa,0.0) for aa in alphabet])
		for (i, aa) in enumerate(seq):
			p[aa] += 1.0
		n = len(seq)
		self._proportion_vec = sp.r_[[float(p[a])/n for a in alphabet]]
		# 
		comp = []
		for a in alphabet:
			d = distributions[a]
			z = (p[a] - d.mean)/d.sd
			comp.append(z)
		self._composition_vec = sp.r_[comp]
		self.normalize()

class ContiguousRegion(Sequence):
	def __init__(self, seq, start=None, end=None):
		self._seq = seq
		#self._max_gap_size = max_gap_size
		self._start = start
		self._end = end
		if not start is None and not end is None:
			assert start <= end
		if start is None:
			self._start = 0
		if end is None:
			self._end = len(seq)
		self._updateComposition()
	
	def _updateComposition(self):
		self._composition = SequenceComposition(self._seq[self._start:self._end], weights=None)
	
	def collides(self, region):
		if self._seq != region._seq:
			return False
		r1 = self
		r2 = region
		if r2.start < r1.start:
			r1 = region
			r2 = self
		return r1.end > r2.start
	
	def avoids(self, region):
		return not self.collides(region)
	
	def copy(self):
		return ContiguousRegion(self._seq, self._start, self._end)
	
	def merge(self, region):
		assert self._seq == region._seq
		
		# DAD: move this to scored region
		#assert self._threshold == region._threshold
		# DAD: check for start/end compatibility?
		merged_region = ContiguousRegion(self._seq, min(self._start, region._start), max(self._end, region._end))
		return merged_region
	
	def interdiff(self, region):
		"""Return difference and intersection regions"""
		left = self
		right = region
		if left.start > right.start:
			right = self
			left = region
		if not self.collides(region):
			# No overlap; just return original regions
			return [left, None, right]
		else:
			# Regions overlap
			left = self
			right = region
			if left.start > right.start:
				right = self
				left = region
			diffleft = ContiguousRegion(left._seq, left.start, right.start)
			if left.end > right.end: # right is a subset of left
				diffright = ContiguousRegion(left._seq, right.end, left.end)
				inter = right
			else: # right is not a subset of left
				diffright = ContiguousRegion(right._seq, left.end, right.end)
				inter = ContiguousRegion(right._seq, right.start, left.end)
			
			return [diffleft, inter, diffright]
	
	def difference(self, region):
		"""Return two regions"""
		[diffleft, inter, diffright] = self.interdiff(region)
		return [diffleft, diffright]
	
	def intersection(self, region):
		[diffleft, inter, diffright] = self.interdiff(region)
		return inter
	
	def distance(self, region):
		res = None
		if self.collides(region):
			res = 0
		else:
			if self._seq == region._seq:
				left = self
				right = region
				if left.start > right.start:
					right = self
					left = region
				res = right.start-left.end
		return res
	
	def isSubset(self, query_region):
		"""Determines whether this region is a subset of the argument (query) region"""
		return self.start >= query_region.start and self.end <= query_region.end
	
	def avoids(self, region):
		return not self.collides(region)
	
	def trimright(self, n):
		self._end = max(0, self._end-n)
		self._end = max(self._start, self._end)
		self._updateComposition()
		

	def trimleft(self, n):
		self._start = min(len(self), self._end+n)
		self._start = min(self._start, self._end)
		self._updateComposition()
	
	@property
	def sequence(self):
		return self._seq[self._start:self._end]
	
	@property
	def composition(self):
		return self._composition
	
	@property
	def start(self):
		return self._start
	
	@property
	def end(self):
		return self._end
	
	def __getitem__(self, it):
		return self._seq[self._start:self._end][it]

	def __len__(self):
		return self._end-self._start
	
	def __str__(self):
		return self._seq[self._start:self._end]


def normScore(n, score):
	res = 0.0
	if n>0:
		#res = sp.sqrt(score)/n
		res = score/(n*n)
	return res

def unnormScore(n, score):
	#return score*score*n
	return score*n*n

class RegionFinderResult(object):
	def __init__(self, region, score):
		self.region = region
		self.score = score

class BinaryCutNode(object):
	def __init__(self):
		self._left_child = None
		self._right_child = None
		self._parent = None
		self._size = 0
		# Data-specific
		self.region = None
		self.score = None
	
	def add_left(self, child):
		self._left_child = child
		self._left_child._parent = self

	def add_right(self, child):
		self._right_child = child
		self._right_child._parent = self
	
	def traversal(self):
		if not self._left_child is None:
			yield self._left_child.descendants()
		if not self._right_child is None:
			yield self._right_child.descendants()
		yield self

class BinaryCutTree(object):
	def __init__(self, root_node):
		self._size = 0
		self._root = root_node
	
	def depth_traverse(self):
		return self._root.traversal()

def pullRegionAndScore(node):
	return node.region, node.score

def accumulateTreeResults(tree):
	regions = []
	scores = []
	for n in tree.depth_traverse():
		reg = n.region
		if not reg is None:
			regions.append(n.region)
			scores.append(n.score)
	return regions, scores
	

# Premature optimization is the root of all evil.
class NormalizedCutRegionFinder(object):
	def __init__(self, aa_similarity_fxn, min_region_size, score_threshold, max_regions=None):
		self._seq = None
		self._min_region_size = min_region_size
		self._aa_sim = aa_similarity_fxn
		#self._merge_similarity = merge_similarity # Merge regions with similarity >= this value (-1,1)
		self._score_threshold = score_threshold
		self._max_regions = max_regions
	
	def find(self, seq):
		# Find minimum normalized cuts
		# From Shi and Malik 2000:
		# Given graph G = (V,E) (vertices, edges) where each vertex is a site,
		# and all sites are connected to all other sites with edge weights w(v_i, v_j).
		# V alone thus fully describes the graph.
		# a cut divides the graph into disjoint vertex sets A and B.
		# Let cut(A,B) = \sum_{u \in A, v \in B} w(u,v), the total weight of edges that are cut
		# Let assoc(A,V) = \sum_{u \in A, t \in V} w(u,t), the total weight of edges connected to vertices in A
		# Define Ncut(A,B) = cut(A,B)/assoc(A,V) + cut(A,B)/assoc(B,V), the normalized cut
		# 
		# We wish to find k cuts which minimize Ncut. This is NP-hard but complete enumeration is not difficult
		# with 1D sequences of typical protein lengths.
		
		# Compute similarity matrix
		D = self._aa_sim.weights(seq)
		n = len(seq)
		max_regions = min(n-1, self._max_regions)
		
		cut_tree_root = BinaryCutNode()
		cut_tree = BinaryCutTree(cut_tree_root)

		self._findRecursive(seq, D, 0, n, cut_tree_root)
		
		# Pull results out of tree
		regions, scores = accumulateTreeResults(cut_tree)
		
		#score_vector = sp.r_[[y[0] for y in sorted(scores, key=lambda x: x[1])]]
		assert len(regions) == len(scores)
		# DAD: limit regions?
		if not self._max_regions is None:
			rs = zip(regions, scores)
			rs.sort(key=lambda x: x[1], reverse=True)
			regions, scores = zip(*(rs[0:self._max_regions]))
		
		res = [RegionFinderResult(r,s) for (r,s) in zip(regions, scores)]
		return res
		
	def _findRecursive(self, seq, distance_matrix, start, end, cut_node):
		n = end-start
		#print seq[start:end], n
		D = distance_matrix

		debug = False
		
		# Initialize values
		max_score = 2.0
		scores = [(max_score,-1)]*n # Scores[i] for cutting before residue i
		scores[0] = (max_score,start,-1,-1,-1)
		# Prepare to find minimum score and index at which it occurs
		min_index = -1
		min_score = max_score
		selfsim_score = self._aa_sim.score(D[start:end,start:end])

		for i in range(start+1,end):
			# Cut before i
			#dAA = D[start:i,start:i]
			#dBB = D[i:end,i:end]
			dAB = D[start:i,i:end]
			dAV = D[start:i,start:end]
			dBV = D[i:end,start:end]
			#ssAA = dAA.sum()    # Self-similarity of A
			#ssBB = dBB.sum()    # Self-similarity of A
			cutAB = dAB.sum()   # Weight of edges to be cut
			assocAV = dAV.sum() # Weight of all edges from A
			assocBV = dBV.sum() # Weight of all edges from B
			ncutAB = max_score
			fcutA = fcutB = 0.0
			if assocAV>0.0:
				fcutA = cutAB/assocAV
			if assocBV>0.0:
				fcutB = cutAB/assocBV
			# The normalized cut
			ncutAB = fcutA + fcutB
			scores[i-start] = (ncutAB, i, cutAB, assocAV, assocBV)
			# Store score and location if it's a potential minimum
			if ncutAB < min_score or (ncutAB==min_score and i>min_index):
				min_score = ncutAB
				min_index = i
				#min_left = normScore(i-start, dAA.sum())
				#min_right = normScore(end-i, dBB.sum())

		if debug:
			min_ind = sorted(scores)[0][1]
			print "---", len(seq), n
			for i in range(n):
				(ncutAB, ind, cutAB, assocAV, assocBV) = scores[i]
				print i+start, seq[i+start], "{:1.4f}\t{:1.4f}\t{:1.4f}\t{:1.4f}\t{:1.4f}\t{:1.4f}".format(
					ncutAB, cutAB, assocAV, assocBV, cutAB/assocAV, cutAB/assocBV),
				if i+start == min_ind:
					print "****"
				else:
					print ""

		if n > self._min_region_size:
			# Binary search each side
			regions_left = []
			regions_right = []
			scores_left = []
			scores_right = []
			noleft = True
			noright = True
			# Recursion condition: if left or right regions are large enough to be split
			if min_index-start > self._min_region_size:
				leftnode = BinaryCutNode()
				cut_node.add_left(leftnode)
				self._findRecursive(seq, D, start, min_index, leftnode)
				noleft = len(regions_left)==0
			if end-min_index > self._min_region_size:
				rightnode = BinaryCutNode()
				cut_node.add_right(rightnode)
				self._findRecursive(seq, D, min_index, end, rightnode)
				noright = len(regions_right)==0
			if noleft and noright:
				# No satisfactory subregions found; return just this region, unbroken
				cut_node.region = ContiguousRegion(seq, start=start, end=end)
				cut_node.score = selfsim_score
		else:
			# No point in splitting further
			cut_node.region = ContiguousRegion(seq, start=start, end=end)
			cut_node.score = selfsim_score
