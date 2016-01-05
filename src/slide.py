#! python

import sys, os, math, random, argparse
import util, protprop, translate, stats
import scipy as sp, numpy as np
import gelscore as gs

# Goal: assess repetitive composition
# Approach: pick composition and window size, slide window over target sequence and record similarity
# Similarity is assessed by comparator objects which compute features of sequences. Sequences are not compared directly.

#prop = protprop.ProteinProperties()

class Comparator(object):
	def __init__(self):
		pass
	
	def compare(self, comp):
		raise Exception("Not implemented")
	def score(self, comp):
		raise Exception("Not implemented")

class SequenceCompositionComparatorFactory(object):
	def __init__(self):
		self._prop = protprop.ProteinProperties()
	
	def make(self, x, source):
		if source =='composition':
			return self.fromComposition(x)
		elif source == 'sequence':
			return self.fromSequence(x)
		else:
			raise ValueError, "Sequence comparator source {} not known; use one of ['composition','sequence']".format(source)
	
	def fromSequence(self, seq):
		aas = translate.AAs()
		composition = self._prop.getComposition(seq, normalize=True, aas=aas)
		return SequenceCompositionComparator(composition)

	def fromComposition(self, comp):
		# Verify normalization
		tot = float(sum([c for (aa,c) in comp]))
		composition = [c/tot for (aa,c) in comp]
		return SequenceCompositionComparator(composition)

# Simplest composition comparison
class SequenceCompositionComparator(Comparator):
	def __init__(self, comp):
		# composition
		self._comp = comp
		self._vec = sp.r_[[f for (aa,f) in self._comp]]
	
	@property
	def composition(self):
		return self._comp
	
	def composition_vector(self):
		return self._vec
	
	def compare(self, scc):
		# Return a score between 0.0 (most dissimilar) and 1.0 (most similar)
		# For normalized histograms, the maximum chiSq distance = 1.0, minimum 0.0
		dist = stats.chiSquaredHistogramDistance(self._vec, scc._vec)
		res = 1.0 - dist
		return res

	def score(self, scc):
		# Return a score between 0.0 (most dissimilar) and 1.0 (most similar)
		# For normalized histograms, the maximum chiSq distance = 1.0, minimum 0.0
		dist = stats.chiSquaredHistogramDistance(self._vec, scc._vec)
		res = 1.0 - dist
		return res

# Simplest composition scorer: fraction of given amino acids in window
class SequenceCompositionScorer(Comparator):
	def __init__(self, aas):
		self._comp = comp
		self._vec = sp.r_[[f for (aa,f) in self._comp]]
	
	@property
	def composition(self):
		return self._comp
	
	def composition_vector(self):
		return self._vec
	
	def compare(self, scc):
		# Return a score between 0.0 (most dissimilar) and 1.0 (most similar)
		# For normalized histograms, the maximum chiSq distance = 1.0, minimum 0.0
		dist = stats.chiSquaredHistogramDistance(self._vec, scc._vec)
		res = 1.0 - dist
		return res

class SearchResult(object):
	def __init__(self, position, score, window, range):
		self.position = position # 1-based position
		self.total_score = score
		self.window = window
		self.sequence = window.currentSequence()
		self.range = range
		self.average_score = float(score)/range

class SlideResult(object):
	def __init__(self, position, score, aa, window):
		self.position = position # 1-based position
		self.score = score
		self.aa = aa
		self.window = window

class SlideResultSet(object):
	def __init__(self, seq, positions, scores, sequence_window):
		self.seq = seq
		self.positions = positions
		self.scores = scores
		self.window = sequence_window
		self.n = len(positions)
		assert(len(positions)==len(scores))
	
	def results(self):
		for xi in range(self.n):
			# Positions are 1-based
			pos = self.positions[xi]
			yield SlideResult(pos, self.scores[xi], self.seq[pos-1], self.seq[(pos-1):((pos-1)+self.window.size)])

class SequenceWindow(object):
	def __init__(self, seq, size, startpos=1):
		"""Create new sequence window. startpos is a 1-based sequence position (first amino acid, startpos = 1)"""
		self._seq = seq
		self._startpos = startpos
		self._curind = startpos-1
		assert(self._curind>=0 and self._curind<len(seq))
		self._size = min(size,len(seq)-self._curind)
		# DAD: Check size
		assert(self._curind+self._size <= len(seq))
		self._comparator_factory = SequenceCompositionComparatorFactory()
	
	def next(self, amt=1):
		"""Move to next position, advancing by amt"""
		newind = self._curind + int(amt)
		self._curind = min(newind, len(self._seq))
	
	def maxNumberOfWindows(self):
		# How many windows are there?
		return len(self._seq)-(self._size-1)
	
	@property
	def position(self):
		# Positions are 1-based (e.g. residue 3 is in position 3)
		return self._curind+1
	
	@property
	def index(self):
		return self._curind
	
	@property
	def size(self):
		return self._size
	
	def isValid(self):
		return self._curind >= 0 and self._curind+self._size <= len(self._seq)
	
	def reset(self):
		self._curind = self._startpos-1
	
	def currentSequence(self):
		newend = min(self._curind+self._size,len(self._seq))
		subseq = self._seq[self._curind:newend]
		return subseq
	
	def comparator(self):
		# Return a comparision object
		subseq = self.currentSequence()
		comp = self._comparator_factory.make(subseq, 'sequence')
		return comp
	
	def __str__(self):
		return str(self.currentSequence())
	
	def slide(self, scorer):
		positions = []
		scores = []
		while self.isValid():
			mycomp = self.comparator()
			score = scorer.score(self.currentSequence())
			#sim = mycomp.compare(scorer)
			positions.append(self.position)
			scores.append(score)
			self.next()
		res = SlideResultSet(self._seq, positions, scores, self)
		return res
	
	def search(self, end_pos=None):
		"""Find windows of specified size that maximize score"""
		# pos = position = 1-based; ind = index = 0-based
		end_ind = None
		if end_pos is None:
			end_ind = len(self._seq)-self.size
			end_pos = len(self._seq)-self.size+1
		else:
			end_ind = end_pos-1
		# Slave window will generate comparators for each position
		slave_window = SequenceWindow(self._seq, self._size, self._startpos)
		max_score = 0.0
		max_ind = 0
		while slave_window.isValid() and slave_window.position<=end_pos:
			comp = slave_window.comparator()
			while self.isValid() and self._curind<=end_ind:
				res = self.slide(comp)
				total_score = sum(res.scores)
				if total_score > max_score:
					max_score = total_score
					max_ind = slave_window.index
				self.next()
			# Start search from the beginning again
			self.reset()
			slave_window.next()
		reswindow = SequenceWindow(self._seq, self._size, max_ind+1)
		res = SearchResult(position=max_ind+1, score=max_score, window=reswindow, range=end_pos-self._startpos)
		return res
			
			
		


