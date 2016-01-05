import scipy as sp
import math
import biofile, translate

def makeAlignment(nrow, ncol, average_gap, average_hamming, gap='-'):
	chars = [c for c in translate.AAs()]
	gaps = [gap]
	alts = [chars, gaps]
	alignment = []
	for x in range(nrow):
		which_kind = sp.random.binomial(1, average_gap, ncol)
		st = ''.join([sp.random.choice(alts[xi]) for xi in which_kind])
		alignment.append(st)
	return alignment

def alignmentToBinaryMatrix(a, gap='-'):
	x = sp.matrix([sp.r_[[int(s[xi]==gap) for xi in range(len(s))]] for s in a])
	return x

def hamming(num1, num2):
	assert len(num1) == len(num2)
	return sp.sum(num1 != num2)

def isGap(x):
	return x == '-'

def entropy(thelist, gap='-', base=20.0):
	counts = [(thelist.count(aa),aa) for aa in translate.AAs()]
	ungapped_len = len(thelist)-thelist.count(gap)
	res = 0.0
	if ungapped_len>0:
		props = [ct/float(ungapped_len) for (ct,aa) in counts if ct>0]
		res = sum([-p*math.log(p,base) for p in props])
	return res

def mostFrequent(thelist, nogap=False, gap='-'):
	counts = sorted([(thelist.count(aa),aa) for aa in [x for x in translate.AAs()]+[gap]], reverse=True)
	res = counts[0][1]
	if nogap and res == gap and len(counts)>1 and counts[1][0]>0:
		res = counts[1][1]
	return res


class ProfileAlignment(object):
	def __init__(self, alignment):
		self._align = alignment
		# Convert alignment to binary matrix
		self._matrix = alignmentToBinaryMatrix(alignment)
		#print self._matrix
		self._L = len(alignment[0])

	def binary(self, i):
		return self._matrix[:,i]

	def column(self, i):
		return [s[i] for s in self._align]

	def distances(self):
		# Iterator
		prevec = self._matrix[:,0]
		for v in self.vectors():
			yield hamming(v, prevec)
			prevec = v

	def vectors(self):
		for i in range(self._L):
			vec = self._matrix[:,i]
			yield vec

	def columns(self):
		for xi in range(self._L):
			yield [s[xi] for s in self._align]

	def entropies(self):
		# Iterator
		for col in self.columns():
			yield entropy(col)

