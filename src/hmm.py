# From https://github.com/phvu/misc/tree/master/viterbi
import numpy as np

'''
N: number of hidden states
'''

class ForwardBackwardResult(object):
	def __init__(self, p_fwd, p_bkw, fwd, bkw, posterior):
		self.p_fwd = p_fwd
		self.p_bkw = p_bkw
		self.fwd = fwd
		self.bkw = bkw
		self.posterior = posterior

class HMM(object):
	def __init__(self, initial_prob, trans_prob, obs_prob):
		# Number of hidden states
		self.N = initial_prob.shape[0]
		# Initial probabilities
		self.initial_prob = initial_prob
		# Transition probabilities
		self.trans_prob = trans_prob
		# Observation probabilities
		self.obs_prob = obs_prob
		assert self.initial_prob.shape == (self.N, 1)
		assert self.trans_prob.shape == (self.N, self.N)
		assert self.obs_prob.shape[0] == self.N
		
	def Obs(self, obs):
		return self.obs_prob[:, obs, None]

	def bestPath(self, obs):
		# Return maximum-likelihood path through the hidden states, given a set of observations
		trellis = np.zeros((self.N, len(obs)))
		backpt = np.ones((self.N, len(obs)), 'int32') * -1
				
		# initialization
		trellis[:, 0] = np.squeeze(self.initial_prob * self.Obs(obs[0]))
				
		for t in xrange(1, len(obs)):
			trellis[:, t] = (trellis[:, t-1, None].dot(self.Obs(obs[t]).T) * self.trans_prob).max(0)
			backpt[:, t] = (np.tile(trellis[:, t-1, None], [1, self.N]) * self.trans_prob).argmax(0)
		# termination
		tokens = [trellis[:, -1].argmax()]
		for i in xrange(len(obs)-1, 0, -1):
			tokens.append(backpt[tokens[-1], i])
		return tokens[::-1]

	# From http://en.wikipedia.org/wiki/Forward-backward_algorithm
	def forwardBackward(self, obs): #, states, a_0, a, e, end_st):
		L = len(obs)
		a_0 = self.initial_prob
		a = self.trans_prob
		e = self.obs_prob
		states = range(self.N)
		end_st = obs[-1]
	 
		fwd = []
		f_prev = {}
		# forward part of the algorithm
		for i, obsi in enumerate(obs):
			f_curr = {}
			for st in states:
				if i == 0:
					# base case for the forward part
					prev_f_sum = a_0[st]
				else:
					prev_f_sum = sum(f_prev[k]*a[k][st] for k in states)
	 
				f_curr[st] = e[st][obsi] * prev_f_sum
	 
			fwd.append(f_curr)
			f_prev = f_curr
	 
		p_fwd = sum(f_curr[k]*a[k][end_st] for k in states)
	 
		bkw = []
		b_prev = {}
		# backward part of the algorithm
		for i, obsi_plus in enumerate(reversed(obs[1:]+[None])):
			b_curr = {}
			for st in states:
				if i == 0:
					# base case for backward part
					b_curr[st] = a[st][end_st]
				else:
					b_curr[st] = sum(a[st][l]*e[l][obsi_plus]*b_prev[l] for l in states)
	 
			bkw.insert(0,b_curr)
			b_prev = b_curr
	 
		p_bkw = sum(a_0[l] * e[l][obs[0]] * b_curr[l] for l in states)
	 
		# merging the two parts
		posterior = []
		for i in range(L):
			posterior.append({st: fwd[i][st]*bkw[i][st]/p_fwd for st in states})
	 
		assert p_fwd == p_bkw
		return ForwardBackwardResult(p_fwd, p_bkw, fwd, bkw, posterior)