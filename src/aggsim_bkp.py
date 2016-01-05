#! python

#import sys, os, string
import math, stats
import scipy as sp

# Goal: simulate assembly
# Outputs:
#   list of entities, with number of monomers and size

'''
Basic algorithm.

Start with n monomers of radius r and diameter d
Monomers can assemble to volume fraction f

For each timestep:
	all monomers deactivated
	monomers activate with propensity p_act(T)
	active monomers nucleate into assemblies with propensity p_nuc 
	active monomers assemble onto assemblies with probability given by number of monomers on surface * p_act(T) ~ surface area (SA) * p_act(T)

	Generate random numbers u
	Choose reaction mu by Pr(mu) = p_mu / sum_i p_i, p_i = propensity of reaction i
	Calculate time of next reaction tau = - (1/p_mu) ln u
	Carry out reaction mu
	Update time from t to t + mu


There are four reactions:
- activation
- deactivation
- nucleation
- assembly onto surface

For surface assembly, just need to calculate the total surface area in the reactor (over all assemblies),
then, if that is the next reaction, choose an actual assembly proportional to its surface area, then update
the surface area of that assembly and update the propensity for surface assembly accordingly.

The number of monomers on the surface
= number of monomers within d of surface
= amount of volume within d of surface * volume fraction of assembled protein / volume per monomer (rounded)
= (4/3) pi (R^3 - (R-d)^3) * f / (4/3) pi r^3
= (R^3 - (R-d)^3) * f / r^3, rounded to nearest integer

Below some threshold number of monomers, say Nmin, we just say that all monomers are on the surface.

'''


N_Avogadro = 6.02214129e23 # Avogadro's constant

def numberOfMolecules(volume_L, concentration_M):
	"""Number of molecules in a given volume at a given concentration"""
	return math.floor(N_Avogadro * concentration_M * volume_L) # molec/mole * mole/L * L



class Reactant(object):
	def __init__(self, id, n):
		self._id = id
		self._num = int(n)

	@property
	def id(self):
		return self._id

	@property
	def number(self):
		return self._num

	def inc(self,n=1):
		self._num += n

	def dec(self,n=1):
		self._num -= n

	@number.setter
	def number(self, num):
		self._num = num

class ActivatableReactant(Reactant):
	"""Rather than tracking activation directly as a reaction, 
	we assume it is so fast that molecules flicker back and forth."""
	def getProportionActivated(self, temperature_C):
		midpoint = 41.0 # degrees C
		cooperativity = 20.0
		return 1/(1 + (midpoint/temperature_C)**cooperativity)

class SphericalAssemblyEnsemble(Reactant):
	"""Rather than track individual assemblies, we track the overall ensemble and 
	its constituents. One and only one ensemble per reaction."""
	def __init__(self, id, monomers, monomer_radius_nm, packing_fraction=sp.pi/(3*sp.sqrt(2))):
		super(SphericalAssemblyEnsemble, self).__init__(id, n=1)
		self._monomer_radius_nm = monomer_radius_nm
		self._monomers = monomers
		self._packing_fraction = packing_fraction
		#self._assembly_list = []
		self._assemblies = {}
		self._next_assembly_id = 0
		#self._weighted_assembly_list = []

	def monomers_at_surface(self):
		return sum([a.monomers_at_surface for a in self.assemblies])

	def newAssembly(self, n=2):
		id = self._next_assembly_id
		asmb = SphericalAssembly(id, n, self._monomer_radius_nm, self._packing_fraction)
		self._next_assembly_id += 1
		#self._assembly_list.append(asmb)
		self._assemblies[id] = asmb
		return asmb

	def removeAssembly(self, asmb):
		del self._assemblies[asmb.id]

	def inc(self,n=1):
		# Choose assembly by surface area
		weights = [(a, a.monomers_at_surface) for a in self.assemblies]
		asmb = stats.weighted_choice(weights)
		# Increment
		asmb.inc(n)

	def dec(self,n=1):
		# Choose assembly by surface area
		weights = [(a, a.monomers_at_surface) for a in self.assemblies if a.number >= n]
		asmb = stats.weighted_choice(weights)
		# Decrement
		diff = asmb.number-n
		assert(diff>=0)
		if diff == 0:
			self.removeAssembly(asmb)
		elif diff == 1:
			# Return to monomer pool
			self.monomers.inc(1)
			self.removeAssembly(asmb)
		else:
			asmb.dec(n)

	def largest(self, fn=lambda x: x.radius):
		res = None
		if self.number>0:
			wlist = [(fn(a),a) for a in self.assemblies]
			res = sorted(wlist, reverse=True)[0][1]
		return res

	def statistic(self, attr_fxn, stat_fxn):
		res = None
		if self.number>0:
			res = stat_fxn([attr_fxn(a) for a in self.assemblies])
		return res

	def __getitem__(self, id):
		return self._assemblies[id]

	@property
	def monomer_radius(self):
		return self._monomer_radius_nm
	@property
	def assemblies(self):
		for a in self._assemblies.values():
			yield a

	@property
	def number(self):
		return len(self._assemblies)


class SphericalAssembly(Reactant):
	def __init__(self, id, n, monomer_radius_nm, packing_fraction=sp.pi/(3*sp.sqrt(2))):
		super(SphericalAssembly, self).__init__(id, n)
		self._monomer_radius_nm = monomer_radius_nm
		self._packing_fraction = packing_fraction
		self._volume = -1
		self._monomers_at_surface = -1
		self._monomer_volume = (4/3.0)*sp.pi*self._monomer_radius_nm**3.0
		self._updateVolume()

	def _sphereVolume(self, radius):
		# (4/3)*pi * r^3
		return 4.1888 * radius*radius*radius

	def _radiusFromVolume(self, volume):
		# r = [(3/4)*(v/pi)]^(1/3)
		return (0.238732414 * volume)**(0.33333333)

	def _calculateMonomersAtSurface(self):
		return int(math.floor((self._radius**3.0 - (self._radius-self._monomer_radius_nm*2)**3.0)*self._packing_fraction/(self._monomer_radius_nm**3)))

	def inc(self,n=1):
		self._num += n
		self._updateVolume()

	def _updateVolume(self):
		self._volume = self._monomer_volume * self._num / self._packing_fraction
		self._radius = self._radiusFromVolume(self._volume)
		self._surface_area = 4 * sp.pi * self._radius * self._radius
		self._monomers_at_surface = self._num
		if self._num < 6: # If fewer than given number of molecules, assume all monomers are surface-accessible
			self._monomers_at_surface = self._num
		elif self._num < 12: # If fewer than given number of molecules, assume all but one monomer is accessible
			self._monomers_at_surface = self._num-1
		else:
			# Actually do a calculation
			self._monomers_at_surface = self._calculateMonomersAtSurface()

	@property
	def volume(self):
		return self._volume

	@property
	def radius(self):
		return self._radius

	@property
	def surface_area(self):
		return self._surface_area

	@property
	def monomers_at_surface(self):
		return self._monomers_at_surface

class Reaction(object):
	def __init__(self):
		self._reactants = {}
		self._reactant_ids = []
		self._products = {}
		self._product_ids = []

	def addReactant(self, id, mult=1):
		self._reactants[id] = mult
		self._reactant_ids.append(id)

	def addProduct(self, id, mult=1):
		self._products[id] = mult
		self._product_ids.append(id)

	def parse(self):
		raise Exception, "not implemented"

	@property
	def rate_constant(self):
		return self._rate_constant

	@rate_constant.setter
	def rate_constant(self, rate_constant):
		self._rate_constant = rate_constant

	@property
	def reactor(self):
		return self._reactor

	@reactor.setter
	def reactor(self, rctr):
		self._reactor = rctr

	@property
	def propensity(self):
		r = self._reactor
		nums = [r.getReactant(id).number for id in self._reactant_ids]
		prop = (self._rate_constant/(r.volume**(len(nums)-1)))*reduce(lambda x,y: x*y, nums)
		return prop

	def react(self):
		raise Exception('Not implemented: must be subclassed')

class DegradationReaction(Reaction):
	def __init__(self, reactant_id, rate_constant):
		super(DegradationReaction, self).__init__()
		self.addReactant(reactant_id)
		self.rate_constant = rate_constant

	def react(self):
		reactor = self._reactor
		reactant = reactor.getReactant(self._reactant_ids[0])
		reactant.dec()

class NucleationReaction(Reaction):
	def __init__(self, monomer_id, assembly_id, rate_constant):
		super(NucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_id = assembly_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant

	@property
	def propensity(self):
		r = self._reactor
		num_monomers = r.getReactant(self._monomer_id).number
		prop = (self._rate_constant/r.volume)*(num_monomers*(num_monomers-1)/2) # second-order
		return prop

	def react(self):
		reactor = self._reactor
		monomer = reactor.getReactant(self._monomer_id)
		monomer.dec(2)
		assembly = reactor.getReactant(self._assembly_id)
		assembly.inc(1)

class EnsembleNucleationReaction(Reaction):
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(EnsembleNucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant

	@property
	def propensity(self):
		r = self._reactor
		num_monomers = r.getReactant(self._monomer_id).number
		prop = (self._rate_constant/r.volume)*(num_monomers*(num_monomers-1)/2) # second-order
		return prop

	def react(self):
		reactor = self._reactor
		monomer = reactor.getReactant(self._monomer_id)
		monomer.dec(2)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)
		nucleus = assembly_ensemble.newAssembly(2)
		#assembly_ensemble.addAssembly(nucleus)

class EnsembleDenucleationReaction(Reaction):
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(EnsembleDenucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant

	def getNuclei(self):
		return [asmb for asmb in self._reactor.getReactant(self._assembly_ensemble_id).assemblies if asmb.number==2]

	@property
	def propensity(self):
		r = self._reactor
		# Get all nuclei
		num_nuclei = len(self.getNuclei())
		prop = (self._rate_constant/r.volume)*num_nuclei # first-order
		return prop

	def react(self):
		reactor = self._reactor
		nucleus = self.getNuclei()[0]
		monomer = reactor.getReactant(self._monomer_id)
		monomer.inc(2)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)
		assembly_ensemble.removeAssembly(nucleus)

class SphericalAssemblyReaction(Reaction):
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(SphericalAssemblyReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		self.addReactant(monomer_id,1)
		self.addProduct(assembly_ensemble_id,1)
		self.rate_constant = rate_constant

	@property
	def propensity(self):
		"""Propensity is given by (k/V) * M_a * S_a, where
			M_a is the number of activated monomers in solution
			S_a is the number of activated monomers on the surface of assemblies
			k is the rate constant
			V is the reaction volume
			"""
		r = self._reactor
		T = r.temperature
		monomers = r.getReactant(self._monomer_id)
		prop_activated = monomers.getProportionActivated(T)
		active_monomers_in_solution = prop_activated * monomers.number
		assembly_ensemble = r.getReactant(self._assembly_ensemble_id)
		active_monomers_at_surface = prop_activated * assembly_ensemble.monomers_at_surface()
		prop = (self._rate_constant/r.volume)*active_monomers_in_solution*active_monomers_at_surface
		return prop

	def react(self):
		reactor = self._reactor
		monomer = reactor.getReactant(self._monomer_id)
		monomer.dec(1)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)
		assembly_ensemble.inc(1)

class SphericalDisassemblyReaction(Reaction):
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(SphericalDisassemblyReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		self.addReactant(monomer_id,1)
		self.addProduct(assembly_ensemble_id,1)
		self.rate_constant = rate_constant

	@property
	def propensity(self):
		"""Propensity is given by (k) * S_i, where
			S_i is the number of inactivated monomers on the surface of assemblies
			k is the rate constant
			"""
		r = self._reactor
		T = r.temperature
		monomers = r.getReactant(self._monomer_id)
		prop_inactivated = 1-monomers.getProportionActivated(T)
		#active_monomers_in_solution = prop_activated * monomers.number
		assembly_ensemble = r.getReactant(self._assembly_ensemble_id)
		inactive_monomers_at_surface = prop_inactivated * assembly_ensemble.monomers_at_surface()
		prop = (self._rate_constant)*inactive_monomers_at_surface
		return prop

	def react(self):
		reactor = self._reactor
		monomer = reactor.getReactant(self._monomer_id)
		monomer.inc(1)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)
		assembly_ensemble.dec(1)


class Event(object):
	def __init__(self):
		self._reactor = None

	@property
	def reactor(self):
		return self._reactor

	@reactor.setter
	def reactor(self, rctr):
		self._reactor = rctr

class TimedAbsoluteTemperatureChangeEvent(Event):
	def __init__(self, time, absolute_temperature_C):
		self._time = time
		self._absolute_temperature_C = absolute_temperature_C

	def trigger(self):
		if not self._reactor is None:
			self._reactor.temperature = self._absolute_temperature_C

	def timeFromNow(self, now):
		return self._time - now


class Reactor(object):
	def __init__(self, volume_uL=1.0, temperature_C=25, logstream=None):
		# Actual time elapsed (sec)
		self._time = 0.0
		# Number of simulation steps executed
		self._num_steps = 0
		self._reactions = []
		# Dictionary of reactants, indexed by their ID
		self._reactants = {}
		# String IDs denoting each reactant type
		self._reactant_ids = []
		# Timed events
		self._timed_events = []
		self._triggered_events = []
		# Volume in uL
		self._volume_uL = volume_uL
		# Temperature in Celsius
		self._temperature_C = temperature_C
		self._exhausted = False
		if not logstream is None:
			self._logger = logstream

	def addEvent(self, event):
		self._events.append(event)
		event.reactor = self

	def addTimedEvent(self, timed_event):
		self._timed_events.append(timed_event)
		now = self._time
		self._timed_events.sort(key=lambda x: x.timeFromNow(now))
		timed_event.reactor = self

	def addReaction(self, rxn):
		self._reactions.append(rxn)
		rxn.reactor = self

	def addReactant(self, reactant):
		self._reactants[reactant.id] = reactant
		self._reactant_ids.append(reactant.id)
		self._reactant_ids.sort()

	def getReactant(self, reactant_id):
		return self._reactants[reactant_id]

	def run(self, time):
		start_time = self._time
		while self._time - start_time < time:
			self._step()

	def step(self, n=1):
		"""Advance simulation by n steps (default 1)"""
		for i in range(n):
			self._step()

	def stepUntil(self, n, time_elapsed):
		"""Advance simulation by n steps or until time_elapsed seconds has passed"""
		start_time = self._time
		elapsed = 0.0
		nsteps = 0
		while nsteps < n and elapsed < time_elapsed:
			self._step()
			elapsed = self._time-start_time
			nsteps += 1

	def _nextTimedEvent(self):
		res = None
		if len(self._timed_events)>0:
			res = self._timed_events[0]
		return res

	def triggerTimedEvent(self, event):
		event.trigger()
		self._triggered_events.append(event)
		self._timed_events.pop()

	def _step(self):
		# Choose next reaction
		(rxn, tau) = self.chooseNextReactionAndTime()
		next_event = self._nextTimedEvent()
		if not next_event is None:
			event_tau = next_event.timeFromNow(self._time)
			if event_tau < tau:
				# Do timed event
				self.triggerTimedEvent(next_event)
				tau = event_tau
		else:
			if not self.exhausted:
				# Execute reaction
				rxn.react()
		# Update time
		self._time += tau
		self._num_steps += 1

	@property
	def exhausted(self):
		"""Is reaction exhausted?"""
		return self._exhausted

	@property
	def volume(self):
		return self._volume_uL

	@property
	def temperature(self):
		return self._temperature_C

	@temperature.setter
	def temperature(self, T):
		self._temperature_C = T

	def _weighted_choice(self, choices, total, r):
		upto = 0
		for c, w in choices:
			if upto + w >= r:
				return c
			upto += w

	def __getitem__(self, key):
		return self._reactants[key]

	def chooseNextReactionAndTime(self):
		propensities = [r.propensity for r in self._reactions]
		#print propensities
		weighted_propensities = zip(self._reactions,propensities)
		total_prop = sum(propensities)
		tau = 0.0
		rxn = None
		if total_prop>0.0:
			u = sp.random.uniform(0.0, total_prop)
			# Reaction is chosen with probability proportional to propensity
			rxn = self._weighted_choice(weighted_propensities, total_prop, u)
			# Time to next reaction is exponential with rate total_prop
			r = sp.random.uniform()
			tau = -(1.0/total_prop) * sp.log(r)
		else:
			# Propensities are zero -- nothing left to do
			self._exhausted = True
		return (rxn, tau)

	@property
	def time(self):
		return self._time

	@property
	def num_steps(self):
		return self._num_steps

	def write(self, stream):
		reactant_str = '\t'.join(['{:d}'.format(self[k].number) for k in self._reactant_ids])
		stream.write("{step}\t{time}\t{reactants:s}\n".format(step=self._num_steps, time=self._time, reactants=reactant_str))

