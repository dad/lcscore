#! python

#import sys, os, string
import math, stats
import scipy as sp


# Goal: simulate assembly
# Outputs:
#   list of entities, with number of monomers and size

"""
Basic algorithm.

Start with n monomers of radius r and diameter d
Monomers can assemble to volume fraction f

For each timestep:
	all monomers deactivated
	monomers activate with propensity p_act(T)
	active monomers nucleate into assemblies with propensity p_nuc 
	active monomers assemble onto assemblies with probability given
		by number of monomers on surface
		* p_act(T) ~ surface area (SA) * p_act(T)
	Generate random numbers u
	Choose reaction mu by
		Pr(mu) = p_mu / sum_i p_i, p_i = propensity of reaction i
	Calculate time of next reaction tau = - (1/p_mu) ln u
	Carry out reaction mu
	Update time from t to t + mu


There are several reactions:
- activation
- deactivation
- nucleation
- denucleation
- assembly onto surface

For surface assembly, just need to calculate the total surface area in the reactor
	(over all assemblies),then, if that is the next reaction, choose an actual
	assembly proportional to its surface area, then update the surface area of
	that assembly and update the propensity for surface assembly accordingly.

The number of monomers on the surface
= number of monomers within d of surface
= amount of volume within d of surface * volume fraction of assembled protein /
	volume per monomer (rounded)
= (4/3) pi (R^3 - (R-d)^3) * f / (4/3) pi r^3
= (R^3 - (R-d)^3) * f / r^3, rounded to nearest integer

Below some threshold number of monomers, say Nmin, we just say that all monomers
are on the surface.

"""

N_Avogadro = 6.02214129e23 # Avogadro's constant

#Compute the number of molecules from the volume and the concentration.
def numberOfMolecules(volume_L, concentration_M):
	"""
	This function takes as input a volume and and a concentration and
	outputs the number of molecules.
	"""
	return math.floor(N_Avogadro * concentration_M * volume_L)


'''
Now we will define a series of classes representing key entities: reactants,
reactions, and a reactor. Reactants are the chemical entities, reactions are
processes which turn one set of reactants into another set, and the reactor is
the vessel in which this occurs, allowing global features like volume,
temperature, and so on to be manipulated.
'''
class Reactant(object):
	"""A class representing any chemically reacting species.

	id: a string (e.g. "A" or "monomer") by which the reactant can be
	referenced

	n: an integer representing the initial number of molecules of the
	reactant.
	"""
	def __init__(self, id, n):
		self._id = id
		self._num = int(n)

	
	@property
	def id(self):
		return self._id

	@property
	def number(self):
		return self._num

	'''
	Numbers of reactant molecules can be incremented, decremented, or
	set to a specific value that the user inputs as the reacion proceeds.
	'''
	
	def inc(self,n=1):
		"""
		Increase the number of molecules of the reactant by an integer n.
		"""
		self._num += n

	def dec(self,n=1):
		"""
		Decrease the number of molecules of the reactant by an integer n.
		"""
		self._num -= n

	@number.setter
	def number(self, num):
		"""
		Set the number of molecules of the reactant to a number num.
		"""
		self._num = num

#class TranslatableReactant(Reactant):
#	"""An mRNA."""


class ActivatableReactant(Reactant):
	"""A class representing an activatable reactant.

	Rather than tracking activation directly as a reaction,
	we assume it is so fast that molecules flicker back and forth."""

	'''The proportion of reactants that are activated at a given temperature
	can be obtained using the function getProportionActivated'''
	def getProportionActivated(self, temperature_C):
		"""
		Takes as input a temperature in degrees Celsius and returns
		the proportion of molecules that are activated at that
		temperature.
		"""
		midpoint = 41.0 # degrees C
		cooperativity = 20.0
		#The proportion activated is a sigmoid.
		return 1/(1 + (midpoint/temperature_C)**cooperativity)

class SphericalAssemblyEnsemble(Reactant):
	"""
	A class representing a spherical assembly ensemble.

	Rather than track individual assemblies, we track the overall ensemble and 
	its constituents. One and only one ensemble per reaction.

	"""
	def __init__(self, id, monomers, monomer_radius_nm,
			 packing_fraction=sp.pi/(3*sp.sqrt(2))):
		super(SphericalAssemblyEnsemble, self).__init__(id, n=1)
		self._monomer_radius_nm = monomer_radius_nm
		self._monomers = monomers
		self._packing_fraction = packing_fraction
		self._assembly_list = []
		self._assemblies = {}
		'''self.assemblies is a dictionary of assemblies within the
		ensemble'''
		self._next_assembly_id = 0
		'''self._next_assembly_id is an integer that keeps track of
		assemblies within the ensemble'''
		#self._weighted_assembly_list = []
	
	'''
	The total number of monomers on the surfaces of all of the assemblies
	within the ensemble can be summed and returned.
	'''
	
	def monomers_at_surface(self):
		"""
		The function monomers_at_surface returns the total number of
		monomers on assembly surfaces within the ensemble.
		"""
		return sum([a.monomers_at_surface for a in self.assemblies])

	'''
	It is possible to introduce a new assembly of any size to the ensemble.
	'''
	
	def newAssembly(self, n=2):
		"""
		The function newAssembly takes as input the number of monomers
		that will compose the newAssembly and it creates a new assembly
		of that size.

		n = number of monomers (default to 2)
		"""

		#The id in the assembly dictionary
		id = self._next_assembly_id 
		asmb = SphericalAssembly(id, n, self._monomer_radius_nm,
					 self._packing_fraction) 

		#Increase the index of the next assembly by 1.
		self._next_assembly_id += 1 
		#
		self._assembly_list.append(asmb)
		
		#Place that assembly in the dictionary of assemblies.
		self._assemblies[id] = asmb 
		return asmb

	'''
	Assemblies can be removed from the assembly dictionary.
	'''
	
	def removeAssembly(self, asmb):
		"""
		The function removeAssembly takes as input an assembly asmb and
		removes that assembly from the assembly dictionary.
		"""
		del self._assemblies[asmb.id]
  
	'''
	Monomers can be added to assemblies.
	'''
		
	def inc(self,n=1):
		"""
		The function inc takes as input a number n which is defaulted to 1.
		It adds one monomer to an assembly that it chooses - where the
		probability of being chosen is proportional to the surface area of
		the assembly.
		"""
		# Choose assembly by surface area
		weights = [(a, a.monomers_at_surface) for a in self.assemblies]
		asmb = stats.weighted_choice(weights)
		# Increment
		asmb.inc(n)

	'''
	Monomers can deleted from assemblies.
	'''
	
	def dec(self,n=1):
		"""
		The function dec takes as input a number n which is defaulted to
		1. It deletes one monomer from an assembly that it chooses -
		where the probability is proportional to the surface area of the
		assembly.
		"""
		weights = [(a, a.monomers_at_surface) for a in self.assemblies if
			   a.number >= n]
		asmb = stats.weighted_choice(weights) # Choose assembly by surface area
		# Decrement
		'''diff is the number of monomers in the assembly minus the number
		you want to remove from the assembly'''
		diff = asmb.number-n

		'''diff must be greater than or equal to 0 because you can't
		remove more than the number of monomers that compose the
		assembly.'''
		assert(diff>=0)

		'''If the number you want to remove = the number of monomers in
		the assembly...'''
		if diff == 0: 
			self.removeAssembly(asmb) #Just remove the assembly.
		elif diff == 1: #If there is only one monomer in the assembly,
			#put the monomer back in solution and...
			self.monomers.inc(1) 
			self.removeAssembly(asmb) #just remove the assembly
		else:
			asmb.dec(n) #Otherwise, decrease number of monomers by 1.

	'''
	The largest assembly within the ensemble can be obtained and returned
	using the following function (largest).
	'''
	def largest(self, fn=lambda x: x.radius):
		"""
		The function largest takes as input a function fn(x) which
		computes the radius of a sphere x. The function returns res,
		the assembly a with the largest radius in the
		SphericalAssemblyEnsemble.
		""" 

		res = None

		# If there is at least one assembly in the ensemble
		if self.number>0:
			'''make a list of (radii, assembly) for each assembly
			in the dictionary of assemblies'''
			wlist = [(fn(a),a) for a in self.assemblies]

			# Sort the list from greates to least.
			res = sorted(wlist, reverse=True)[0][1] 
			return res #Return the assembly with the largest radius.
		

	'''
	It is useful to be able to apply an attribute function (such as finding
	the radius) and then a statistics function (i.e. computing the median)
	on the assemblies that compose the assembly ensemble. This next function
	performs these two tasks successively.
	'''
	
	def statistic(self, attr_fxn, stat_fxn):
		"""
		The function statistic takes as input attr_fxn (an attribute
		function) and stat_fxn (a statistics function) and runs the
		attribute function on each assembly in the ensemble and then
		the statistics function on the output of the attribute function.
		"""
		res = None

		# If there is at least one assembly in the ensemble...
		if self.number>0:
			
			''' Then apply the attribute function to each assembly.
				Apply the statistics function to the outpur of the
				attribute function. '''
			res = stat_fxn([attr_fxn(a) for a in self.assemblies])
		return res #Return the output of the statistics function.

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
	"""
	A class representing a spherical assembly.
	"""
	def __init__(self, id, n, monomer_radius_nm,
			 packing_fraction=sp.pi/(3*sp.sqrt(2))):
		super(SphericalAssembly, self).__init__(id, n)
		self._monomer_radius_nm = monomer_radius_nm
		self._packing_fraction = packing_fraction
		self._volume = -1
		self._monomers_at_surface = -1
		self._monomer_volume = (4/3.0)*sp.pi*self._monomer_radius_nm**3.0
		self._updateVolume()

	'''
	The volume of a sphere (i.e. spherical assembly) of a given radius
	can be computed.
	'''
	def _sphereVolume(self, radius):
		"""
		The function sphereVolume takes as input a radius and returns
		the volume of a sphere of that radius using the formula
		volume = (4/3)*pi*(r^3).
		"""
		return 4.1888 * radius*radius*radius


	'''
	The radius of a sphere of a given volume can be computed.
	'''
	
	def _radiusFromVolume(self, volume):
		"""
		The function radiusFromVolume takes as input a volume and computes
		the radius of a sphere from that volume using the formula
		radius = [(3/4)*(v/pi)]^(1/3).
		"""
		return (0.238732414 * volume)**(0.33333333)

	'''
	The number of monomers on the surface of a spherical assembly can be
	computed and returned.
	'''
	
	def _calculateMonomersAtSurface(self):
		"""
		The function calcculateMonomersAtSurface calculates and returns
		the number of monomers on the surface of a spherical assembly.
		"""
		return int(math.floor((self._radius**3.0 -
					   (self._radius-self._monomer_radius_nm*2)
					   **3.0)*self._packing_fraction/
					  (self._monomer_radius_nm**3)))

	'''
	The number of monomers that compose a spherical assembly can be
	incremented by an integer number.
	'''
	
	def inc(self,n=1):
		"""
		The function inc takes as input an integer n and adds n monomers
		to the spherical assembly, and updates its volume accordingly.

		n is an integer number that is set to 1 as a default.
		"""
		self._num += n
		self._updateVolume()


	'''
	The volume of a spherical assembly can be recalculated and updated
	after changing a parameter such as the number of monomers composing the
	assembly.
	'''
	
	def _updateVolume(self):
		"""
		The function updateVolume recalculates and updates the volume of
		a spherical assembly.
		"""
		
		self._volume = (self._monomer_volume * self._num
				/self._packing_fraction)

		'''Get the radius of the SphericalAssembly from the volume
		of the SphericalAssembly.'''
		self._radius = self._radiusFromVolume(self._volume)
		
		#Compute the surface area of the SphericalAssembly.
		self._surface_area = 4 * sp.pi * self._radius * self._radius
		
		#num is the same as the number of monomers at surface
		self._monomers_at_surface = self._num

		'''If there are fewer than 6 monomers, assume all monomers
		are surface-accessible''' 
		if self._num < 6:
			self._monomers_at_surface = self._num

		elif self._num < 12:
			self._monomers_at_surface = self._num-1

		else:
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

	"""
	A class representing a reaction.

	reactants: a dictionary of reactants within the reaction
	reactant_ids: a list of reactant ids (strings) to keep track of the
			 reactants
	products: a dictionary of products within the reaction
	product_ids: a list of product ids (strings) to keep track of the
			products
	"""
	def __init__(self):
		self._reactants = {}
		self._reactant_ids = []
		self._products = {}
		self._product_ids = []


	'''
	Reactants and products can be added to reactions.
	'''
	
	def addReactant(self, id, mult=1):
		"""
		The function addReactant takes as input an id and a mult (a
		multiplicity - the number of a reactant or product involved in a
		given reaction) and it appends the id (the reactant) to the
		reaction's list of reactants. The integer input mult is the
		multiplicity of that reactant.
		"""
 
		self._reactants[id] = mult
		self._reactant_ids.append(id)

	
	def addProduct(self, id, mult=1):
		"""
		The function addProduct takes as input an id and mult
		(a multiplicity - the number of a reactant or product involved
		in a given reaction) and adds that product with its given
		multiplicity to the reaction.
		"""
		self._products[id] = mult
		self._product_ids.append(id)

	def parse(self):
		raise Exception, "not implemented"

	'''
	Reactions have several important qualities, including rate constants
	(which can be set), reactors (which are the vessels in which the reactions
	occur), and propensities (the likelihood of the reaction occurring).
	'''
	@property
	def rate_constant(self):
		return self._rate_constant

	
	@rate_constant.setter
	def rate_constant(self, rate_constant):
		"""
		The function rate_constant takes as input a rate constant and sets
		the Reaction's rate constant to that input.
		"""
		self._rate_constant = rate_constant

	@property
	def reactor(self):
		return self._reactor

	
	@reactor.setter
	def reactor(self, rctr):
		"""
		The function reactor takes as input a reactor (rctr) and sets
		the Reaction's reactor to that input.
		"""
		self._reactor = rctr

	@property
	def propensity(self):
		"""
		Reactions have propensities. The function propensity  outputs
		the propensity for a given reaction to occur.
		"""
		r = self._reactor #r is the reaction's reactor
		nums = [r.getReactant(id).number for id in self._reactant_ids]

		'''The propensity is given by the reaction's rate constant k/V
		(reactor volume)^(number of reactants? - 1)*some function's output
		of nums'''
		prop = (self._rate_constant/(r.volume**(len(nums)-1)))*reduce(
			lambda x,y: x*y, nums)
		return prop

	def react(self):
		raise Exception('Not implemented: must be subclassed')

class DegradationReaction(Reaction):
	"""
	A class representing a degradation reaction.
	"""
	
	def __init__(self, reactant_id, rate_constant):
		super(DegradationReaction, self).__init__()
		self.addReactant(reactant_id)
		self.rate_constant = rate_constant

	
	def react(self):
		"""
		The function react first identifies the DegradationReaction's
		reactor. It sets reactant to the first reactant in the dictionary
		of the reaction's reactants. It deletes that reactant from the
		reaction.
		"""
		reactor = self._reactor
		reactant = reactor.getReactant(self._reactant_ids[0])
		reactant.dec()


class NucleationReaction(Reaction):
	"""
	A class representing a nucleatiokn reaction.
	"""
	
	def __init__(self, monomer_id, assembly_id, rate_constant):
		super(NucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_id = assembly_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant


	'''We can obtain the propensity for a nucleation reaction to occur.'''
	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity for
		a nucleation reaction to occur.
		"""
		r = self._reactor #r is the nucleation reaction's reactor

		#num_monomers is the total number of monomers in the reactor
		num_monomers = r.getReactant(self._monomer_id).number

		'''The propensity is the Nucleation Reaction's rate constant
		(k/V (reactor volume))*(# monomers*(# monomers - 1)/2). This is
		second-order.'''
		prop = (self._rate_constant/r.volume)*(num_monomers*
							   (num_monomers-1)/2) 
		return prop
	  
	def react(self):
		"""
		The function react creates a new assembly composed of 2 monomers.
		"""

		#reactor is the nucleation reaction's reactor
		reactor = self._reactor

		monomer = reactor.getReactant(self._monomer_id)
		
		#Decrease the number of monomers by 2
		monomer.dec(2)
		
		#Create an assembly.
		assembly = reactor.getReactant(self._assembly_id)

		#Increase the number of assemblies by 1.
		assembly.inc(1)


class EnsembleNucleationReaction(Reaction):
	"""
	A class representing an ensemble nucleation reaction.
	"""
	
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(EnsembleNucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant

	'''
	We can find the propensity for an ensemble nucleation reaction to occur.
	'''
	
	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity
		for an ensemble nucleation reaction to occur.
		"""

		#r is the ensemble nucleation reaction's reactor
		r = self._reactor
		
		#Let num_monomers is the number of monomers in the reactor
		num_monomers = r.getReactant(self._monomer_id).number

		#prop is the propensity, which we compute here (second-order)
		prop = (self._rate_constant/r.volume)*(num_monomers*
							   (num_monomers-1)/2)
		return prop


	'''
	Nuclei can be added to the assembly ensemble.
	'''
	
	def react(self):
		"""
		The function react adds a nucleus to the assembly ensemble.
		"""

		#reactor is the ensemble nucleation reaction's reactor
		reactor = self._reactor

		monomer = reactor.getReactant(self._monomer_id)
		
		#Decrease the number of monomers by 2.
		monomer.dec(2)
		
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)

		'''nucleus is a new assembly added to the assembly ensemble -
		composed of 2 monomers.'''

		nucleus = assembly_ensemble.newAssembly(2)
		#assembly_ensemble.addAssembly(nucleus)


class EnsembleDenucleationReaction(Reaction):
	"""
	A class representing an ensemble denucleation reaction.

	An ensemble denucleation reaction involves breaking apart a nucleus (an
	assembly composed of 2 monomers).
	"""
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(EnsembleDenucleationReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		#self.addReactant(monomer_id,2)
		#self.addProduct(assembly_id,1)
		self.rate_constant = rate_constant

	'''
	We can obtain a list of nuclei (assemblies with only 2 monomers) in the
	assembly ensemble.
	'''
	
	def getNuclei(self):
		"""
		The function getNuclei returns a list of nuclei in the assembly
		ensemble.
		"""
		return [asmb for asmb in self._reactor.getReactant
			(self._assembly_ensemble_id).assemblies if asmb.number==2]


	'''
	The propensity for a denucleation reaction to occur can be computed and
	returned.
	'''

	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity for
		an EnsembleDenucleationReaction to occur.
		"""

		# r is the ensemble denucleation reaction's reactor
		r = self._reactor
		
		# num_nuclei is a list of nuclei in the ensemble
		num_nuclei = len(self.getNuclei())

		#prop is the propensity for this reaction to occur (first-order)
		prop = (self._rate_constant/r.volume)*num_nuclei 
		return prop

	def react(self):
		"""
		The function react breaks breaks apart a nucleus, removing that
		assembly (nucleus) from the assembly ensemble, and returning
		the two monomers that composed the nucleus to solution.
		"""
		reactor = self._reactor
		nucleus = self.getNuclei()[0]
		monomer = reactor.getReactant(self._monomer_id)
		monomer.inc(2)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)
		assembly_ensemble.removeAssembly(nucleus)

"""Introduce a new class called SphericalAssemblyReaction of type Reaction.  A SphericalAssemblyReaction has a monomer id, 
   an assembly ensemble id, and a rate constant."""
class SphericalAssemblyReaction(Reaction):
	"""
	A class that represents a spherical assembly reaction.

	A spherical assembly reaction involves adding a monomer to an assembly
	within the assembly ensemble.
	"""
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(SphericalAssemblyReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		self.addReactant(monomer_id,1)
		self.addProduct(assembly_ensemble_id,1)
		self.rate_constant = rate_constant

	'''
	The propensity for a spherical assembly reaction to occur can be
	computed and returned.
	'''
	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity
		for a spherical assembly reaction to occur.
		
		Propensity is given by (k/V) * M_a * S_a, where
		M_a is the number of activated monomers in solution
		S_a is the number of activated monomers on the surface
			of assemblies
		k is the rate constant
		V is the reaction volume
		"""

		#r is the spherical assembly reaction's reactor
		r = self._reactor

		#T is the temperature of the reactor
		T = r.temperature

		monomers = r.getReactant(self._monomer_id)

		#prop_activated is the proportion of monomers that are activated
		prop_activated = monomers.getProportionActivated(T)

		#active_monomers = number of activated monomers in the reactor.
		active_monomers_in_solution = prop_activated * monomers.number
		
		assembly_ensemble = r.getReactant(self._assembly_ensemble_id)

		'''active_monomers_at_surface is the number of active monomers
		on the surface of all assemblies within the assembly ensemble'''
		active_monomers_at_surface = (prop_activated *assembly_ensemble.
						  monomers_at_surface())

		#prop is the propensity for a monomer to be added to the surface
		prop = ((self._rate_constant/r.volume)*active_monomers_in_solution
			*active_monomers_at_surface)
		return prop

	def react(self):
		"""
		The function react adds one monomer to the assembly ensemble.
		"""
		
		 #Let reactor be the SphericalAssemblyReaction's reactor.
		reactor = self._reactor

		#Let monomer be a monomer in the reactor.
		monomer = reactor.getReactant(self._monomer_id)

		#Decrease the number of monomers by 1
		monomer.dec(1)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)

		#Add the monomer to the assembly ensemble.
		assembly_ensemble.inc(1)


class SphericalDisassemblyReaction(Reaction):
	"""
	A class representing a spherical disassembly reaction.

	A spherical disassembly reaction removes a monomer from the assembly
	ensemble and places it back in solution.
	"""
	
	def __init__(self, monomer_id, assembly_ensemble_id, rate_constant):
		super(SphericalDisassemblyReaction, self).__init__()
		self._monomer_id = monomer_id
		self._assembly_ensemble_id = assembly_ensemble_id
		self.addReactant(monomer_id,1)
		self.addProduct(assembly_ensemble_id,1)
		self.rate_constant = rate_constant


	'''
	The propensity for a spherical disassembly reaction to occur can be
	computed and returned.
	'''
	
	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity
		for a spherical disassembly reaction to occur.
		
		Propensity is given by (k) * S_i, where
		S_i is the number of inactivated monomers on the surface of
			assemblies
		k is the rate constant
		"""

		#r is the spherical disassembly reaction's reactor
		r = self._reactor

		#T is the spherical disassembly reaction reactor's temperature
		T = r.temperature

		#monomers is the monomers in the reactor
		monomers = r.getReactant(self._monomer_id)

		#prop_inactivated = proportion of monomers that are inactivated
		prop_inactivated = 1-monomers.getProportionActivated(T)
		
		#active_monomers_in_solution = prop_activated * monomers.number
		assembly_ensemble = r.getReactant(self._assembly_ensemble_id)

		'''inactive monomers at surface is the number of inactivated
		monomers on the surface of assemblies in the ensemble'''
		inactive_monomers_at_surface = prop_inactivated * assembly_ensemble.monomers_at_surface()
		prop = (self._rate_constant)*inactive_monomers_at_surface
		return prop

	def react(self):
		"""
		The function react causes the spherical disassembly reaction to
		occur.
		"""

		#reactor is the spherical disassembly reaction's reactor
		reactor = self._reactor
		
		monomer = reactor.getReactant(self._monomer_id)

		#add one monomer to the monomers in the reactor
		monomer.inc(1)
		assembly_ensemble = reactor.getReactant(self._assembly_ensemble_id)

		#remove a monomer from the assembly ensemble
		assembly_ensemble.dec(1)

class Event(object):

	"""
	A class that represents an event.

	Events are actions such as temperature changes that occur within reactors.
	"""
	def __init__(self):
		self._reactor = None
	 
	@property
	def reactor(self):
		return self._reactor

	'''
	Events are associated with reactors - the vessel in which the event
	occurs. The reactor associated with the event can be set.
	'''
	@reactor.setter
	def reactor(self, rctr):
		"""
		The function reactor takes as input a reactor (rctr) and sets
		the reactor associated with the event to rctr.
		"""
		self._reactor = rctr

"""Introduce a new class called TimedAbsoluteTemperatureChangeEvent of type Event. A TimedAbsoluteTemperatureChangeEvent has 
a time and an absolute temperature."""
class TimedAbsoluteTemperatureChangeEvent(Event):
	"""
	A class representing a timed absolute temperature change event.

	time: the time at which this event occurs
	
	absolute_temperature_C: the new temperature of the reactor that will
		be set at time
	"""
	def __init__(self, time, absolute_temperature_C):
		self._time = time
		self._absolute_temperature_C = absolute_temperature_C

	'''
	A TimedAbsoluteTemperatueChangeEvent can be triggered.
	'''

	def trigger(self):
		"""
		The function trigger sets the temperature of the reactor
		to absolute_temperature_C.
		"""
		if not self._reactor is None:
			self._reactor.temperature = self._absolute_temperature_C
	  
	'''
	The time until the TimedAbsoluteTemperatureChangeEvent takes place can
	be computed and returned.
	'''
	
	def timeFromNow(self, now):
		"""
		The function timeFromNow takes as input a time (now) and computes
		the time until the TimedAbsoluteTemperatureChangeEvent will
		occur.
		"""
		return self._time - now

class TranslationReaction(Reaction):
	"""
	A class representing translation of an mRNA to produce a protein.
	"""
	def __init__(self, mrna_id, protein_id, rate_constant):
		super(Reaction, self).__init__()
		self._mrna_id = mrna_id
		self._protein_id = protein_id
		self.rate_constant = rate_constant

	'''
	The propensity for the translation reaction to occur.
	'''
	@property
	def propensity(self):
		"""
		The function propensity computes and returns the propensity
		for an ensemble nucleation reaction to occur.
		"""
		r = self._reactor
		
		# Number of mRNAs in the reactor
		num_mrnas = r.getReactant(self._mrna_id).number

		#prop is the propensity, which we compute here (first-order)
		prop = (self._rate_constant/r.volume)*num_mrnas
		return prop

	'''
	Translation produces a protein.
	'''
	def react(self):
		"""
		The function react adds a nucleus to the assembly ensemble.
		"""
		reactor = self._reactor
		protein = reactor.getReactant(self._protein_id)
		
		#Translation of the mRNA increases the number of proteins by 1.
		protein.inc(1)


class Reactor(object):

	"""
	A class representing a reactor (the vessel in which everything occurs).

	"""
	
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

	'''
	Events can be added to reactors.
	'''
	
	def addEvent(self, event):
		"""
		The function addEvent takes as input an event and appends that
		event to the reactor's list of events.
		"""
		self._events.append(event)
		event.reactor = self

	'''
	Timed events can be added to reactors.
	'''

	def addTimedEvent(self, timed_event):
		"""
		The function addTimedEvent is a function that takes as input a
		timed_event, appends the timed_event to its list of timed
		events, and sorts these events in the order in which they should
		occur.
		"""

		self._timed_events.append(timed_event)
		now = self._time
		self._timed_events.sort(key=lambda x: x.timeFromNow(now))
		timed_event.reactor = self

	'''
	Reactions and reactants can be added to reactors.
	'''
	
	def addReaction(self, rxn):
		"""
		The function addReaction takes as input a rxn and adds it to the
		list of reactions.
		"""
		self._reactions.append(rxn)
		rxn.reactor = self
	 
	def addReactant(self, reactant):
		"""
		The function addReactant takes as input a reactant and appends it
		to the list of reactants.
		"""
		self._reactants[reactant.id] = reactant
		self._reactant_ids.append(reactant.id)
		self._reactant_ids.sort()

	def getReactant(self, reactant_id):
		"""
		The function getReactant takes as input a reactant_id and returns
		that reactant from the list of reactants.
		"""
		return self._reactants[reactant_id]

	def run(self, time):
		"""
		The function run takes as input time and steps the time as long as
		the reactor's time (start time) is less than that time.
		"""
		start_time = self._time
		while self._time - start_time < time:
			self._step()

	def step(self, n=1):
		"""
		The function step advances the simulation by n steps (default 1).
		"""
		for i in range(n):
			self._step()

	def stepUntil(self, n, time_elapsed):
		"""
		The function stepUntil advances the simulation by n steps or until
		time_elapsed seconds has passed.
		"""
		start_time = self._time
		elapsed = 0.0
		nsteps = 0
		while nsteps < n and elapsed < time_elapsed:
			self._step()
			elapsed = self._time-start_time
			nsteps += 1

	
	def _nextTimedEvent(self):
		"""
		The function nextTimedEvent returns the next timed event that is
		supposed to occur if there are still events left in the reactor.
		"""
		res = None
		if len(self._timed_events)>0:
			res = self._timed_events[0]
		return res

	def triggerTimedEvent(self, event):
		"""
		The function triggerTimedEvent takes as input an event and
		triggers it. It then adds that event to the list of triggered
		events and takes that event away from the list of timed events.
		"""
		event.trigger()
		self._triggered_events.append(event)
		self._timed_events.pop()
	 
	
	def _step(self):
		"""
		The function step causes the next reaction to react, and advances
		the time and number of steps.
		"""
		
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
		"""
		The function exhausted checks whether a reaction is exhausted.
		"""
		return self._exhausted

	@property
	def volume(self):
		"""
		The function volume returns the reactor's volume in microliters.
		"""
		return self._volume_uL

	@property
	def temperature(self):
		"""
		The function temperature returns the reactor's temperature in
		Celsius.
    	"""
		return self._temperature_C

	@temperature.setter
	def temperature(self, T):
		"""
		The function temperature takes as input a temperature T and sets
		the reactor's temperature to that input.
		"""
		self._temperature_C = T

	def _weighted_choice(self, choices, total, r):
		"""
		The function weighted_choice takes as input choices, total, and r
		and returns the first c for which the index upto+c's corresponding
		w >= the input r.
		"""
		upto = 0
		for c, w in choices:
			if upto + w >= r:
				return c
			upto += w

	def __getitem__(self, key):
		return self._reactants[key]

	'''
	The next reaction to occur within the reactor and the time at which this
	reaction occurs can be chosen.
	'''
	def chooseNextReactionAndTime(self):
		"""
		The function chooseNextReactionAndTime computes the propensities
		for each reaction to occur and creates a list of these
		propensities.  It returns the next reaction and the time
		at which that reaction will occur.
		"""
		propensities = [r.propensity for r in self._reactions]
		#print propensities
		weighted_propensities = zip(self._reactions,propensities)
		total_prop = sum(propensities)
		tau = 0.0
		rxn = None
		if total_prop>0.0:
			u = sp.random.uniform(0.0, total_prop)
			'''Reaction is chosen with probability proportional to
			propensity'''
			rxn = self._weighted_choice(weighted_propensities,
							total_prop, u)
			'''Time to next reaction is exponential with rate
			total_prop'''
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

