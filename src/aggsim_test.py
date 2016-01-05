#! python

import os, sys, math, unittest
import aggsim
import scipy as sp
# from [library] import [function]

class test001(unittest.TestCase):
	def test_deg(self):
		"""degradation"""
		reactor = aggsim.Reactor()
		A = aggsim.Reactant('A', 1000)
		Deg = aggsim.DegradationReaction('A', 2.0)
		reactor.addReactant(A)
		reactor.addReaction(Deg)
		while not reactor.exhausted:
			reactor.step()
		self.assertTrue(reactor['A'].number == 0)

	def test_steps(self):
		"""advancing n steps"""
		reactor = aggsim.Reactor()
		A = aggsim.Reactant('A', 1000)
		Deg = aggsim.DegradationReaction('A', 2.0)
		reactor.addReactant(A)
		reactor.addReaction(Deg)
		reactor.step(10)
		self.assertTrue(reactor.num_steps == 10)

class test002(unittest.TestCase):
	def test_even(self):
		"""nucleation, even monomers"""
		reactor = aggsim.Reactor()
		n = 1000
		mon = aggsim.Reactant('monomer', n)
		nuc = aggsim.Reactant('nucleus', 0)
		nucr = aggsim.NucleationReaction('monomer','nucleus',2.0)
		reactor.addReactant(mon)
		reactor.addReactant(nuc)
		reactor.addReaction(nucr)
		while not reactor.exhausted:
			reactor.step()
		self.assertTrue(reactor['nucleus'].number == n/2)
		self.assertTrue(reactor['monomer'].number == 0)

	def test_odd(self):
		"""nucleation, odd number of monomers"""
		reactor = aggsim.Reactor()
		n = 1001
		mon = aggsim.Reactant('monomer', n)
		nuc = aggsim.Reactant('nucleus', 0)
		nucr = aggsim.NucleationReaction('monomer','nucleus',2.0)
		reactor.addReactant(mon)
		reactor.addReactant(nuc)
		reactor.addReaction(nucr)
		while not reactor.exhausted:
			reactor.step()
		self.assertTrue(reactor['nucleus'].number == (n-1)/2)
		self.assertTrue(reactor['monomer'].number == 1)

	def test_asmb(self):
		"""nucleation, into assembly"""
		reactor = aggsim.Reactor()
		n = 1001
		mon = aggsim.Reactant('monomer', n)
		nuc = aggsim.SphericalAssembly('nucleus', 0, 1.0, 0.74)
		nucr = aggsim.NucleationReaction('monomer','nucleus',2.0)
		reactor.addReactant(mon)
		reactor.addReactant(nuc)
		reactor.addReaction(nucr)
		while not reactor.exhausted:
			reactor.step()
		self.assertTrue(reactor['nucleus'].number == (n-1)/2)
		self.assertTrue(reactor['monomer'].number == 1)

	def test_ens(self):
		"""nucleation, into ensemble"""
		reactor = aggsim.Reactor()
		n = 1001
		mon = aggsim.Reactant('monomer', n)
		nuc = aggsim.SphericalAssemblyEnsemble('nucleus', 1.0, 0.74)
		nucr = aggsim.EnsembleNucleationReaction('monomer','nucleus',2.0)
		reactor.addReactant(mon)
		reactor.addReactant(nuc)
		reactor.addReaction(nucr)
		while not reactor.exhausted:
			#print reactor['monomer'].number, reactor['nucleus'].number
			reactor.step()
		self.assertTrue(reactor['nucleus'].number == (n-1)/2)
		self.assertTrue(reactor['monomer'].number == 1)

class test004(unittest.TestCase):
	def test_single_assembly(self):
		"""monomers at surface, single assembly"""
		reactor = aggsim.Reactor()
		n = 9
		mon = aggsim.Reactant('monomer', n)
		nuc = aggsim.Reactant('nucleus', 0)
		amb = aggsim.SphericalAssembly('assembly', n, 1.0)
		self.assertTrue(amb.monomers_at_surface == n-1)
		amb.inc(10)
		self.assertTrue(amb.monomers_at_surface > n-1)

class test005(unittest.TestCase):
	def test_assembly_by_area(self):
		"""add monomers to largest assembly"""
		reactor = aggsim.Reactor()
		large_n = 10000
		#mon = aggsim.Reactant('monomer', n)
		ens = aggsim.SphericalAssemblyEnsemble('nucleus', 1.0, 0.74)
		small = ens.newAssembly(1)
		large = ens.newAssembly(large_n)
		# Add a single monomer
		ens.inc(1)
		for a in ens.assemblies:
			#print a.id, a.number, a.monomers_at_surface
			if a.id == large.id:
				self.assertTrue(a.number == large_n+1)
			elif a.id == small.id:
				self.assertTrue(a.number == 1)

class test006(unittest.TestCase):
	def setUp(self):
		self.reactor = aggsim.Reactor()
		self.large_n = 10000
		self.ens = aggsim.SphericalAssemblyEnsemble('nucleus', 1.0, 0.74)

	def test_largest_assembly(self):
		"""retrieve largest assembly"""
		ens = self.ens
		small = ens.newAssembly(1)
		large = ens.newAssembly(self.large_n)
		a = ens.largest()
		self.assertTrue(a.number == self.large_n)
		a = ens.largest(fn=lambda x: x.number)
		self.assertTrue(a.number == self.large_n)

	def test_del(self):
		"""delete largest assembly"""
		a = self.ens.newAssembly(500)
		self.ens.removeAssembly(a)
		self.assertTrue(True)

class test007(unittest.TestCase):
	def setUp(self):
		reactor = aggsim.Reactor(volume_uL=1e-9*45, temperature_C=40)
		self.reactor = reactor
		n_monomers = aggsim.numberOfMolecules(1e-15 * 45, 1e-6 * 10)
		monomers = aggsim.ActivatableReactant('monomer', n_monomers)
		self.monomer = monomers
		reactor.addReactant(monomers)
		#print monomer.getProportionActivated(reactor.temperature)
		ensemble = aggsim.SphericalAssemblyEnsemble('ensemble', monomers, 1.0, packing_fraction=sp.pi/(3*sp.sqrt(2)))
		self.ens = ensemble
		reactor.addReactant(ensemble)
		# Nucleation
		nuc_reaction = aggsim.EnsembleNucleationReaction('monomer','ensemble', 1e-13 * 1e-6) # convert rate to (uM s)^-1
		reactor.addReaction(nuc_reaction)
		# Nucleation
		denuc_reaction = aggsim.EnsembleDenucleationReaction('monomer','ensemble', 1e-3 * 1e-6) # convert rate to (uM s)^-1
		reactor.addReaction(denuc_reaction)
		# Assembly
		asmb_reaction = aggsim.SphericalAssemblyReaction('monomer','ensemble', 1e-8 * 1e-6) # convert rate to (uM s)^-1
		reactor.addReaction(asmb_reaction)
		# Disassembly
		disasmb_reaction = aggsim.SphericalDisassemblyReaction('monomer','ensemble', 1e-6 * 1e-6) # convert rate to (uM s)^-1
		reactor.addReaction(disasmb_reaction)

	def total_molecules(self):
		return self.monomer.number + sum([a.number for a in self.ens.assemblies])

	def test_conservation(self):
		tot = self.total_molecules()
		for xi in xrange(1000):
			self.reactor.step(1)
			#print self.monomer.number, self.ens.number
			self.assertTrue(self.total_molecules()==tot)

# Regulation
class test008(unittest.TestCase):
	def setUp(self):
		self.n_initial_mrnas = 10
		self.n_initial_prots = 0
		reactor = aggsim.Reactor(volume_uL=1e-9*45, temperature_C=25)
		self.reactor = reactor
		n_prot_A = self.n_initial_prots #aggsim.numberOfMolecules(1e-15 * 45, 1e-6 * 10)
		prot_A = aggsim.ActivatableReactant('prot-A', n_prot_A)
		reactor.addReactant(prot_A)
		n_mrna_A = self.n_initial_mrnas #aggsim.numberOfMolecules(1e-15 * 45, 1e-6 * 10)
		mrna_A = aggsim.Reactant('mrna-A', n_mrna_A)
		#self.prot_A = monomers
		reactor.addReactant(mrna_A)
		transl_reaction = aggsim.TranslationReaction('mrna-A','prot-A', 1e-7)
		reactor.addReaction(transl_reaction)

	def test_production(self):
		reactor = self.reactor
		for xi in xrange(1000):
			#print(reactor.time, reactor['mrna-A'].number, reactor['prot-A'].number)
			#self.assertTrue(reactor['monomer'].number == 1)
			self.reactor.step(1)
		self.assertTrue(reactor['mrna-A'].number==self.n_initial_mrnas)
		self.assertTrue(reactor['prot-A'].number>self.n_initial_prots)

# Repressible translation
class test008(unittest.TestCase):
	def setUp(self):
		self.n_initial_mrnas = 10
		self.n_initial_prots = 0
		reactor = aggsim.Reactor(volume_uL=1e-9*45, temperature_C=25)
		self.reactor = reactor
		n_prot_A = self.n_initial_prots #aggsim.numberOfMolecules(1e-15 * 45, 1e-6 * 10)
		prot_A = aggsim.ActivatableReactant('prot-A', n_prot_A)
		reactor.addReactant(prot_A)
		n_mrna_A = self.n_initial_mrnas #aggsim.numberOfMolecules(1e-15 * 45, 1e-6 * 10)
		mrna_A = aggsim.Reactant('mrna-A', n_mrna_A)
		#self.prot_A = monomers
		reactor.addReactant(mrna_A)
		transl_reaction = aggsim.TranslationReaction('mrna-A','prot-A', 1e-7)
		reactor.addReaction(transl_reaction)

	def test_production(self):
		reactor = self.reactor
		for xi in xrange(1000):
			#print(reactor.time, reactor['mrna-A'].number, reactor['prot-A'].number)
			#self.assertTrue(reactor['monomer'].number == 1)
			self.reactor.step(1)
		self.assertTrue(reactor['mrna-A'].number==self.n_initial_mrnas)
		self.assertTrue(reactor['prot-A'].number>self.n_initial_prots)



if __name__=="__main__":
	unittest.main(verbosity=2)
