#! python

import sys, os, math, random, argparse
import util, aggsim, stats, na
import scipy as sp


def intensityWeightedRadius(n_monomers, monomer_radius, assembly_radii):
	# <R_h>_z = (sum_i n_i r_i^6)/(sum_i n_i r_i^5)
	num = sum([n_monomers * monomer_radius**6.0] + [r**6.0 for r in assembly_radii])
	den = sum([n_monomers * monomer_radius**5.0] + [r**5.0 for r in assembly_radii])
	return num/den

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest='initial_concentration_uM', type=float, default=None, help="initial concentration of monomers (uM)")
	parser.add_argument(dest='volume_fL', type=float, default=None, help="volume of reactor (fL)")
	# Optional arguments
	parser.add_argument("-T", "--temperature", dest='initial_temperature_C', type=float, default=30, help='initial temperature (Celsius)')
	parser.add_argument("--k-nucleation", dest='k_nucleation', type=float, default=1.0, help='rate constant for nucleation [1/(uM sec)]')
	parser.add_argument("--k-denucleation", dest='k_denucleation', type=float, default=0.0, help='rate constant for denucleation [1/(uM sec)]')
	parser.add_argument("--k-assembly", dest='k_assembly', type=float, default=1.0, help='rate constant for assembly [1/(uM sec)]')
	parser.add_argument("--k-disassembly", dest='k_disassembly', type=float, default=1.0, help='rate constant for disassembly [1/(uM sec)]')
	parser.add_argument("-r", "--monomer-radius", dest='monomer_radius_nm', type=float, default=1.0, help='radius of monomer (nm)')
	parser.add_argument("--jump-temp", dest='jump_temperature_C', type=float, default=42, help='jump temperature (Celsius)')
	parser.add_argument("--jump-time", dest='jump_time', type=float, default=None, help='time of temperature jump (sec)')
	parser.add_argument("--output-step", dest='output_step_number', type=int, default=1, help='max number of steps between output')
	parser.add_argument("--output-time", dest='output_step_time', type=int, default=1, help='max number of seconds between output')
	parser.add_argument("-t", "--time", dest="simulation_time", type=float, default=1e4, help="time, in seconds, to simulate")
	parser.add_argument("-s", "--seed", dest="random_seed", type=float, default=111, help="seed for random number generator")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Set seed
	sp.random.seed(options.random_seed)

	# Set up reaction
	reactor = aggsim.Reactor(volume_uL=1e-9*options.volume_fL, temperature_C=options.initial_temperature_C)
	n_monomers = aggsim.numberOfMolecules(1e-15 * options.volume_fL, 1e-6 * options.initial_concentration_uM)
	monomers = aggsim.ActivatableReactant('monomer', n_monomers)
	reactor.addReactant(monomers)
	#print monomer.getProportionActivated(reactor.temperature)
	ensemble = aggsim.SphericalAssemblyEnsemble('ensemble', monomers, options.monomer_radius_nm, packing_fraction=sp.pi/(3*sp.sqrt(2)))
	reactor.addReactant(ensemble)
	# Nucleation
	nuc_reaction = aggsim.EnsembleNucleationReaction('monomer','ensemble', options.k_nucleation * 1e-6) # convert rate to (uM s)^-1
	reactor.addReaction(nuc_reaction)
	# Nucleation
	denuc_reaction = aggsim.EnsembleDenucleationReaction('monomer','ensemble', options.k_denucleation * 1e-6) # convert rate to (uM s)^-1
	reactor.addReaction(denuc_reaction)
	# Assembly
	asmb_reaction = aggsim.SphericalAssemblyReaction('monomer','ensemble', options.k_assembly * 1e-6) # convert rate to (uM s)^-1
	reactor.addReaction(asmb_reaction)
	# Disassembly
	disasmb_reaction = aggsim.SphericalDisassemblyReaction('monomer','ensemble', options.k_disassembly * 1e-6) # convert rate to (uM s)^-1
	reactor.addReaction(disasmb_reaction)

	# DAD: implement temp changes as timed reactions. Let reactor pick when to execute them.
	# E.g. next reaction at tau, but timed change at change_time - cur_time < tau, then execute change instead of reaction.
	if not options.jump_time is None:
		tempchange = aggsim.TimedAbsoluteTemperatureChangeEvent(options.jump_time, options.jump_temperature_C) # convert rate to (uM s)^-1
		reactor.addTimedEvent(tempchange)

	# Write output headers
	dout = util.DelimitedOutput()
	dout.addHeader('step','Simulation step','d')
	dout.addHeader('time','Time', 'f')
	dout.addHeader('temperature','Temperature (Celsius)', 'f')
	dout.addHeader('prop.activated','Proportion of activated monomers', 'f')
	dout.addHeader('n.monomers', 'Number of monomers', 'd')
	dout.addHeader('prop.assembled', 'Proportion of monomers in assemblies', 'f')
	dout.addHeader('n.assemblies', 'Number of assemblies', 'd')
	dout.addHeader('largest.radius', 'Radius of largest assembly (nm)', 'f')
	dout.addHeader('mean.radius', 'Mean particle radius (nm)', 'f')
	dout.addHeader('intensity.weighted.radius', 'Intensity-weighted assembly radius (nm) as expected in dynamic light-scattering experiments', 'f')
	dout.addHeader('numbers.in.assemblies', 'Numbers of monomers in each assembly', 's')
	dout.describeHeader(data_outs)

	def computeAndWrite(reactor):
		"""The per-step work of reading reactor state and computing figures of interest is done here."""
		ens = reactor['ensemble']
		mon = reactor['monomer']
		numbers = [a.number for a in ens.assemblies]
		numbers.sort(reverse=True)
		if len(numbers)==0:
			sizestr = '0'
		else:
			sizestr = ','.join(['{:d}'.format(n) for n in numbers])
		largest_assembly = ens.largest()
		if largest_assembly is None:
			largest_radius = ens.monomer_radius
		else:
			largest_radius = largest_assembly.radius
		mean_radius = options.monomer_radius_nm
		if len(numbers)>0:
			mean_assembly_radius = ens.statistic(lambda x: x.radius, stats.mean)
			mean_radius = (mon.number * options.monomer_radius_nm + mean_assembly_radius * ens.number)/float(mon.number + ens.number)
		intwt_radius = intensityWeightedRadius(mon.number, options.monomer_radius_nm, [a.radius for a in ens.assemblies])

		prop_active = mon.getProportionActivated(reactor.temperature)
		prop_assembled = 1.0 - float(mon.number)/n_monomers

		line = ("{r.num_steps:d}\t{time:1.2f}\t{r.temperature}\t{prop_active:1.1E}\t{m.number:d}\t{prop_assembled:1.1E}\t{a.number:d}\t" + \
			"{largest_radius:1.2f}\t{mean_radius:s}\t{intwt_radius:s}\t{sizes:s}\n").format(
			r=reactor, time=min(options.simulation_time, reactor.time), prop_active=prop_active, m=mon, prop_assembled=prop_assembled, a=ens, 
			largest_radius=largest_radius, mean_radius=na.formatNA(mean_radius, '{:1.2f}'), 
			intwt_radius=na.formatNA(intwt_radius, '{:1.2f}'), sizes=sizestr)
		data_outs.write(line)

	# Main loop
	dout.writeHeader(data_outs)
	n_written = 0
	jumped = False
	do_jump = not options.jump_time is None
	while (reactor.time < options.simulation_time) and not reactor.exhausted:
		# Compute figures, write out
		computeAndWrite(reactor)
		n_written += 1
		#if do_jump and reactor.time > options.jump_time:
		#	# Change temperature and capture the change
		#	reactor.temperature = options.jump_temperature_C
		#	computeAndWrite(reactor)
		#	n_written += 1
		# Advance the simulation
		reactor.stepUntil(options.output_step_number, options.output_step_time)
	# Final step
	computeAndWrite(reactor)
	n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

