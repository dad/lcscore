#! python

import sys, os, math, random, argparse
import util, biofile, muscle, translate

'''
Prune out:
1) 
Choose 
'''

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prune alignment to a single species and maximum percent identity")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")
	parser.add_argument(dest="in_filter_fname", help="input species-to-keep filename")
	# Optional arguments
	parser.add_argument("--gap-thresh", dest="gap_threshold", type=float, default=0.20, help="maximum proportion gaps to contribute to consensus")
	parser.add_argument("--id-thresh", dest="identity_threshold", type=float, default=0.20, help="minimum proportion identity to consensus")
	parser.add_argument("--coverage-thresh", dest="coverage_threshold", type=float, default=0.50, help="minimum proportion of coverage of consensus")
	#parser.add_argument("-f", "--fasta", dest="other_fasta_fnames", action='append', default=[], help="other FASTA files")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)
	if not options.fasta_out_fname is None:
		outf = file(options.fasta_out_fname,'w')
		fasta_outs.addStream(outf)
	else:
		# By default, write to stdout
		fasta_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)

	pref_ids = {}
	if not os.path.isfile(options.in_filter_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_filter_fname))
	with open(options.in_filter_fname,'r') as inf:
		tab = util.readTable(inf, header=True)
		pref_ids = dict(zip(tab['species'],tab['orthodb.name']))

	# Now go through headers, find multiples, and select one from each.
	selected_indices = [] # index into headers and sequences
	new_headers = []
	new_seqs = []

	orthodb_ids = [h.split()[-1] for h in headers]
	species_names = [h.split()[0].split('_')[0] for h in headers]
	for species_name in list(set(species_names)):
		dupe_indices = [xi for (xi,spec) in enumerate(species_names) if spec==species_name]
		if len(dupe_indices)==1:
			# No problem, no duplicate
			selected_indices.append(dupe_indices[0])
			continue
		# There are duplicates. Let's see if we already have a preferred strain.
		dupe_odbid = dict([(orthodb_ids[xi],xi) for xi in dupe_indices])
		resolved = False
		try:
			# Look up preferred strain.
			odbid = pref_ids[species_name]
			ind = dupe_odbid[odbid]
			selected_indices.append(ind)
			#print "Resolved", species_name
			resolved = True
		except KeyError:
			pass
		if not resolved:
			# Duplicates, but no preferred strain. Must do something more quantitative.
			selected_indices += dupe_indices

	(headers, seqs) = zip(*[(headers[i],seqs[i]) for i in selected_indices])

	thresh_consensus = translate.ConsensusSequence(options.gap_threshold, '.')
	thresh_consensus.computeFrom(seqs)
	data_outs.write("# Consensus sequence: {:s}\n".format(thresh_consensus))

	# Identify sites with >X identity
	# Prune out duplicates which have the lowest identity to the consensus (or lowest median identity) at those sites
	selected_indices = []
	orthodb_ids = [h.split()[-1] for h in headers]
	species_names = [h.split()[0].split('_')[0] for h in headers]
	for species_name in list(set(species_names)):
		dupe_indices = [xi for (xi,spec) in enumerate(species_names) if spec==species_name]
		indices = sorted([(xi,thresh_consensus.identityTo(seqs[xi])) for xi in dupe_indices], key=lambda x:x[1], reverse=True)
		# Keep index with highest identity to consensus
		# if it's got high enough identity
		if indices[0][1] >= options.identity_threshold:
			coverage = thresh_consensus.coverage(seqs[indices[0][0]])
			if coverage >= options.coverage_threshold:
				selected_indices.append(indices[0][0])
			else:
				data_outs.write("# {} does not have sufficient coverage {:1.3f}\n".format(species_name, coverage))
		else:
			data_outs.write("# {} not identical enough to consensus {:1.3f}\n".format(species_name, indices[0][1]))

	(headers, seqs) = zip(*[(headers[i],seqs[i]) for i in selected_indices])

	# Streamline the headers: add leading ID without strain information
	new_headers = []
	for h in headers:
		species_name = h.split()[0]
		short_species_name = h.split()[0].split('_')[0]
		if species_name != short_species_name:
			new_header = '{} {}'.format(short_species_name, h)
		else:
			new_header = h
		new_headers.append(new_header)
	headers = new_headers

	# Degap the remaining alignment
	# DAD: implement

	# Write output
	n_written = 0
	for (hdr, seq) in zip(headers,seqs):
		line = ">{hdr}\n{seq}\n".format(hdr=hdr, seq=seq)
		fasta_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

