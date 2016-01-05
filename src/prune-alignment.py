#! python

import sys, os, math, random, argparse
import util, biofile, muscle, translate

'''
Prune out:
1) 
Choose 
'''

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prune alignment to a single sequence per species")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")
	#parser.add_argument(dest="in_filter_fname", help="input species-to-keep filename")
	# Optional arguments
	parser.add_argument("--anchor", dest="anchor", action='append', help="identifiers to serve as the anchor(s) for the alignment")
	parser.add_argument("--gap-thresh", dest="gap_threshold", type=float, default=0.20, help="maximum proportion gaps to contribute to consensus/anchor")
	parser.add_argument("--identity", dest="identity_threshold", type=float, default=0.0, help="minimum proportion identity to consensus/anchor")
	parser.add_argument("--coverage", dest="coverage_threshold", type=float, default=0.50, help="minimum proportion of coverage of consensus/anchor")
	parser.add_argument("-1", "--one-per-species", dest="one_per_species", action='store_true', default=False, help="choose only one member of each species?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove all-species gaps from resulting alignment?")
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
	info_outs.write("# Read {:d} sequences\n".format(len(seqs)))

	# Either get the anchor sequence, or construct the consensus as the anchor
	anchor = translate.ConsensusSequence(options.gap_threshold, '.')
	if not options.anchor is None:
		anchor_seqs = []
		for (xi, h) in enumerate(headers):
			for a in options.anchor:
				if a in h:
					anchor_seqs.append(seqs[xi])
		if len(anchor_seqs) == 0:
			info_outs.write("# Could not find sequences '{}'; exiting\n".format(','.join(options.anchor)))
			sys.exit()
		anchor.computeFrom(anchor_seqs)
	else:
		# Build consensus.
		# DAD: implement
		anchor.computeFrom(seqs)
		data_outs.write("# Consensus sequence: {:s}\n".format(anchor))

	#print anchor

	# Identify sites with >X identity
	# Prune out duplicates which have the lowest identity to the anchor (or lowest median identity) at those sites
	selected_indices = []
	orthodb_ids = [h.split()[-1] for h in headers]
	species_names = [h.split()[0].split('_')[0].split("/")[0] for h in headers]
	for species_name in sorted(list(set(species_names))):
		dupe_indices = [xi for (xi,spec) in enumerate(species_names) if spec==species_name]
		indices = sorted([(xi,anchor.identityTo(seqs[xi])) for xi in dupe_indices], key=lambda x:x[1], reverse=True)
		# Keep index with highest identity to anchor
		# if it's got high enough identity
		#print species_name, indices[0][1]
		tups = indices
		if options.one_per_species:
			# Retain only the highest-identity ortholog
			tups = [indices[0]]
		for tup in tups:
			if tup[1] >= options.identity_threshold:
				coverage = anchor.coverage(seqs[tup[0]])
				if coverage >= options.coverage_threshold:
					selected_indices.append(tup[0])
				else:
					data_outs.write("# {} does not have sufficient coverage {:1.3f}\n".format(species_name, coverage))
			else:
				data_outs.write("# {} not identical enough to anchor/consensus {:1.3f}\n".format(species_name, tup[1]))

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
	if options.degap:
		gap = '-'
		new_seqs = []
		n = len(seqs)
		# Find positions with at least one non-gap character
		nongap = [xi for xi in range(len(seqs[0])) if (''.join([s[xi] for s in seqs])).count(gap)<n]
		#print nongap
		for s in seqs:
			new_s = ''.join([s[xi] for xi in nongap])
			new_seqs.append(new_s)
		seqs = new_seqs

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

