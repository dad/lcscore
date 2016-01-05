#! python

import sys, os, math, random, argparse, re
import util, gelscore, biofile, translate, protprop, stats

# Objective: Quantify usage and ordering of amino acids.
#
# Example: given sequence QQQQFQRNQPPPY, quantify FY and R usage
#       a = FY
#       b = R
# So there are:
#     2 a
#     1 b
# And motif is:
#     aba
#
# We are not dealing with motifs here, just single-AA usage, so there can be
# no overlap between groups of AAs.

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="slide window")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional arguments
	parser.add_argument("--sentinel-species", dest="sentinel_species", type=str, default=None, help="species identifier for sequence")
	parser.add_argument("--start-sequence", dest="start_sequence", type=str, default=None, help="sequence that starts the domain")
	parser.add_argument("--end-sequence", dest="end_sequence", type=str, default=None, help="sequence that ends the domain")
	parser.add_argument("--start-position", dest="start_position", type=int, default=None, help="starting sequence position for search (1-based)")
	parser.add_argument("--end-position", dest="end_position", type=int, default=None, help="ending sequence position for search (1-based, inclusive)")
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("-a", "--aas", dest="aa_groups", action='append', default=[], help="amino-acid groups to quantify as units, e.g. FYW, DE, KRH, CLIMFV")
	parser.add_argument("-s", "--singles", dest="single_aas", type=str, default='', help="amino acids to count singly")
	parser.add_argument("--nearest", dest="nearest_distances", action='store_true', default=False, help="compute only nearest-neighbor distances?")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--hist-out", dest="histogram_out_fname", default=None, help="histogram output filename")

	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	motif_outs = None

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

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	(headers, seqs) = biofile.readFASTA(file(options.in_fname, 'r')) #, key_fxn=biofile.secondField)
	if options.translate_sequences:
		seqs = [translate.translate(s) for s in seqs]
	zhs = [(h,s) for (h,s) in zip(headers,seqs) if not s is None]
	all_keys = [biofile.firstField(h) for (h,s) in zhs]
	(headers, seqs) = zip(*zhs)
	prot_dict = dict([(biofile.firstField(h), s) for (h,s) in zhs])
	gene_orf_dict = dict([(biofile.secondOrFirstField(h), biofile.firstField(h)) for h in headers])
	orf_gene_dict = dict([(v,k) for (k,v) in gene_orf_dict.items()])

	# Select which genes to process
	query_keys = []
	if not options.query_orf is []:
		# Specific ORF(s)
		query_keys += options.query_orf
	if not options.query_gene is []:
		# Specific gene(s)
		query_keys += [gene_orf_dict[k] for k in options.query_gene]
	if len(query_keys) == 0:
		# Go through all proteins in database
		query_keys = all_keys
	start_position = 0
	end_position = len(seqs[0])

	def found_sentinel(x):
		return options.sentinel_species is None or (options.sentinel_species in x)
	
	def gapfind(seq, substring, start=True, gap='-'):
		gap_pattern = ('[{}]*'.format(gap)).join([a for a in substring])
		pat = re.compile(gap_pattern)
		res = re.search(pat, seq)
		ind = -1
		if not res is None:
			if start:
				ind = res.start()
			else:
				ind = res.end()
		return ind

	if not options.start_sequence is None:
		# Find this sequence
		target_seq = options.start_sequence
		for (h,s) in zip(headers,seqs):
			if found_sentinel(h):
				ind = gapfind(s, target_seq, start=True)
				if ind>0:
					start_position = ind
					# find ending
					if not options.end_sequence is None:
						endpos = gapfind(s, options.end_sequence, start=False)
						assert endpos > 0
						end_position = endpos
					#print s[start_position:end_position].replace('-','')
					# Found the target sequence; bail.
					break

	if not options.start_position is None:
		start_position = options.start_position
	if not options.end_position is None:
		end_position = options.end_position


	# Select out the subset of sequences
	seqslice = slice(start_position-1,end_position)
	for k in query_keys:
		prot_dict[k] = prot_dict[k][seqslice]

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:1]

	# Groups/classes of amino acids to quantify (e.g. A, H, FY, NQ, STY)
	aa_groups = [''.join(sorted(group)) for group in options.aa_groups] + [x for x in options.single_aas]

	comp = protprop.Composition()
	pprop = protprop.ProteinProperties()
	# Compute two distributions:
	# Distance between nearest elements
	# Distance between all elements
	all_distances = {}
	nearest_distances = {}
	distances = {}
	max_length = 0
	for id in query_keys:
		seq = prot_dict[id]
		if len(seq) > max_length:
			max_length = len(seq)
		if options.nearest_distances:
			distances[id] = pprop.nearestDistances(seq, aa_groups)
		else:
			distances[id] = pprop.allDistances(seq, aa_groups)

	# Convert into histograms
	fac = stats.HistogramFactory()
	dist_hists = {}
	for id in query_keys:
		dist_hists[id] = {}
		for aag in aa_groups:
			hist = fac.integerHistogram(0,max_length)
			hist.add(distances[id][aag])
			dist_hists[id][aag] = hist
	
	# Write output
	n_written = 0
	dout = util.DelimitedOutput()
	dout.addHeader('distance','amino-acid distance','d')
	nearest_str = ""
	if options.nearest_distances:
		nearest_str = "nearest-neighbor "
	for id in query_keys:
		for group in aa_groups:
			dout.addHeader('{id:s}.{gp:s}.count'.format(gp=group, id=id),'Distribution of {nearest:s}sequence distances for {gp:s} residues in {id:s}'.format(nearest=nearest_str, gp=group, id=id),'d')
	dout.describeHeader(data_outs)

	# Write out results
	dout.writeHeader(data_outs)
	for dist in range(0,max_length):
		line = "{dist:d}".format(dist=dist)
		for id in query_keys:
			for aag in aa_groups:
				h = dist_hists[id][aag]
				line += "\t{hc:d}".format(hc=h[dist].count)
		line += '\n'
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

