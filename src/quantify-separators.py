#! python

import sys, os, math, random, argparse
import util, gelscore, biofile, translate, protprop, stats

# Objective: Quantify distribution of usage of amino acid X between instances of amino acid Y.
#   e.g. with sequence XXyyXyXyyyX
#   output is 0, 2, 1, 3
# and accounting for other residues
#   e.g. XaXyayXyaXyyyaaX
#   output is 1,3,2,5,  0,2,1,3

def maskSequence(s, posts, pickets, postrep='p', picketrep='k', otherrep='x'):
	res = ''
	for aa in s:
		if aa in posts:
			res += postrep
		elif aa in pickets:
			res += picketrep
		else:
			res += otherrep
	return ''.join(res)

def picketDistribution(s, posts, pickets, include_others=False):
	seq = s
	target_aas = posts + pickets
	if not include_others:
		# eliminate everything that's not a post or picket
		seq = [aa for aa in seq if aa in target_aas]
	post_inds = [xi for (xi,aa) in enumerate(seq) if aa in posts]
	# Get distribution
	res = []
	for xi in range(len(post_inds)-1):
		res.append(post_inds[xi+1]-post_inds[xi]-1)
	return res


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Compute distribution of amino acids between 'fencepost' residues")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input filename")
	parser.add_argument(dest="posts", type=str, help="boundary residues")
	parser.add_argument(dest="pickets", type=str, help="intervening residues")
	# Optional arguments
	parser.add_argument("--start", dest="start_position", type=int, default=1, help="starting sequence position for search (1-based)")
	parser.add_argument("--end", dest="end_position", type=int, default=None, help="ending sequence position for search (1-based, inclusive)")
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("--max-dist", dest="max_dist", type=int, default=None, help="maximum distance to track")
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

	# Select out the subset of sequences
	seqslice = slice(options.start_position-1,options.end_position)
	for k in query_keys:
		prot_dict[k] = prot_dict[k][seqslice]

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:1]

	# Compute distribution for all separators
	# Compute distribution for picket residues
	# Track maximum distance
	dists_dict = {}
	max_dist = 0
	post_counts = {}
	picket_counts = {}
	for id in query_keys:
		seq = prot_dict[id]
		all_dist = picketDistribution(seq, options.posts, options.pickets, include_others=True)
		#print id, sorted(all_dist)
		max_dist = max(max_dist, max(all_dist)) # Must be >= max picket dist
		picket_dist = picketDistribution(seq, options.posts, options.pickets, include_others=False)
		#print id, sorted(picket_dist)
		dists_dict[id] = (all_dist, picket_dist)
		# Counts
		post_counts[id] = len([x for x in seq if x in options.posts])
		picket_counts[id] = len([x for x in seq if x in options.pickets])

	if not options.max_dist is None:
		max_dist = min(options.max_dist, max_dist)

	# Convert into histograms
	fac = stats.HistogramFactory()
	all_dist_hists = {}
	picket_dist_hists = {}
	for id in query_keys:
		(all_dist, picket_dist) = dists_dict[id]
		hist = fac.integerHistogram(0,max_dist)
		hist.add(all_dist)
		all_dist_hists[id] = hist
		hist = fac.integerHistogram(0,max_dist)
		hist.add(picket_dist)
		picket_dist_hists[id] = hist

	dist_range = range(max_dist+1)
	
	# Write output
	n_written = 0
	dout = util.DelimitedOutput()
	dout.addHeader('id','protein identifier','s')
	dout.addHeader('{}.count'.format(options.posts), '# instances of {} in region'.format(options.posts), 'd')
	dout.addHeader('{}.count'.format(options.pickets), '# instances of {} in region'.format(options.pickets), 'd')
	for prefix in ['all',options.pickets]:
		for xi in dist_range:
			dout.addHeader('{}.{:d}'.format(prefix, xi),'\t# of times there are {:d} {} residues between {} residues'.format(xi, prefix, options.posts), 'd')
	
	dout.describeHeader(data_outs)

	# Write out results
	dout.writeHeader(data_outs)
	for id in query_keys:
		line = "{id:s}\t{postcount:d}\t{picketcount:d}".format(id=id, postcount=post_counts[id], picketcount=picket_counts[id])
		h = all_dist_hists[id]
		#print h
		for xi in dist_range:
			dist_count = h[xi].count
			line += "\t{:d}".format(dist_count)
		h = picket_dist_hists[id]
		for xi in dist_range:
			dist_count = h[xi].count
			line += "\t{:d}".format(dist_count)
		line += '\n'
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

