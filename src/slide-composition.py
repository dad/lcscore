#! python

import sys, os, math, random, argparse
import util, slide, gelscore, biofile

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="slide window")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	#parser.add_argument(dest="motif", default=None, help="sequence motif to score")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("--motif", "-m", dest="motif", default=None, help="sequence motif to score")
	group.add_argument("--composition-file", "-c", dest="composition_fname", default=None, help="filename containing composition information")
	# Optional arguments
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
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
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	def secondField(h):
		f = None
		try:
			f = biofile.secondField(h)
		except:
			f = biofile.firstField(h)
		return f

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	(headers, seqs) = biofile.readFASTA(file(options.in_fname, 'r'))
	if options.translate_sequences:
		seqs = [translate.translate(s) for s in seqs]
	zhs = [(h,s) for (h,s) in zip(headers,seqs) if not s is None]
	all_keys = [biofile.firstField(h) for (h,s) in zhs]
	(headers, seqs) = zip(*zhs)
	prot_dict = dict([(biofile.firstField(h), s) for (h,s) in zhs])
	gene_orf_dict = dict([(secondField(h), biofile.firstField(h)) for h in headers])
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

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:min(len(query_keys,100))]
	
	# Set up motif to compare
	fac = slide.SequenceCompositionComparatorFactory()
	comparator = fac.make(options.motif, 'sequence')
	window_size = len(options.motif)
	
	# Write output
	n_written = 0
	for orf in query_keys:
		seq = gelscore.Sequence(prot_dict[orf])
		sw = slide.SequenceWindow(seq, window_size)
		# Slide window and collect results
		ress = sw.slide(comparator)
		# Write out results
		data_outs.write("pos\taa\tscore\twindow\n")
		for (xi, res) in enumerate(ress.results()):
			line = "{res.position}\t{res.aa}\t{res.score}\t{res.window}\n".format(res=res)
			data_outs.write(line)
		data_outs.write("\n") # empty line between orfs
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

