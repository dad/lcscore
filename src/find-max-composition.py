#! python

import sys, os, math, random, argparse
import util, slide, gelscore, biofile

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="slide window")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	parser.add_argument(dest="window_size", type=int, default=None, help="length of window to search for")
	# Optional arguments
	parser.add_argument("--start", dest="start_position", type=int, default=1, help="starting sequence position for search (1-based)")
	parser.add_argument("--end", dest="end_position", type=int, default=None, help="ending sequence position for search (1-based, inclusive)")
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
		query_keys = query_keys[0:min(len(query_keys,100))]
	
	# Write output
	n_written = 0
	dout = util.DelimitedOutput()
	dout.addHeader('orf','S. cerevisiae systematic name','s')
	dout.addHeader('max.position','1-based sequence position of window (start of window) having the maximum score (1 - chi-squared histogram distance on normalized aa-composition histograms)')
	dout.addHeader('total.score','Total score for window sequence across whole sequence')
	dout.addHeader('average.score','Average score (total score/length of whole sequence)')
	dout.addHeader('seq','Sequence of window having maximum score')
	dout.addHeader('sorted.seq','Sorted sequence of window having maximum score')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	for orf in query_keys:
		seq = gelscore.Sequence(prot_dict[orf])
		sw = slide.SequenceWindow(seq, options.window_size)
		# Slide window and collect results
		res = sw.search()
		# Write out results
		line = "{orf:s}\t{r.position:d}\t{r.total_score:1.4f}\t{r.average_score:1.4f}\t{r.sequence:s}\t{srted:s}\n".format(
			orf=orf, r=res, srted=''.join(sorted(res.sequence)))
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

