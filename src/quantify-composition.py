#! python

import sys, os, math, random, argparse, re
import util, gelscore, biofile, translate, protprop, scriptutil

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
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--isolate-out", dest="isolate_out_fname", default=None, help="output filename")

	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	isolate_outs = util.OutStreams()
	params_outs = util.OutStreams([data_outs])
	motif_outs = None

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		if options.isolate_out_fname is None:
			data_outs.addStream(sys.stdout)
	if not options.isolate_out_fname is None:
		isolate_outf = file(options.isolate_out_fname,'w')
		isolate_outs.addStream(isolate_outf)

	# Write out parameters
	params_outs.write("# Run started {}\n".format(util.timestamp()))
	params_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	params_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		params_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

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

	start_position, end_position = scriptutil.findSequencePositions(options.start_position, options.end_position, 
		options.start_sequence, options.end_sequence, options.sentinel_species, headers, seqs)

	data_outs.write("# Start = {}, end = {}\n".format(start_position, end_position))
	# Select out the subset of sequences
	seqslice = slice(start_position-1,end_position)
	for k in query_keys:
		prot_dict[k] = prot_dict[k][seqslice]

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:min(len(query_keys,100))]

	# Groups/classes of amino acids to quantify (e.g. A, H, FY, NQ, STY)
	aa_groups = [''.join(sorted(group)) for group in options.aa_groups] + [x for x in options.single_aas]
	all_aas = ''.join(aa_groups)
	# Write output
	comp = protprop.Composition()
	pprop = protprop.ProteinProperties()
	n_written = 0
	n_isolate_written = 0
	dout = util.DelimitedOutput()
	dout.addHeader('id','protein identifier','s')
	for group in aa_groups:
		sorted_group = ''.join(sorted(group))
		dout.addHeader('{gp:s}.count'.format(gp=sorted_group),'Number of instances of {gp:s}'.format(gp=sorted_group),'d')
	dout.addHeader('total.length','Total amino acids','d')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	for id in query_keys:
		seq = prot_dict[id]

		counts = pprop.counts(seq, aa_groups)
		count_str = '\t'.join(['{:d}'.format(c) for c in counts])
		# Write out results
		line = "{id:s}\t{counts:s}\t{totlength:d}\n".format(id=id, counts=count_str, totlength=len(seq))
		data_outs.write(line)
		n_written += 1
		squeezed = ''.join([s for s in seq if s in all_aas])
		line = ">{id:s}\n{squeezed:s}\n".format(id=id, squeezed=squeezed)
		isolate_outs.write(line)
		n_isolate_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

	if not options.isolate_out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_isolate_written, options.isolate_out_fname))
		isolate_outf.close()

