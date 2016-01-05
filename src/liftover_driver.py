#! python

import sys, os, math, random, argparse
import util, biofile, translate, muscle
import liftover

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Lift over domain annotations by alignment")
	# Required arguments
	parser.add_argument(dest="in_alignment_fname", help="input filename")
	parser.add_argument(dest="in_annotation_fname", help="input filename")
	parser.add_argument(dest="alignment_key", help="identifier for source sequence in the alignment")
	parser.add_argument(dest="annotation_key", help="identifier for source sequence in the annotations")
	# Optional arguments
	parser.add_argument("--noalign", dest="dont_align_sequences", action='store_false', default=True, help="don't align sequences (e.g. because previously aligned)?")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
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
	for (k,v) in sorted(optdict.items()):
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

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
	(headers, seqs) = biofile.readFASTA(file(options.in_fname, 'r')) #, key_fxn=biofile.secondField)
	if options.translate_sequences:
		seqs = [translate.translate(s) for s in seqs]
	if not options.dont_align_sequences:
		aligned_seqs = muscle.alignSequences(seqs)
		seqs = aligned_seqs
	zhs = [(h,s) for (h,s) in zip(headers,seqs) if not s is None]
	all_keys = [biofile.firstField(h) for (h,s) in zhs]
	(headers, seqs) = zip(*zhs)
	prot_dict = dict([(biofile.firstField(h), s) for (h,s) in zhs])
	gene_orf_dict = dict([(secondField(h), biofile.firstField(h)) for h in headers])
	orf_gene_dict = dict([(v,k) for (k,v) in gene_orf_dict.items()])
	
	# Write output
	n_written = 0
	data_outs.write("header\n")
	for orf in query_keys:
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

