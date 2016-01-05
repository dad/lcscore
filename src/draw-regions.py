#! python

import sys, os, math, random, argparse, bisect
import util, na, stats, gelscore, biofile
import drawdomain

'''
Goal: Factor self-similar sequences by composition and score them.

Procedure:
	- Collect regions
	- Apply cutoff for self-similarity
	- Take composition vectors and compute correlation matrix
	- Factor correlation matrix by eigendecomposition
	- Write out ED
	- For all sequences, write out factorization of composition along with start, end, etc.

'''

class DomainInfo(object):
	def __init__(self, bstart, bend, color):
		self.bstart = bstart
		self.bend = bend
		self.color = color

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Draw protein regions")
	# Required arguments
	parser.add_argument(dest="in_fasta_fname", help="input FASTA filename")
	parser.add_argument(dest="in_domain_fname", help="input domain filename")
	# Optional arguments
	parser.add_argument("-r", "--max-regions", dest="max_regions", type=int, default=None, help="maximum number of regions per protein")
	parser.add_argument("-t", "--threshold", dest="score_threshold", type=float, default=0.01, help="minimum P-value to draw a region")
	parser.add_argument("-z", "--min-size", dest="min_region_size", type=int, default=20, help="minimum size of a region in residues")
	parser.add_argument("--row-height", dest="row_height", type=float, default=1.2, help="height of a single protein diagram in cm")
	parser.add_argument("-x", "--width", dest="width", type=float, default=20.0, help="width of a protein diagram in cm")
	parser.add_argument("-y", "--height", dest="height", type=float, default=None, help="height of a protein diagram in cm")
	parser.add_argument("-d", "--domain-height", dest="domain_height", type=float, default=0.5, help="height of a protein domain in cm")
	parser.add_argument("--label-width", dest="label_width", type=float, default=2.0, help="width of a protein label in cm")
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	data_outs.addStream(sys.stdout)

	# Write out parameters
	info_outs.write("# Run started {}\n".format(util.timestamp()))
	info_outs.write("# Parameters:\n")
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
	prot_dict = {}
	# Read input
	if not os.path.isfile(options.in_fasta_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fasta_fname))
	(headers, seqs) = biofile.readFASTA(file(options.in_fasta_fname, 'r'))
	all_keys = [biofile.firstField(h) for h in headers]
	if options.translate_sequences:
		seqs = [translate.translate(s) for s in seqs]
	zhs = [(h,s) for (h,s) in zip(headers,seqs) if not s is None]
	(headers, seqs) = zip(*zhs)
	prot_dict = dict([(biofile.firstField(h), s) for (h,s) in zhs])
	gene_orf_dict = dict([(secondField(h), biofile.firstField(h)) for h in headers])
	orf_gene_dict = dict([(o,g) for (g,o) in gene_orf_dict.items()])
	
	query_keys = []
	if not options.query_orf is []:
		# Specific ORF(s)
		query_keys += options.query_orf
	if not options.query_gene is []:
		# Specific gene(s)
		query_keys += [gene_orf_dict[k] for k in options.query_gene] # if k.startswith(options.query_gene)])
	if len(query_keys) == 0:
		# Go through all proteins in database
		query_keys = all_keys
	
	# Load domain/color boundaries
	if not os.path.isfile(options.in_domain_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_domain_fname))
	dlr = util.DelimitedLineReader(file(options.in_domain_fname, 'r'), header=True)
	domain_boundaries = util.listdict()
	domain_colors = util.listdict()
	for flds in dlr.dictentries:
		orf = flds['orf']
		score = flds['p.value']
		#if orf in query_keys:
		#	print orf, score, flds['start'], flds['end'], flds['sequence']
		if score < options.score_threshold:
			domain_boundaries[orf].append((flds['start'], flds['end']))
			domain_colors[orf].append(flds['color'].replace("0x","#"))
	
	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	# Find maximum length
	max_length = 0
	for k in query_keys:
		seq = prot_dict[k]
		if len(seq)>max_length:
			max_length = len(seq)
	
	if options.debugging:
		query_keys = query_keys[0:min(len(query_keys,100))]

	xmargin = 0.5
	ymargin = 1.0
	num_rows = len(query_keys)
	if options.height is None:
		options.height = 2*ymargin + options.row_height * num_rows
	
	# Protein starts at 2*xmargin, ends at width - xmargin
	diagram_width = (options.width-options.label_width-3.0*xmargin)
	cm_per_aa = diagram_width/max_length
		
	# New drawing
	draw = drawdomain.SequenceDrawing(options.width, options.height, 
		cm_per_aa=cm_per_aa, startx=xmargin*2 + options.label_width, starty=ymargin, label_startx=options.label_width+xmargin,
		domain_height=options.domain_height, row_height=options.row_height)

	n_written = 0
	for orf in query_keys:
		seq = gelscore.Sequence(prot_dict[orf])
		gene_name = orf_gene_dict[orf]
		#print gene_name
		boundaries = domain_boundaries[orf]
		colors = domain_colors[orf]
		# Draw the sequence
		draw.cartoon(seq, gene_name, boundaries, colors)
		draw.label(drawdomain.proteinName(gene_name))
		draw.nextRow()
		n_written += 1

	if not options.out_fname is None:
		draw.write(options.out_fname)

	# Write out stopping time
	info_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	info_outs.write("# Wrote {} diagrams to {}\n".format(n_written, options.out_fname))

