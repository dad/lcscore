#! python

import sys, os, math, random, argparse
import util, biofile, translate
import gelscore

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Segment sequences by self-similarity")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input filename")
	# Optional arguments
	parser.add_argument("-r", "--max-regions", dest="max_regions", type=int, default=None, help="maximum number of regions per protein")
	parser.add_argument("-t", "--threshold", dest="score_threshold", type=float, default=0.9, help="minimum score to include a residue in a region")
	parser.add_argument("-z", "--min-size", dest="min_region_size", type=int, default=20, help="minimum size of a region in residues")
	#parser.add_argument("-m", "--merge-similarity", dest="merge_similarity", type=float, default=0.99, help="similarity threshold (-1,1) for merging identified regions")
	parser.add_argument("-b", "--base", dest="logarithm_base", type=int, default=20, help="base of logarithm for entropy calculations")
	parser.add_argument("-d", "--max-distance", dest="max_distance", type=int, default=None, help="maximum distance over which to score similarity")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
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
	
	
	# Distance definition
	dist = gelscore.SimilarityWeight(sim_add=1.0, max_dist=options.max_distance)
	# Region finder
	reg_finder = gelscore.NormalizedCutRegionFinder(dist, min_region_size=options.min_region_size, 
		score_threshold=options.score_threshold, max_regions=options.max_regions)

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:min(len(query_keys,100))]

	region_dict = {}
	n_processed = 0
	for orf in query_keys:
		seq = gelscore.Sequence(prot_dict[orf])
		#print str(seq)
		# Find the regions
		res = reg_finder.find(seq)
		info_outs.write("# Found {:d} regions for {:s}\n".format(len(res), orf))
		region_dict[orf] = res
		n_processed += 1
		if n_processed % 200 == 0:
			info_outs.write("{:d} ".format(n_processed))
	info_outs.write("\n")
	
	# Write output
	quant = gelscore.EntropyQuant()
	n_written = 0
	data_outs.write("id\torf\tgene\tnum.regions\tstart\tend\tlength\tsim.score\tentropy\tavg.aas\tsequence\n")
	for orf in query_keys:
		res = region_dict[orf]
		seq = prot_dict[orf]
		nreg = 0
		res.sort(key=lambda x: x.region.start)
		for r in res:
			nreg += 1
			id = "{}.{}".format(orf, str(nreg).zfill(2))
			entropy = quant.quant(r.region.sequence, base=options.logarithm_base)
			line = "{id}\t{orf}\t{gene}\t{nreg:d}\t{start:d}\t{end:d}\t{len:d}\t{score:1.3f}\t{ent:1.3f}\t{chars:1.3f}\t{seq:s}\n".format(
				id=id, orf=orf, gene=orf_gene_dict[orf], nreg=len(res), start=r.region.start+1, end=r.region.end, len=len(r.region), 
				score=r.score, ent=entropy, chars=options.logarithm_base**entropy, seq=r.region.sequence)
			data_outs.write(line)
			n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

