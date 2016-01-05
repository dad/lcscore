#! python

import sys, os, math, random, argparse, bisect
import util, na, stats, gelscore
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

def selfSimilarity(seq, dist):
	pass

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Factor and score self-similar regions by composition")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input filename")
	# Optional arguments
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
	data_outs.write("# Run started {}\n".format(util.timestamp(timeformat="%a %b %d %H:%M:%S %Y")))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))
	
	# Read input
	data = {}
	headers = None
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	inf = file(options.in_fname, 'r')
	dlr = util.DelimitedLineReader(inf, header=True)
	headers = dlr.headers
	for flds in dlr.dictentries:
		data[flds['id']] = flds
	inf.close()
	
	# Distance measure
	dist = gelscore.SimilarityWeight(sim_add=1.0, max_distance=200)
	
	# Arrange the data
	keys = sorted(data.keys())
	
	for k in keys:
		dat = data[k]
		seq = dat['sequence']
		sim = selfSimilarity(seq, dist)
	
	# Write output
	n_written = 0
	# Tack sequence on the end, for better readability
	columns = headers[:]
	columns.remove('sequence')
	data_outs.write('\t'.join(columns) + "\tp.value\tp.value.adj\tsequence\n")
	for key in keys:
		flds = data[key]
		#entropy = flds['entropy']
		#length = flds['length'] # length of the region
		pval = p_value_dict[key]
		adj_pval = None
		if options.adjust:
			adj_pval = adjusted_p_value_dict[key]
		line = '\t'.join(['{}'.format(flds[col]) for col in columns])
		line += "\t{pval}\t{adj_pval}".format(pval=na.formatNA(pval, format="{:1.3E}"), adj_pval=na.formatNA(adj_pval, format="{:1.3E}"))
		line += '\t' + flds['sequence'] + '\n'
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

