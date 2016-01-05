#! python

import sys, os, math, random, argparse
import util, biofile
import structevol

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
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

	# Quant the alignment
	pal = structevol.ProfileAlignment(seqs)

	
	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('position','Column in the alignment (1-based)','d')
	dout.addHeader('hamming.diff','Hamming distance between this and previous column','d')
	dout.addHeader('hamming.prop','Hamming distance as a proportion of maximum distance between this and previous column','1.4f')
	dout.addHeader('entropy','Sequence entropy for this column','1.4f')
	dout.addHeader('pr.ungapped','Proportion ungapped in this column','1.4f')
	dout.addHeader('consensus','Consensus amino acid for this column','s')
	dout.addHeader('consensus.nogaps','Consensus amino acid for this column, excluding gaps','s')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	format = dout.getFormat(named=True)
	n_written = 0

	for (pos, ham) in enumerate(pal.distances()): #range(len(seqs[0])):
		col = pal.column(pos)
		ent = structevol.entropy(col)
		N = float(len(col))
		cov = (N-col.count('-'))/N
		line = format.format(position=pos, hamming_diff=ham, hamming_prop=ham/N, entropy=ent, pr_ungapped=cov, consensus=structevol.mostFrequent(col), consensus_nogaps=structevol.mostFrequent(col, nogap=True))
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

