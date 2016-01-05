#! python

import sys, os, math, random, argparse
import util, biofile, stats
import gelscore

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Quantify distribution of sequence entropy of random sequences with a fixed window size and print histogram")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input filename")
	parser.add_argument(dest="window_size", type=int, help="size of window, in residues")
	# Optional arguments
	parser.add_argument("-n", "--num-samples", dest="num_samples", type=int, default=1e4, help="number of samples to draw")
	parser.add_argument("-m", "--num-bins", dest="num_bins", type=int, default=1000, help="number of bins for histogram")
	parser.add_argument("-u", "--upper-window-size", dest="upper_window_size", type=int, default=None, help="number of samples to draw")
	parser.add_argument("-b", "--base", dest="logarithm_base", type=int, default=20, help="base of logarithm for entropy calculations")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Check
	assert options.window_size > 0
	if not options.upper_window_size is None:
		assert options.upper_window_size >= options.window_size
	else:
		options.upper_window_size = options.window_size
	assert options.window_size > 0

	prot_dict = {}
	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	prot_dict = biofile.readFASTADict(file(options.in_fname, 'r'))
	
	# Generate sequence windows and quantify them
	seq_weights = [(s,1.0/len(s)) for s in prot_dict.values()]
	
	window_sizes = range(options.window_size, options.upper_window_size+1, 1)
	for window_size in window_sizes:
		# Start up output
		if not options.out_fname is None:
			if len(window_sizes)>1:
				# Use formatted filename for each window size
				fname = "{}-{:d}mers.txt".format(options.out_fname, window_size)
			else:
				# Use filename as given, for single file
				fname = options.out_fname
			outf = file(fname,'w')
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
		data_outs.write("# Window size = {}\n".format(window_size))
	
		n_samples = 0
		scores = []
		min_ent = 1e6
		min_seq = None
		max_ent = 0
		max_seq = None
		quant = gelscore.EntropyQuant()
		while n_samples < options.num_samples:
			# Pick sequence at random
			seq = stats.weighted_choice(seq_weights)
			if len(seq) >= window_size:
				# Pick a subsequence within sequence at random
				startpos = random.randint(0, len(seq)-window_size)
				subseq = seq[startpos:(startpos+window_size)]
				if not '*' in subseq:
					score = quant.quant(subseq, base=options.logarithm_base)
					if score < min_ent:
						min_ent = score
						min_seq = subseq
					if score > max_ent:
						max_ent = score
						max_seq = subseq
					scores.append(score)
					n_samples += 1

		data_outs.write("# Minimum-entropy sequence = {}\n".format(min_seq))
		data_outs.write("# Maximum-entropy sequence = {}\n".format(max_seq))

		hist = stats.Histogram([20**s for s in scores], options.num_bins)
		data_outs.write("ent.lower\tent.mid\tent.upper\tcount\tdensity\tcum.count\tcum.density\n")
		n_written = 0
		for b in hist.bins:
			line = "{me.lower:1.4f}\t{me.mid:1.4f}\t{me.upper:1.4f}\t{me.count:d}\t{me.density:1.3E}\t{me.cum:d}\t{me.cumdensity:1.3E}\n".format(me=b)
			data_outs.write(line)
			n_written += 1

		# Write out stopping time
		data_outs.write("# Run finished {}\n".format(util.timestamp()))

		# Shut down output
		if not options.out_fname is None:
			info_outs.write("# Wrote {} lines to {}\n".format(n_written, fname))
			data_outs.removeStream(outf)
			outf.close()
		else:
			data_outs.removeStream(sys.stdout)

