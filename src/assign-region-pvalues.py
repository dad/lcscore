#! python

import sys, os, math, random, argparse, bisect
import util, na, stats
from numpy.lib.scimath import logn

# Goal: read table of called regions, assign P-value for those regions' entropies based on empirical probabilities, optionally adjust P values for multiple testing.
# Write out list of (possibly adjusted) P-values.

def findPValue(x, entropy_list, pval_list):
	# plist consists of (x, probability_that_random_seq_has_value_less_than_x)
	# assume plist is sorted
	#print entropy_list[0:10]
	pos = bisect.bisect_left(entropy_list, x)
	n = len(pval_list)
	if pos == n:
		pos = n-1
	return pval_list[pos]

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Compute P-values for entropy")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input filename")
	parser.add_argument(dest="in_distribution_fname_pattern", help="input distribution filename pattern")
	# Optional arguments
	parser.add_argument("-a", "--adjust", dest="adjust", action='store_true', default=False, help="adjust P-values")
	parser.add_argument("-l", "--lower", dest="min_window_size", default=5, help="minimum window size")
	parser.add_argument("-u", "--upper", dest="max_window_size", default=300, help="maximum window size")
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
	
	# Read densities
	# Density files should be specified with a regexp pattern matching only the length of the window.
	assert '*' in options.in_distribution_fname_pattern
	entropies = {}
	densities = {}
	for winsize in range(options.min_window_size, options.max_window_size,1):
		fname = options.in_distribution_fname_pattern.replace('*', '{:d}'.format(winsize))
		#if not os.path.isfile(fname):
		#	raise IOError("# Error: file {} does not exist".format(fname))
		try:
			inf = file(fname, 'r')
			dlr = util.DelimitedLineReader(inf, header=True)
			dens = []
			ent = []
			for flds in dlr.dictentries:
				bin_upper = flds['ent.upper']
				p_entropy_lower = flds['cum.density']
				ent.append(bin_upper)
				dens.append(p_entropy_lower)
			densities[winsize] = dens
			entropies[winsize] = ent
		except IOError:
			info_outs.write("# File {} not found\n".format(fname))
			continue
		except util.ReaderEOFError:
			info_outs.write("# EOF encountered in file {}\n".format(fname))
			continue

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
	
	# Arrange the data
	keys = sorted(data.keys())
	
	# Adjust P values? If so, need to get all P values together
	p_value_dict = {}
	length_dict = {}
	adjusted_p_value_dict = {}
	p_value_list = []
	for key in keys:
		flds = data[key]
		entropy = flds['entropy']
		length = flds['length'] # length of the region
		if length > options.max_window_size:
			length = options.max_window_size
		elif length < options.min_window_size:
			length = options.min_window_size
		try:
			d = densities[length]
			e = entropies[length]

			pval = findPValue(20**entropy, e, d)
			p_value_dict[key] = pval
			p_value_list.append(pval)
			try:
				length_dict[length].append(key)
			except KeyError:
				length_dict[length] = [key]
		except KeyError:
			print "Erroring...this should not continue"
			p_value_dict[key] = 1.0
			p_value_list.append(1.0)
			continue
	if options.adjust:
		# Each distribution, based on length, gets its own correction
		for length in length_dict.keys():
			len_keys = length_dict[length]
			p_values = [p_value_dict[lk] for lk in len_keys]
			adjusted_p_values = stats.adjustPValues(p_values, method='FDR')
			for (i,k) in enumerate(len_keys):
				adjusted_p_value_dict[k] = adjusted_p_values[i]
		#adjusted_p_value_dict = dict([(k, v) for (k,v) in zip(keys, adjusted_p_values)])
	
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

