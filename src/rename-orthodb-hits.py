#! python

import sys, os, math, random, argparse
import util, biofile, muscle, translate

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Process and name OrthoDB downloads")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	parser.add_argument(dest="in_names_fname", default=None, help="input filename")
	# Optional arguments
	parser.add_argument("-f", "--fasta", dest="other_fasta_fnames", action='append', default=[], help="other FASTA files")
	parser.add_argument("-x", "--remove-x", dest="remove_x", action='store_true', default=False, help="remove sequences containing X's?")
	parser.add_argument("-a", "--align", dest="align", action='store_true', default=False, help="align sequences?")
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
	for fname in options.other_fasta_fnames:
		if not os.path.isfile(fname):
			raise IOError("# Error: file {} does not exist".format(fname))
		with open(fname,'r') as inf:
			# Read a FASTA file?
			(new_headers, new_seqs) = biofile.readFASTA(inf)
			headers = headers + new_headers
			seqs = seqs + new_seqs

	if not os.path.isfile(options.in_names_fname):
	 	raise IOError("# Error: file {} does not exist".format(options.in_names_fname))
	with open(options.in_names_fname,'r') as inf:
		species = util.readTable(inf, header=True)

	def shorten(x):
		flds = x.split(' ')
		res = '{g}.{s}'.format(g=flds[0][0], s='_'.join(flds[1:]))
		return res

	species_name_lookup = dict(zip(species['ODB_code'], [shorten(x) for x in species['Organism']]))
	#print species_name_lookup

	# Preserve all the sequences
	# Align them

	# Store output
	named_seqs = {}
	named_headers = {}
	new_headers = []
	new_seqs = []
	for (hdr, seq) in zip(headers,seqs):
		odb_code = hdr.strip().split()[-1]
		assert len(odb_code)==5
		species_name = species_name_lookup[odb_code]
		store_seq = True
		if options.remove_x:
			store_seq = not 'X' in seq
		if '_' in seq:
			store_seq = False
		if store_seq:
			new_header = '{sname} {hdr}'.format(sname=species_name, hdr=hdr)
			#named_headers[species_name] = new_header
			#named_seqs[species_name] = seq
			new_headers.append(new_header)
			new_seqs.append(seq)

	#(headers, seqs) = zip(*[(named_headers[x],named_seqs[x]) for x in sorted(named_seqs.keys())])
	headers = new_headers
	seqs = new_seqs

	if options.align:
		#print os.path.expanduser('/cygdrive/f/develop/muscle3.8.31/muscle')
		aligned_seqs = muscle.alignSequences(seqs) #, exepath='~\\develop\\muscle3.8.31\\muscle')
		seqs = aligned_seqs



	# Write output
	n_written = 0
	for (hdr, seq) in zip(headers,seqs):
		line = ">{hdr}\n{seq}\n".format(hdr=hdr, seq=seq)
		data_outs.write(line)
		n_written += 1


	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

