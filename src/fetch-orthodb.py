#! python

import sys, os, math, random, argparse
import util, biofile, muscle, translate
import urllib


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Fetch OrthoDB alignment")
	# Required arguments
	parser.add_argument(dest="orthodb_id", type=str, help="OrthoDB ID")
	# Optional arguments
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	#parser.add_argument("-o", "--out-dir", dest="out_dir", default='', help="output file directory")
	options = parser.parse_args()

	info_outs = util.OutStreams()
	#data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		info_outs.addStream(outf)
	else:
		# By default, write to stdout
		info_outs.addStream(sys.stdout)
	if not options.fasta_out_fname is None:
		outf = file(options.fasta_out_fname,'w')
		fasta_outs.addStream(outf)
	else:
		# By default, write to stdout
		fasta_outs.addStream(sys.stdout)

	# Write out parameters
	info_outs.write("# Run started {}\n".format(util.timestamp()))
	info_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	info_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		info_outs.write("#\t{k}: {v}\n".format(k=k, v=v))


	local_fname = options.fasta_out_fname
	if local_fname is None:
		local_fname = "tmp{:d}".format(random.randint(0, 1e20))
	# Fetch file from OrthoDB
	if not options.orthodb_id is None:
		#local_fname = "uniprot-yeast.txt"
		remote_fname = "http://cegg.unige.ch/orthodb7/fasta.fasta?ogs={:s}".format(options.orthodb_id)
		urllib.urlretrieve(remote_fname, local_fname)
		print "# Downloaded {} to {}".format(remote_fname, local_fname)
		info_outs.write("# Downloaded {} to {}\n".format(remote_fname, local_fname))

	# Read input
	if not os.path.isfile(local_fname):
		raise IOError("# Error: file {} does not exist".format(local_fname))
	with open(local_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
		if options.fasta_out_fname is None:
			# Write data
			biofile.writeFASTA(seqs, fasta_outs, headers=headers)

	# Write out stopping time
	info_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Fetched {} sequences to {}\n".format(len(headers), options.fasta_out_fname))
		outf.close()

