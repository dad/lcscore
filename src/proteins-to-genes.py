#! python

import sys, os, math, random, argparse
import util, biofile, translate

def makeMutantFromSequence(target_protein_seq, base_dna_seq):
	codons = [x for x in translate.codons(base_dna_seq)]
	base_prot_seq = translate.translate(base_dna_seq)
	assert len(base_prot_seq) == len(target_protein_seq)
	mutant_dna_seq = ''
	for (i, aa) in enumerate(target_protein_seq):
		if aa == base_prot_seq[i]:
			mutant_dna_seq += codons[i]
		else:
			mutant_dna_seq += translate.randomReverseTranslate(aa)
	assert translate.translate(mutant_dna_seq) == target_protein_seq
	return mutant_dna_seq

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input FASTA filename")
	# Optional arguments
	parser.add_argument("--suggest", dest="suggest_sequences", default=None, help="gene sequence to prepend")
	parser.add_argument("-s", "--seed", dest="random_seed", type=float, default=111, help="seed for random number generator")
	parser.add_argument("--repeat", dest="repeat", type=int, default=1, help="number of randomizations to perform")
	parser.add_argument("-r", "--randomize", dest="randomize", action="store_true", default=False, help="randomize input aa sequence?")
	parser.add_argument("--prefix", dest="prefix", default='', help="gene sequence to prepend")
	parser.add_argument("--suffix", dest="suffix", default='', help="gene sequence to append")
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

	# Set seed
	random.seed(options.random_seed)

	# Read input
	if not os.path.isfile(options.in_fname):
	 	raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
	 	# Read a FASTA file?
	 	(headers, seqs) = biofile.readFASTA(inf)
	
	sug_dict = {}
	if not options.suggest_sequences is None:
		if not os.path.isfile(options.suggest_sequences):
		 	raise IOError("# Error: file {} does not exist".format(options.in_fname))
		with open(options.suggest_sequences,'r') as inf:
		 	# Read a FASTA file?
		 	sug_dict = biofile.readFASTADict(inf)
	
	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('name','name of construct')
	dout.addHeader('sequence','sequence')
	dout.addHeader('notes','notes')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	n_written = 0
	mutant_seqs = {}

	def parseHeader(x):
		name = biofile.firstField(x)
		property_entries = [tuple(y.split('=')) for y in x.split() if '=' in y]
		props = dict(property_entries)
		return name, props

	def randomizeSequence(x):
		return ''.join(random.sample(x,len(x)))

	if options.randomize and options.repeat>1:
		new_headers = []
		new_seqs = []
		for xi in range(len(headers)):
			for rep in range(options.repeat):
				flds = headers[xi].split()
				new_name = '{}_rand{}'.format(flds[0],rep+1)
				new_header = '{}'.format(new_name)
				if len(flds)>1:
					new_header = ' '.join([new_header]+flds[1:])
				new_headers.append(new_header)
				s = randomizeSequence(seqs[xi])
				#print s.replace(" ",'')
				new_seqs.append(s)
		headers = new_headers
		seqs = new_seqs

	#for xi in range(len(headers)):
	#	print headers[xi]
	#	print seqs[xi]

	for (hdr,seq) in zip(headers,seqs):
		seq = seq.replace(' ','')
		seq = seq.replace('-','')
		
		(name, props) = parseHeader(hdr)
		mutantof = None
		try:
			mutantof = props['mutant.of']
			baseseq = sug_dict[mutantof]
			dnaseq = makeMutantFromSequence(seq, baseseq)
			#print "Used suggestion"
		except KeyError:
			if not mutantof is None:
				raise Exception, "Asked to make mutant of {} but sequence not found in suggestions".format(mutantof)
			dnaseq = translate.randomReverseTranslate(seq)
		#dnaseq = translate.reverseTranslate(seq)
		assert(translate.translate(dnaseq)==seq)
		fullseq = options.prefix + dnaseq + options.suffix
		mutant_seqs[name] = (dnaseq, fullseq)
		#name = biofile.firstField(hdr)
		line = "{name:s}\t{dna:s}\tL={length:d}bp, {desc:s}\n".format(name=name, dna=fullseq, length=len(fullseq), desc=hdr)
		data_outs.write(line)
		n_written += 1

	data_outs.write("\n\n# Confirmation details:\n")
	for (hdr,seq) in zip(headers,seqs):
		(name, props) = parseHeader(hdr)
		(mutant_seq, fullseq) = mutant_seqs[name]
		prot = translate.translate(mutant_seq)
		fullprots = [translate.translateRaw(fullseq[i:]) for i in range(3)]
		allnames = ['core','frame 1','frame 2','frame 3']
		data_outs.write("# {name:s} core and 3-frame full translation\n".format(name=name))
		for (i,p) in enumerate([prot] + fullprots):
			line = "# {name:s} {prot:s}\n".format(prot=p, name=allnames[i])
			data_outs.write(line)
		data_outs.write("#\n")

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

