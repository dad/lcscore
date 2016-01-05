#! python

import sys, os, math, random, argparse
import util, biofile, muscle, translate

'''
Prune out:
1) 
Choose 
'''

def shorten(x):
	flds = x.split(' ')
	res = '{g}.{s}'.format(g=flds[0][0], s='_'.join(flds[1:]))
	return res

def stringdiff(s1, s2):
	if len(s1) != len(s2):
		return [-1]
	diffs = []
	for (xi, s) in enumerate(s1):
		if s2[xi] != s:
			diffs.append((xi,s,s2[xi]))
	return diffs

def dedupe(headers, sequences, fn=biofile.firstField):
	assert len(headers) == len(sequences)
	new_headers = []
	new_seqs = []
	ids = [fn(h) for h in headers]
	for (xi, id) in enumerate(ids):
		dupe_ids = [yi for yi in range(len(ids)) if ids[yi]==id]
		# De-gap the duplicate-id sequences
		seqs = list(set([sequences[yi].replace('-','') for yi in dupe_ids]))
		if len(seqs)>1:
			#print id, len(seqs)
			# test to see whether sequences are identical.
			for i in range(len(seqs)-1):
				dupe = False
				for j in range(i+1,len(seqs)):
					#print "({},{}) ".format(i,j),
					if seqs[i] == seqs[j]:
						dupe = True
				if not dupe:
					#print "not a dupe", len(seqs[0]), len(seqs[1])
					#print stringdiff(seqs[0], seqs[1])
					#if len(seqs[0]) != len(seqs[1]):
					#	print seqs[0]
					#	print seqs[1]
					new_headers.append(headers[dupe_ids[i]])
					new_seqs.append(sequences[dupe_ids[i]])
		else:
			new_headers.append(headers[xi])
			new_seqs.append(sequences[xi])
	return new_headers, new_seqs


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Prune FASTA file to a single sequence per species")
	# Required arguments
	parser.add_argument(dest="in_fname", help="input FASTA filename")

	# Optional arguments
	parser.add_argument("--orthodb-names", dest="in_names_fname", default=None, help="naming filename")
	parser.add_argument("-i", "--include", dest="include_identifiers", action='append', default=[], help="identifiers to look for (OR logic)")
	parser.add_argument("-x", "--exclude", dest="exclude_identifiers", action='append', default=[], help="identifiers to exclude (AND logic)")
	parser.add_argument("-s", "--sequence", dest="require_sequence", action='append', default=[], help="sequences to look for ")

	parser.add_argument("--sentinel", dest="sentinels", action='append', default=[], help="identifiers to absolutely retain")
	parser.add_argument("--remove-x", dest="remove_x", action='store_true', default=False, help="remove sequences containing X's?")
	parser.add_argument("--min-length", dest="minimum_length", type=int, default=None, help="minimum length")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)
	if not options.fasta_out_fname is None:
		outf = file(options.fasta_out_fname,'w')
		fasta_outs.addStream(outf)
	else:
		# By default, write to stdout
		fasta_outs.addStream(sys.stdout)

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
		info_outs.write("# Read {:d} sequences\n".format(len(headers)))

	if not options.in_names_fname is None:
		if not os.path.isfile(options.in_names_fname):
		 	raise IOError("# Error: file {} does not exist".format(options.in_names_fname))
		with open(options.in_names_fname,'r') as inf:
			species = util.readTable(inf, header=True)
		
		species_name_lookup = dict(zip(species['ODB_code'], [shorten(x) for x in species['Organism']]))
		new_headers = []
		new_seqs = []
		for (hdr, seq) in zip(headers,seqs):
			odb_code = hdr.strip().split()[-1]
			if len(odb_code) == 5:
				try:
					# Rename
					species_name = species_name_lookup[odb_code]
					new_headers.append('{} {}'.format(species_name, hdr))
				except KeyError:
					new_headers.append(hdr)
			else:
				new_headers.append(hdr)
			new_seqs.append(seq)
		headers = new_headers
		seqs = new_seqs

	(headers, seqs) = zip(*[(h,seqs[xi]) for (xi,h) in enumerate(headers) if (not 'X' in seqs[xi] and not '_' in seqs[xi])])
	info_outs.write("# Pruned X and _; {:d} sequences\n".format(len(headers)))

	# De-dupe the sequences
	#(headers, seqs) = dedupe(headers, seqs)

	# Find sentinels, and save them
	sentinel_hdrs = []
	sentinel_seqs = []
	new_headers = []
	new_seqs = []
	for sent in options.sentinels:
		for (xi, h) in enumerate(headers):
			if sent in h:
				sentinel_hdrs.append(h)
				sentinel_seqs.append(seqs[xi])

	# Reassemble the rest of the sequences, without the sentinels
	(headers, seqs) = zip(*[(h,seqs[xi]) for (xi,h) in enumerate(headers) if not h in sentinel_hdrs])
	info_outs.write("# Removed sentinels; {:d} sequences\n".format(len(headers)))
	# Filter for length
	if not options.minimum_length is None:
		(headers, seqs) = zip(*[(h,seqs[xi]) for (xi,h) in enumerate(headers) if len(seqs[xi]) >= options.minimum_length])
		info_outs.write("# Pruned minimum length; {:d} sequences\n".format(len(headers)))



	# Now go through headers, find multiples, and select one from each.
	selected_indices = [] # index into headers and sequences
	new_headers = []
	new_seqs = []

	orthodb_ids = [h.split()[-1] for h in headers]
	species_names = [h.split()[0].split('_')[0] for h in headers]
	# Sentinel species means we automatically reject all other sequences from this species
	sentinel_species_names = [h.split()[0].split('_')[0] for h in sentinel_hdrs]
	for species_name in list(set(species_names)):
		if species_name in sentinel_species_names:
			# Skip it -- we'll add the sentinel sequence and header later.
			continue
		dupe_indices = [xi for (xi,spec) in enumerate(species_names) if spec==species_name]
		#print species_name, len(dupe_indices)
		if len(dupe_indices)==1:
			# No problem, no duplicate
			selected_indices.append(dupe_indices[0])
			continue
		# Duplicates of the same OrthoDB ID?
		odbids = [orthodb_ids[xi] for xi in dupe_indices]
		odbids_set = list(set(odbids))
		if len(odbids)>1:
			# Duplicate sequences from the same strain/species.
			# Check for header strings
			included_indices = []
			for iid in options.include_identifiers:
				hdrs = [(headers[xi],xi) for xi in dupe_indices]
				iid_hdrs = [xi for (h,xi) in hdrs if iid in h]
				included_indices += iid_hdrs
			excluded_indices = []
			for iid in options.exclude_identifiers:
				hdrs = [(headers[xi],xi) for xi in dupe_indices]
				iid_hdrs = [xi for (h,xi) in hdrs if iid in h]
				excluded_indices += iid_hdrs
			if len(included_indices)>0:
				satisfactory_indices = included_indices
			else:
				satisfactory_indices = dupe_indices[:]
			satisfactory_indices = list(set(satisfactory_indices).difference(set(excluded_indices)))
			#print species_name, len(dupe_indices), len(included_indices), len(excluded_indices), len(satisfactory_indices)
			#print species_name, len(satisfactory_indices)
			if len(satisfactory_indices)>0:
				selected_indices += satisfactory_indices #.append(satisfactory_indices[0])
			else:
				print "Excluded everything for ", species_name

	(headers, seqs) = zip(*[(headers[i],seqs[i]) for i in selected_indices])
	headers = list(headers)
	seqs = list(seqs)
	info_outs.write("# Pruned duplicate-ID sequences; {:d} sequences\n".format(len(headers)))
	if len(options.require_sequence)>0:
		new_headers = []
		new_seqs = []
		for (h,seq) in zip(headers,seqs):
			include_sequence = True
			for s in options.require_sequence:
				#print s, seq.find(s)
				include_sequence = include_sequence and seq.find(s)>-1
			if include_sequence:
				new_headers.append(h)
				new_seqs.append(seq)
		headers = new_headers
		seqs = new_seqs
	info_outs.write("# Pruned required sequences; {:d} sequences\n".format(len(headers)))

	###
	# Done with filtering
	###

	# Add back in the sentinels
	headers += sentinel_hdrs
	seqs += sentinel_seqs
	info_outs.write("# Added back in sentinels; {:d} sequences\n".format(len(headers)))

	# Streamline the headers: add leading ID without strain information
	new_headers = []
	for h in headers:
		species_name = h.split()[0]
		short_species_name = h.split()[0].split('_')[0]
		if species_name != short_species_name:
			new_header = '{} {}'.format(short_species_name, h)
		else:
			new_header = h
		new_headers.append(new_header)
	headers = new_headers

	# Write output
	n_written = 0
	for (hdr, seq) in zip(headers,seqs):
		line = ">{hdr}\n{seq}\n".format(hdr=hdr, seq=seq)
		fasta_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

