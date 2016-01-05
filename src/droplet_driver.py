#! python

import sys, os, math, random, argparse
import util, biofile, protprop, stats, scriptutil

'''
Given a sequence, and a set of residue classes, residue volumes, and packing density:
Compute the volume of residues in each class
Compute the radius of 
Compute the outer radius of the sphere that residues would cover, given an inner radius
'''

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input FASTA filename")
	# Optional arguments
	parser.add_argument("--sentinel-species", dest="sentinel_species", type=str, default=None, help="species identifier for sequence")
	parser.add_argument("--start-sequence", dest="start_sequence", type=str, default=None, help="sequence that starts the domain")
	parser.add_argument("--end-sequence", dest="end_sequence", type=str, default=None, help="sequence that ends the domain")
	parser.add_argument("--start-position", dest="start_position", type=int, default=None, help="starting sequence position for search (1-based)")
	parser.add_argument("--end-position", dest="end_position", type=int, default=None, help="ending sequence position for search (1-based, inclusive)")
	parser.add_argument("-d", "--packing-density", dest="packing_density", type=float, default=0.61, help="partial specific volume, proportion of a given volume occupied by amino acids, if packed together")
	parser.add_argument("--query-orf", dest="query_orf", action='append', default=[], help="starting letters of protein systematic name to process")
	parser.add_argument("--query-gene", dest="query_gene", action='append', default=[], help="starting letters of gene name to process")
	parser.add_argument("--translate", dest="translate_sequences", action='store_true', default=False, help="translate sequences?")
	parser.add_argument("--degap", dest="degap", action='store_true', default=False, help="remove gaps from input FASTA file?")
	parser.add_argument("-a", "--aas", dest="aa_groups", action='append', default=[], help="amino-acid groups to quantify as units, e.g. FYW, DE, KRH, CLIMFV, in order from center to outside of droplet")
	parser.add_argument("-g", "--debug", dest="debugging", action='store_true', default=False, help="execute debugging code?")
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

	# Read background information
	aa_volumes = {}
	vols = util.readTable(open(os.path.expanduser('~/research/lib/data/harpaz-aa-volumes.txt')))
	#print vols['volume.a3']
	aa_volumes = dict(zip(vols['aa'], [x/1000.0 for x in vols['mean.volume.a3']]))
	#print aa_volumes

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
	gene_orf_dict = dict([(biofile.secondOrFirstField(h), biofile.firstField(h)) for h in headers])
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

	# Find starting positions
	start_position, end_position = scriptutil.findSequencePositions(options.start_position, options.end_position, 
		options.start_sequence, options.end_sequence, options.sentinel_species, headers, seqs)
	data_outs.write("# start = {:d}, end = {:d}\n".format(start_position, end_position))


	# Select out the subset of sequences
	seqslice = slice(start_position-1,end_position)
	for k in query_keys:
		prot_dict[k] = prot_dict[k][seqslice]

	# Remove gaps?
	if options.degap:
		for k in query_keys:
			prot_dict[k] = prot_dict[k].replace("-",'')
	
	if options.debugging:
		query_keys = query_keys[0:min(len(query_keys,100))]

	# Groups/classes of amino acids to quantify (e.g. A, H, FY, NQ, STY)
	aa_groups = [''.join(sorted(group)) for group in options.aa_groups] # + [x for x in options.single_aas]
	all_aas = ''.join(aa_groups)
	# Write output
	comp = protprop.Composition()
	pprop = protprop.ProteinProperties()
	n_written = 0
	pi = 3.14159265

	def radius_from_volume(v):
		return (v*0.75/pi)**(0.333333)
	
	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('id','Protein identifier','s')
	dout.addHeader('length','Total number of amino acids in sequence','d')
	dout.addHeader('total.volume','Total droplet volume (nm^3)','f')
	dout.addHeader('aas','Amino acids in class','s')
	dout.addHeader('num.aas','Number of amino acids belonging to this class','d')
	dout.addHeader('mean.aa.radius','Mean radius of amino acids in this class (nm)','f')
	dout.addHeader('aa.volume','Volume of material composed of amino acids in this class (nm^3)','f')
	dout.addHeader('cumulative.volume','Total volume of droplet so far (nm^3)','f')
	dout.addHeader('inner.radius','Radius of droplet before adding this class (nm)','f')
	dout.addHeader('outer.radius','Radius of droplet after adding this class (nm)','f')
	dout.addHeader('thickness','Thickness of layer (nm)','f')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	format = dout.getFormat(named=True)
	n_written = 0
	for id in query_keys:
		seq = prot_dict[id]
		cur_radius = 0.0
		cur_volume = 0.0
		total_volume = sum([aa_volumes[aa] for aa in seq])/options.packing_density
		for aas in aa_groups:
			target_aas = [aa for aa in seq if aa in aas]
			vols = [aa_volumes[aa] for aa in target_aas]
			vol = sum(vols)/options.packing_density
			#print vol
			if vol>0.0:
				avg_radius = stats.mean([radius_from_volume(v) for v in vols])
			else:
				avg_radius = 0.0
			new_volume = cur_volume + vol
			new_radius = radius_from_volume(new_volume)
			line = format.format(id=id, length=len(seq), total_volume=total_volume,
				aas=aas, num_aas=len(target_aas), mean_aa_radius=avg_radius, 
				aa_volume=vol, cumulative_volume=new_volume, inner_radius=cur_radius, 
				outer_radius=new_radius, thickness=new_radius-cur_radius)
			data_outs.write(line)
			n_written += 1
			cur_radius = new_radius
			cur_volume = new_volume

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

