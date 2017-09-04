#!/usr/bin/env python

import sys,os,go_lib,generic_lib,getopt

# build a generic coarse model based on a sequence and specified 
# pair potential
# ======================================================================

Usage="""

Usage:
	generic_builder.py [-h] [--help] [-k key] [--key=key]
	[-p potential] [-a] [-q] seq1 seq2 seq3 ...

	in which:
	        -q add charges
		-k key is a handy key for naming the output files (default = junk)
		-s seq gives the sequence to build
		-S scale factor to apply to pair potential (default=1.0)
		-a harmonic angle potential (default: statistical potential)
		-p : the pair potential to use, either 'mj' (default) or 'bln' 

	Only the sequence is a required argument. 

""" 

# parse input
# ======================================================================

if (len(sys.argv)==1):
	print Usage
	sys.exit(0)

defb = 3.81
defa = 130.0	# extended
defr = 6.0
sequence = ""
pairpot = "mj"
pairscale = 1.0
key = "junk"
addcharges = 0
uniform_mass = 0
genangle = 1
optlist, rawsequences = getopt.getopt(sys.argv[1:], 'hqk:p:a', \
		['help','key=','pairpot='])

usechains = []
for opt, val in optlist:
	if opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in [ "-k", "--key" ]:
		key = val
	elif opt in [ "-q", "--charge" ]:
		addcharges = 1
	#elif opt in [ "-s", "--seq" ]:
	#	seq = val.upper()
	elif opt in [ "-S", "--pairscale" ]:
		pairscale = val
	elif opt in [ "-a" ]:
		genangle = 0
	elif opt in [ "-p", "--pairpot" ]:
		pairpot = val
	else:
		print "\nUNKNOWN OPTION: %s\n\n" % opt
		print Usage
		sys.exit(1)

if len(rawsequences) == 0:
	print Usage
	sys.exit(1)

if pairpot == "kh":
	addcharges = 1

sequences = map(lambda x: x.upper(), rawsequences)
#sequence = sequences[0].upper()

# create output directories and base file names
# ======================================================================

outp_base = "%s_%s" % (pairpot,key)

outp_dir = outp_base
if not os.path.exists(outp_dir):
	os.mkdir(outp_dir)

outp_base = outp_dir + "/" + outp_base

# print handy information to *_info.dat file and stdout
# ======================================================================

info = []
info.append("=========================================================\n")
command = reduce(lambda x, y: x+" "+y,sys.argv)
command = command[command.find("generic_builder.py"):]
info.append("This generic model generated with the command:\n\"%s\"\n" % (command))

info.append("Building generic model from sequences:\n")
for s in sequences:
	info.append("\t%s\n" % (s))

info.append("Key for naming output files = %s\n" % (key))
if pairpot == "mj":
	info.append("Using simple miyazawa-jernigan pair potential\n")
if pairpot == "kh":
	info.append("Using kim-hummer potential\n")
elif pairpot == "mjbln":
	info.append("Using BLN pair potential derived from mj\n")
elif pairpot == "bln":
	info.append("Using BLN pair potential\n")

info.append("Base for output file names = %s\n" % (outp_base))
info.append("=========================================================\n")

info_out = open(outp_base+"_info.dat","w")
for i in info:
	sys.stdout.write("%s" %(i))
	info_out.write("%s" %(i))
info_out.close()

# translate one letter code into mj or bln residues
residue_seqs = map(lambda y: map(lambda x: generic_lib.resmap[x], y), sequences)

generic_lib.write_sequence(outp_base+"_seq.inp", key, residue_seqs)

generic_lib.write_topology(outp_base+"_top.inp", addcharges, uniform_mass)

generic_lib.write_parms(outp_base+"_parm.inp",defb,defa,defr,pairpot,genangle,residue_seqs)

crd_dat, xyz = generic_lib.make_sequ_coords(residue_seqs,defb,defa)
generic_lib.write_CA_pdb(outp_base+".pdb", xyz, crd_dat, ['A'])


sys.exit(0)












