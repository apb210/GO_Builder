#!/usr/bin/env python

# build a go model for a homo-oligomeric protein in which no distinction
# is made between interactions within and between subunits
# 

import sys,os,go_lib,getopt

# ======================================================================

Usage="""

Usage:
	go_build_sym_nmer.py [-h] [--help] [-t Tf] [--tf=Tf] [-k key] [--key=key]
	[-g gamma] [-d nongo_dist] [-r repdist] [-s] file.pdb 

	in which:
		Tf is the desired folding temperature (default = 350 K)
		key is a handy key for naming the output files (default = junk)
		nrep is the number of concatenated repeats to make (default = 1)
		gamma is the fraction of the native contact energy assigned to
			non-native contacts (default = 0.0)
		-s indicates that interactions should be symmetrized (i.e. all
		bond and angle minima will be averaged and non-bonded interactions
		will be cumulatively merged (i.e. Ai--Bj => Bi--Aj)
		repdist is a fixed repulsive distance to use for all CA's (rather
			than those taken from the native structure). 
	
	Only the pdb file is a required argument. It must have all backbone heavy
	atoms and amide hydrogens defined, though side-chain heavy atoms are
	used too.

"""

# parse input
# ======================================================================

if (len(sys.argv)==1):
	print Usage
	sys.exit(0)

Tf = 350.0 	# Kelvin
n_repeat = 1
key = "junk"
gamma = 0.0
symm = 0
ngdist = 5.5	# default distance (in A) for non-native CA-CA contacts
repdist = -1.0
optlist, files = getopt.getopt(sys.argv[1:], 'ht:k:g:d:sr:', \
		['help','tf=','key=','gamma=','nongodist=','symmetrize','repdist='])
for opt, val in optlist:
	if opt in [ "-t", "--tf" ]:
		Tf = float(val)
	elif opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in [ "-k", "--key" ]:
		key = val
	elif opt in [ "-s", "--symmetrize" ]:
		symm = 1
	elif opt in [ "-g", "--gamma" ]:
		gamma = float(val)
	elif opt in [ "-d", "--nongodist" ]:
		ngdist = float(val)
	elif opt in [ "-r", "--repdist" ]:
		repdist = float(val)
	else:
		print "\nUNKNOWN OPTION: %s\n\n" % opt
		print Usage
		sys.exit(1)
if len(files) != 1:
	print "\n\nExactly one pdb file must be specified on command line\n\n"
	print Usage
	sys.exit(1)

sys.stdout.write("=========================================================\n")
sys.stdout.write("Building go model dimer from PDB file %s\n" % (files[0]))
sys.stdout.write("Target folding temperature = %8.3f K \n" % (Tf))
sys.stdout.write("Key for naming output files = %s\n" % (key))
if symm:
	sys.stdout.write("Will symmetrize interactions in oligomer\n")
else:
	sys.stdout.write("Will not symmetrize interactions in oligomer\n")
sys.stdout.write("Non-go weighting ('gamma') = %8.3f\n" % (gamma))
sys.stdout.write("Default non-go contact distance = %8.3f\n" % (ngdist))

if gamma < 0.00001:
	outp_base = "go_%s_tf%s" % (key,str(Tf))
else:
	outp_base = "go_%s_tf%s_ng%s" % (key,str(Tf),str(gamma))

if symm:
	outp_base += "_sym"

outp_dir = outp_base
if not os.path.exists(outp_dir):
	os.mkdir(outp_dir)

outp_base = outp_dir + "/" + outp_base

sys.stdout.write("Base for output file names = %s\n" % (outp_base))
sys.stdout.write("=========================================================\n")

info_out = open(outp_base+"_info.dat","w")
info_out.write("Building go model from PDB file %s\n" % (files[0]))
info_out.write("Target folding temperature = %8.3f K \n" % (Tf))
info_out.write("Key for naming output files = %s\n" % (key))
if symm:
	info_out.write("Will symmetrize interactions in dimer\n")
else:
	info_out.write("Will not symmetrize interactions in dimer\n")
info_out.write("Non-go weighting ('gamma') = %8.3f\n" % (gamma))
info_out.write("Base for output file names = %s\n" % (outp_base))
info_out.write("Default non-go contact distance = %8.3f\n" % (ngdist))
info_out.close()

# write out CA-only coords 
# ======================================================================

meta_crd, xyz = go_lib.make_coords(open(files[0]).readlines())
chains = meta_crd.keys()
chains.sort()
nchain = len(chains)

sys.stdout.write("%i chains found in PDB file:\n" % (len(chains)))
for c in chains:
	sys.stdout.write("\t%s\n" % (c))

#if len(chains) >2:
#	sys.stderr.write("\n******************************************************\n")
#	sys.stderr.write("This version can only handle dimers\n")
#	sys.stderr.write("******************************************************\n")
#	sys.stderr.write("%s\n" % (Usage))

go_lib.write_CA_pdb(outp_base + ".pdb", xyz, chains)

# write out sequences
# ======================================================================
seqs = {}
if symm:
	for c in chains:
		seqs[c] = map(lambda x: "G%s"%(x), range(1,len(xyz[c])+1))
	nres_tot = len(xyz[chains[0]])
else:
	off = 0
	nres_tot = 0
	for c in chains:
		seqs[c] = map(lambda x: "G%s"%(x), range(1+off,len(xyz[c])+1+off))
		off += len(xyz[c])
		nres_tot += len(xyz[c])

go_lib.write_sequences(outp_base+"_seq.inp", seqs, key)

# write out topology
# ======================================================================

if symm:
	chain = chains[0]
	nres = len(xyz[chain])
	residues = []
	for i in range(1,nres+1):
		name = meta_crd[chain][i]["name"]
		mass = go_lib.massmap[name]
		residues.append(("G%i"%(i),mass))
else:
	residues = []
	idx = 0
	for c in chains:
		nres = len(xyz[c])
		for i in range(1,nres+1):
			idx += 1
			name = meta_crd[c][i]["name"]
			mass = go_lib.massmap[name]
			residues.append(("G%i"%(idx),mass))

go_lib.write_topology(outp_base+"_top.inp",residues,key)

# calculate parameters
# ======================================================================

link = 0
eps_res = 0.0054*Tf

bonds = {}
angles = {}
dihedrals = {}
nonbonded = {}

for i in range(nchain):
	chaini = chains[i]
	bonds[chaini] = go_lib.calc_bonds(meta_crd[chaini],link)
	angles[chaini] = go_lib.calc_angles(meta_crd[chaini],link)
	dihedrals[chaini] = go_lib.setup_dihe(meta_crd[chaini],link)
	nonbonded[chaini] = {}
	nonbonded[chaini][chaini] = go_lib.calc_ncon_intra(meta_crd[chaini],ngdist,gamma)
	for j in range(i+1,nchain):
		chainj = chains[j]
		nonbonded[chaini][chainj] = go_lib.calc_ncon_inter(meta_crd[chaini],meta_crd[chainj],ngdist,gamma)

#bonds = go_lib.average_bonded_list(bonds)
	
#print nonbonded.keys()
#print nonbonded['A'.keys()
#print nonbonded.keys()

#print "inter"
#x = nonbonded['A   ']['B   ']
#nbx = x['nc']
#for k in nbx.keys():
#	print k[0],k[1],nbx[k][0], nbx[k][1]
#
#print "intra"
#x = nonbonded['A   ']['A   ']
#nbx = x['nc']
#for k in nbx.keys():
#	print k[0],k[1],nbx[k][0], nbx[k][1]
#
#interkeys = nonbonded['A   ']['B   ']['nc'].keys()
#intrakeys = nonbonded['A   ']['A   ']['nc'].keys()
#for k in interkeys:
#	krev = (k[1],k[0])
#	if k in intrakeys or krev in intrakeys:
#		print k
#sys.exit(0)

if symm:
	bonds = go_lib.average_bonded_list(bonds)
	angles = go_lib.average_bonded_list(angles)
	nonbonded =go_lib.average_nonbonded_matrix(nonbonded)
	#dihedrals = go_lib.merge_bonded(dihedrals_a,dihedrals_b)
	dihedrals = dihedrals[chains[0]]
else:
	bonds = go_lib.merge_bonded(bonds_a, bonds_b)
	angles = go_lib.merge_bonded(angles_a, angles_b)
	nonbonded =go_lib.merge_nonbonded(nonbonded_a, nonbonded_b,nonbonded_ab)
	dihedrals = go_lib.merge_bonded(dihedrals_a,dihedrals_b)

go_map,nongo_map, qlist, nat_scale_fac = go_lib.make_go_nongo_map(nonbonded, \
		eps_res, nres_tot)

go_and_nongo = go_lib.mergemap(go_map,nongo_map)
GOMODEL_nbfix = go_lib.gotozero(go_map,nongo_map,nres_tot)

# write linear pdb for each chain.
# these are identical coordinates, so will need to be translated/rotated
# before energy/force evaluation
# ======================================================================

for chain in chains:
	go_lib.write_CA_pdb_linear(outp_base+"_linear_%s.pdb"%(chain.strip()),\
			bonds,angles,chain)


# write out parameters
# ======================================================================
go_lib.write_parms(outp_base+"_parm.inp",key,bonds,angles,dihedrals,
		nonbonded['sig'],go_and_nongo,nat_scale_fac,Tf,repdist)
go_lib.write_parms(outp_base+"_GOMODEL_parm.inp",key,bonds,angles,dihedrals,
		nonbonded['sig'],GOMODEL_nbfix,nat_scale_fac,Tf,repdist)
go_lib.write_golist(outp_base+"_GOMODEL_golist.dat", go_map, nat_scale_fac,
		nres_tot, n_repeat)

# write 'qdetails' file
# ======================================================================
#go_lib.write_qdetails(outp_base+"_qdetails.dat",meta_crd[chains[0]],nonbonded['details'],
#		go_map)

#go_lib.write_qlist(outp_base+"_qlist.dat",qlist,nres_tot,n_repeat)
