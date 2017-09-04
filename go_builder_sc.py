#!/usr/bin/env python

import sys,os,go_lib, getopt

# ======================================================================

Usage="""

	Builds Karanicolas-Brooks Go model [ + SIDE-CHAINS!!! ]

Usage:
	go_builder_sc.py [-h] [--help] [-t Tf] [--tf=Tf] [-k key] [--key=key]
	[-g gamma] [-n nrep] [-d nongo_dist] [--nrepeat=nrep] [-l list] 
	[-e energies.dat] [--ncon] [-w N|--water N] [-S scfac] [-B bbfac] [-G]
	[-S]
	file.pdb 

	in which:
		Tf is the desired folding temperature (default = 350 K)
		key is a handy key for naming the output files (default = junk)
		nrep is the number of concatenated repeats to make (default = 1)
		gamma is the fraction of the native contact energy assigned to
			non-native contacts (default = 0.0)
		nongo_dist is the sigma to use for non-go contacts when gamma > 0.0
		list is a list of go contacts and relative weights to be 
			used instead of the list that would usually be
			constructed.
		energies.dat is a list of relative interaction energies between
			different residues that will be used to scale native
			contact strengths (i1, j1, eij1; i2, j2, eij2; ...)
		--ncon will use the number of all-atom contacts to weight 
		side-chain interactions
		--water=N | -w N will put the protein in the centre of a
		box of "water" made of NxN cubes of %8.3f A each containing
		%i W4 particles, representing four waters (an equilibrated
		box I happen to have!). If you want to use a different size
		box that is relatively simple to do by hand.
		-S scfac will scale side-chain interactions by scfac
		-B bbfac will scale backbone interactions by bbfac
		-G use geometry to define hydrogen bonding
		-A write statistical angle potentials
	
	ONLY THE PDB FILE IS A REQUIRED ARGUMENT. It must have all backbone heavy
	atoms and amide hydrogens defined, though side-chain heavy atoms are
	used too.

""" % (go_lib.std_boxlen, go_lib.std_nwat)

# parse input
# ======================================================================

if (len(sys.argv)==1):
	print Usage
	sys.exit(0)

Tf = 350.0 	# Kelvin
n_repeat = 1
regularize = 1
key = "junk"
use_list = 0
gamma = 0.0
average = 0
symm = 0
water = 0
defbon = 0
scscale = 1.0
bbscale = 1.0
ncon_weight = 0
geomh = 0
stata = 0
charmmpdb = 1
energyfile = "none"
ngdist = 5.5	# default distance (in A) for non-native CA-CA contacts
optlist, files = getopt.getopt(sys.argv[1:], 'ht:k:n:g:d:l:w:sS:B:e:GAp', \
		['help','tf=','key=','nrepeat=','gamma=','list=',\
		'symmetrize','nongodist=','water=','scscale=','energies=',\
		'defbon','ncon','geomh','stat-angle'])
for opt, val in optlist:
	if opt in [ "-t", "--tf" ]:
		Tf = float(val)
	elif opt in [ "-h", "--help" ]:
		print Usage
		sys.exit(0)
	elif opt in [ "-k", "--key" ]:
		key = val
	elif opt in [ "-n", "--nrepeat" ]:
		n_repeat = int(val)
	elif opt in [ "-s", "--symmetrize" ]:
		symm = 1
	elif opt in [ "-g", "--gamma" ]:
		gamma = float(val)
	elif opt in [ "-d", "--nongodist" ]:
		ngdist = float(val)
	elif opt in [ "-l", "--list" ]:	
		go_list_file = val
		use_list = 1
	elif opt in [ "-w", "--water" ]:	
		water = int(val)
	elif opt in [ "-p" ]:	
		charmmpdb = 0
	elif opt in [ "-S", "--scscale" ]:	
		scscale = float(val)
	elif opt in [ "-B", "--bbscale" ]:	
		bbscale = float(val)
	elif opt in [ "-G", "--geomh" ]:	
		geomh = 1
	elif opt in [ "-A", "--stat-angle" ]:	
		stata = 1
	elif opt in [ "-e", "--energies" ]:	
		energyfile = val
	elif opt in ["--defbon" ]:	
		defbon = 1
	elif opt in ["--ncon" ]:	
		ncon_weight = 1
	else:
		print "\nUNKNOWN OPTION: %s\n\n" % opt
		print Usage
		sys.exit(1)
if len(files) < 1:
	#print "\n\nExactly one pdb file must be specified on command line\n\n"
	print Usage
	sys.exit(1)
elif len(files) > 1:
	average = 1

# create output directories and base file names
# ======================================================================

#if gamma < 0.00001:
#	if n_repeat == 1:
#		outp_base = "go_sidec_%s_tf%s" % (key,str(Tf))
#	else:
#		outp_base = "go_sidec_%s_tf%s_nr%i" % (key,str(Tf),n_repeat)
#else:
#	if n_repeat == 1:
#		outp_base = "go_sidec_%s_tf%s_ng%s" % (key,str(Tf),str(gamma))
#	else:
#		outp_base = "go_sidec_%s_tf%s_ng%s_nr%i" % (key,str(Tf),str(gamma),n_repeat)

if gamma < 0.00001:
	if n_repeat == 1:
		outp_base = "go_sidec_%s" % (key)
	else:
		outp_base = "go_sidec_%s_nr%i" % (key,n_repeat)
else:
	if n_repeat == 1:
		outp_base = "go_sidec_%s_ng%s" % (key,str(gamma))
	else:
		outp_base = "go_sidec_%s_ng%s_nr%i" % (key,str(gamma),n_repeat)

if symm:
	outp_base += "_sym"
if water>0:
	outp_base += "_w"

outp_dir = outp_base
if not os.path.exists(outp_dir):
	os.mkdir(outp_dir)

outp_base = outp_dir + "/" + outp_base

# print handy information to *_info.dat file and stdout
# ======================================================================

info = []
info.append("=========================================================\n")
command = reduce(lambda x, y: x+" "+y,sys.argv)
command = command[command.find("go_builder_sc.py"):]
info.append("This go model generated with the command\n%s\n" % (command))
info.append("=========================================================\n")
info.append("Building go model from PDB file %s\n" % (files[0]))
info.append("Target folding temperature = %8.3f K \n" % (Tf))
info.append("Key for naming output files = %s\n" % (key))
info.append("Number of NC-linked repeats = %i\n" % (n_repeat))
info.append("Non-go weighting ('gamma') = %8.3f\n" % (gamma))
info.append("Default non-go contact distance = %8.3f\n" % (ngdist))
info.append("Scaling side-chain contacts by %8.3f\n" % (scscale))
info.append("Scaling backbone contacts by %8.3f\n" % (bbscale))
if symm:
	info.append("Will attempt to symmetrize interactions in oligomer\n")
if average:
	info.append("Will average over multiple structures\n")
else:
	info.append("Will build model from single structure\n")
if water>0:
	info.append("Will place protein in water box of %ix%i cubes of\n"%(water,water))
	info.append("dimension %8.3f A containing %i pseudo waters\n"\
			%(go_lib.std_boxlen,go_lib.std_nwat))
	info.append("i.e. total box side = %8.3f A\n" % (float(water)*go_lib.std_boxlen))

if use_list:
	info.append("Taking go contacts from the following file = %s\n" \
			% (go_list_file))

info.append("Base for output file names = %s\n" % (outp_base))
info.append("=========================================================\n")

info_out = open(outp_base+"_info.dat","w")
for i in info:
	sys.stdout.write("%s" %(i))
	info_out.write("%s" %(i))
info_out.close()

# parse pdb and write out CA-only coords 
# ======================================================================

# check if NMR file
#pdblines = open(files[0]).readlines()
#if len(filter(lambda x: x.find("MODEL")==0, pdblines)) != 0:
#	nmrfile = 1
#else:
#	nmrfile = 0
#
#models = go_lib.split_NMR_pdb(pdblines)
#print models

if average:
	meta_list = []
	xyz_list = []
	for file in files:
		meta_crd, xyz = go_lib.make_coords(open(file).readlines(),charmmpdb)
		meta_list.append(meta_crd)
		xyz_list.append(xyz)

	#meta_crd, xyz = go_lib.average_coords(meta_list,xyz_list)
else:
	meta_crd, xyz = go_lib.make_coords(open(files[0]).readlines(),charmmpdb)

chains = meta_crd.keys()
chains.sort()
nchain = len(chains)

sys.stdout.write("Chains found in PDB file:\n")
for c in chains:
	sys.stdout.write("\t%s\n" % (c))

if len(chains) >1 and n_repeat>1:
	sys.stdout.write("\n******************************************************\n")
	sys.stdout.write("WARNING!\n")
	sys.stdout.write("There is more than one chain in your file and you are\n")
	sys.stdout.write("building an N-C linked repeat protein;\n")
	sys.stdout.write("only the first chain will be used!\n")
	sys.stdout.write("******************************************************\n")

if (n_repeat==1):
	oxyz = xyz
	ochains = chains
else:
	oxyz = { chains[0]: go_lib.make_CA_nclink(xyz,n_repeat) }
	ochains = [ chains[0] ]

go_lib.write_CA_SC_pdb(outp_base + ".pdb", oxyz, meta_crd, ochains, regularize)

if symm:
	# all chains must AT LEAST be the same length!
	l0 = len(xyz[chains[0]])
	for chaini in chains:
		li = len(xyz[chaini])
		for chainj in chains:
			lj = len(xyz[chainj])
			if li != lj:
				sys.stderr.write("  !!! ERROR !!!  \n")
				sys.stderr.write("Cannot symmetrize contacts!\n")
				sys.stderr.write("Length of chain %s = %i\n"\
						% ( chaini, li ))
				sys.stderr.write("Length of chain %s = %i\n"\
						% ( chainj, lj ))
				sys.exit(1)

# solvate protein and write coords if necc.
# ======================================================================

if water > 0:
	# (i) centre protein coordinates
	xyz,meta_crd = go_lib.center_coords(xyz,meta_crd)
	#go_lib.write_CA_pdb(outp_base + "_center.pdb", xyz, chains)
	# (ii) create BIG water box!
	wbox = go_lib.make_wbox(water)
	# (iii) attempt to solvate by deleting waters
	go_lib.solvate(xyz,wbox)
	go_lib.write_CA_SC_pdb(outp_base + "_solv.pdb", xyz, meta_crd, \
			xyz.keys(), regularize)


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

if water == 0:
	nsolv = 0
else:
	nsolv = len(xyz["SOLV"])

go_lib.write_sequences(outp_base+"_seq.inp", seqs, "", key, nsolv)

# write out topology
# ======================================================================
if symm:
	chain = chains[0]
	nres = len(xyz[chain])
	residues = []
	for i in range(1,nres+1):
		name = meta_crd[chain][i]["name"]
		#mass = go_lib.massmap[name]
		residues.append(("G%i"%(i),name))
else:
	residues = []
	idx = 0
	for c in chains:
		nres = len(xyz[c])
		for i in range(1,nres+1):
			idx += 1
			name = meta_crd[c][i]["name"]
			#mass = go_lib.massmap[name]
			#residues.append(("G%i"%(idx),mass))
			residues.append(("G%i"%(idx),name))

#chain = chains[0]
#nres = len(xyz[chain])
#residues = []
#for i in range(1,nres+1):
#	name = meta_crd[chain][i]["name"]
#	mass = go_lib.massmap[name]
#	residues.append(("G%i"%(i),mass))

go_lib.write_sc_topology(outp_base+"_top.inp",residues,key)

# calculate parameters
# ======================================================================

if n_repeat > 1:
	link = 1
else:
	link = 0

eps_res = go_lib.big_fat_fudge_factor*Tf

bonds = go_lib.calc_sc_bonds(meta_crd[chains[0]],link)
angles = go_lib.calc_sc_angles(meta_crd[chains[0]],link,stata)
dihedrals = go_lib.setup_sc_dihe(meta_crd[chains[0]],link)
nonbonded = go_lib.calc_sc_ncon_intra(meta_crd[chains[0]],ngdist,gamma,geomh)
if energyfile != "none":	# read relative contact energies from file
	nonbonded['nc'] = go_lib.ene_sub(nonbonded['nc'],energyfile)
if ncon_weight == 1:
	nonbonded['nc'] = go_lib.ncon_sub(nonbonded['nc'],nonbonded['ncon'])

if use_list:
	nonbonded['nc'] = go_lib.make_nc_from_list(meta_crd[chains[0]],go_list_file)
#go_map,nongo_map, qlist, nat_scale_fac = go_lib.make_go_nongo_map(nonbonded, \
#		eps_res, nres,scscale,bbscale)
#go_and_nongo = go_lib.mergemap(go_map,nongo_map)
#GOMODEL_nbfix = go_lib.gotozero(go_map,nongo_map,nres)
atomind = go_lib.compute_indices(meta_crd[chains[0]])
nbfix, GOMODEL_nbfix, GOMODEL_golist, qlist, nat_scale_fac, allqlist = go_lib.sc_mergemaps(nonbonded,
		eps_res, nres,atomind,scscale,bbscale)
sig_rep = go_lib.make_sc_sigrep(meta_crd[chains[0]])
if water > 0:
	winte = go_lib.calc_nbond_w_sc(residues)
else:
	winte = {}


# write linear pdb for each chain.
# these will probably overlap, so will need to be translated/rotated
# before energy/force evaluation
# ======================================================================

#for chain in chains:
#	tmpxyz = { chain: go_lib.make_CA_linear(bonds,angles) }
#	go_lib.write_CA_pdb(outp_base+"_linear_%s.pdb"%(chain),tmpxyz, [chain])

# write out parameters
# ======================================================================
go_lib.write_sc_parms(outp_base+"_parm.inp",key,bonds,angles,dihedrals,
		sig_rep,nbfix,nat_scale_fac,Tf,winte)
go_lib.write_sc_parms(outp_base+"_gomodel_parm.inp",key,bonds,angles,dihedrals,
		sig_rep,GOMODEL_nbfix,nat_scale_fac,Tf,winte)
go_lib.write_sc_golist(outp_base+"_gomodel_golist.dat", GOMODEL_golist, nat_scale_fac,
		nres, n_repeat)

# write 'qdetails' file
# ======================================================================
#go_lib.write_qdetails(outp_base+"_qdetails.dat",meta_crd[chains[0]],nonbonded['details'],
#		go_map)

# write GROMACS topology
# ======================================================================

go_lib.write_qlist(outp_base+"_qlist.dat",qlist,nres,n_repeat)
go_lib.write_qlist(outp_base+"_allqlist.dat",allqlist,nres,n_repeat)

sys.exit(0)
#if n_repeat > 1:
#	go_lib.write_qlist_intra(outp_base+"_intra_qlist.dat",qlist,nres,n_repeat)
#	go_lib.write_qlist_inter(outp_base+"_inter_qlist.dat",qlist,nres,n_repeat)
#	for r in range(n_repeat):
#		go_lib.write_module_qlist("%s_module%i_qlist.dat" % (outp_base,r),
#				qlist,nres,r)
