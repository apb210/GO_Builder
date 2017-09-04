#!/usr/bin/env python

import sys,os,math

statangle = [106.4, 87.7, 26.3, 130.0, 0.08,  5.0]

#DATA_DIR = "/Users/best/Programs/Scratch_cvs/go_builder"
DATA_DIR = "/home/rbb24/sw/scratch_cvs/go_builder"
res2charge = { 'ALA': 0, 'GLY': 0, 'THR': 0, 'TYR':0,
	  'VAL':0, 'LEU':0, 'ILE':0, 'TRP':0,
	  'GLU':-1, 'ASP':-1, 'SER':0, 'ASN':0,
	  'GLN':0, 'PRO':0, 'PHE':0, 'ARG':+1,
	  'CYS':0, 'HIS':+0.5, 'LYS':+1, 'MET':0 }

res2khsig = { 'ALA': 5.0, 'GLY': 4.5, 'THR': 5.6, 'TYR': 6.5,
	  'VAL': 5.9, 'LEU': 6.2, 'ILE': 6.2, 'TRP': 6.8,
	  'GLU': 5.9, 'ASP':5.6, 'SER': 5.2, 'ASN': 5.7,
	  'GLN': 6.0, 'PRO': 5.6, 'PHE': 6.4, 'ARG': 6.6,
	  'CYS': 5.5, 'HIS': 6.1, 'LYS': 6.4, 'MET': 6.2 }

aalist = res2charge.keys()

resmap = { 'A':'ALA', 'G':'GLY', 'T':'THR', 'Y':'TYR',
	  'V':'VAL', 'L':'LEU', 'I':'ILE', 'W':'TRP',
	  'E':'GLU', 'D':'ASP', 'S':'SER', 'N':'ASN',
	  'Q':'GLN', 'P':'PRO', 'F':'PHE', 'R':'ARG',
	  'C':'CYS', 'H':'HIS', 'K':'LYS', 'M':'MET' }
#resmap['bln'] = { 'A':'B', 'G':'N', 'T':'L', 'Y':'B',
#	  'V':'B', 'L':'B', 'I':'B', 'W':'B',
#	  'E':'L', 'D':'L', 'S':'N', 'N':'L',
#	  'Q':'L', 'P':'N', 'F':'B', 'R':'L',
#	  'C':'B', 'H':'L', 'K':'L', 'M':'B' }
#
#resmap['mj'] = { 'A':'ALA', 'G':'GLY', 'T':'THR', 'Y':'TYR',
#	  'V':'VAL', 'L':'LEU', 'I':'ILE', 'W':'TRP',
#	  'E':'GLU', 'D':'ASP', 'S':'SER', 'N':'ASN',
#	  'Q':'GLN', 'P':'PRO', 'F':'PHE', 'R':'ARG',
#	  'C':'CYS', 'H':'HIS', 'K':'LYS', 'M':'MET' }

blnmap = { 'ALA':'B', 'GLY':'N', 'THR':'L', 'TYR':'B',
	  'VAL':'B', 'LEU':'B', 'ILE':'B', 'TRP':'B',
	  'GLU':'L', 'ASP':'L', 'SER':'N', 'ASN':'L',
	  'GLN':'L', 'PRO':'N', 'PHE':'B', 'ARG':'L',
	  'CYS':'B', 'HIS':'L', 'LYS':'L', 'MET':'B' }

massmap = { 'ALA': 71.0, 'GLY': 57.0,'THR': 101.0,'TYR': 163.0,
	'VAL': 99.0, 'LEU': 113.0,'ILE': 113.0,'TRP': 186.0,
	'GLU': 128.0, 'ASP': 114.0,'SER': 87.0,'ASN': 114.0,
	'GLN': 128.0, 'PRO': 97.0,'PHE': 147.0,'ARG': 157.0,
	'CYS': 103.0, 'HIS': 138.0,'LYS': 128.0,'MET': 131.0 }

def read_karanicolas_dihe():
	"""read in and parse Karanicolas parameters for CA-CA-CA-CA pseudo-dihedral
	angles
	REF: Karanicolas&Brooks, Prot. Sci., 11, 2351-2361"""
	inp = map(lambda x: x.split(), open("%s/karanicolas_dihe_parm.dat"%(DATA_DIR)).readlines())
	kdihe = {}
	for line in inp:
		if line[0] not in kdihe.keys():
			kdihe[line[0]] = {}
		if line[1] not in kdihe[line[0]].keys():
			kdihe[line[0]][line[1]] = []
		kdihe[line[0]][line[1]].append((float(line[2]),int(line[3]),float(line[4])))
	return kdihe

def read_miyazawa_jernigan():
	"""read in and parse attractive Miyazawa-Jernigan parameters
	REF: Miyazawa&Jernigan, JMB, 256, 623-644 (1996)"""
	mjdat = map(lambda x: x.split(), open("%s/miyazawa_jernigan.dat"%(DATA_DIR)).readlines())
	mjmap = {}
	ave_mj = 0.0
	for d in mjdat:
		i,j,mij = d[0],d[1],float(d[2])
		if i not in mjmap.keys():
			mjmap[i] = {}
		if j not in mjmap.keys():
			mjmap[j] = {}
		mjmap[i][j] = mij
		mjmap[j][i] = mij
		ave_mj += mij
	ave_mj = abs(ave_mj/float(len(mjdat)))
	return mjmap,ave_mj

def write_sequence(filename,key,seqs):
	oseq = open(filename,"w")
	oseq.write("* This CHARMM .seq file describes a generic model of %s\n*\n\n"%(key))
	alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	for s in range(len(seqs)):
		seq = seqs[s]
		chain = alpha[s]
		nres = len(seq)
		oseq.write("read sequence card\n")
		oseq.write("* Sequence \n*\n%i\n" % (nres))
		oline = ""
		for r in seq:
			oline += "%s " % (r)
			if len(oline) > 70:
				oseq.write("%s\n"%(oline))
				oline = ""
			#oseq.write( "%s " % (r))
			#cres += 1
			#if cres % 8 == 0:
			#	oseq.write("\n")
		if len(oline) > 0:
			oseq.write("%s\n"%(oline))
		oseq.write("\n\n")
			
		oseq.write("generate %s setup\n\n"%(chain))

	oseq.close()

def write_topology(filename,addcharges,uniform_mass):
	"""write topology for streaming into CHARMM
	residues is a list of residues and their masses to write"""
	otop = open(filename,"w")
	otop.write("* This CHARMM .top file describes a generic model\n")
	otop.write("*\n")
	otop.write("\n")
	otop.write("read rtf card\n")
	otop.write("* Topology for generic model\n")
	otop.write("*\n")
	otop.write("   20   1\n")

	i = 0
	for aa in aalist:
		i+=1
		mass = massmap[aa]
		otop.write("MASS %-5i %4s  %12.6f\n" % (i,aa,mass))
	
	otop.write("\nDECL +CA\n\nAUTOGENERATE ANGLES DIHEDRAL\n\n")

	for aa in aalist:
		if addcharges:
			charge = res2charge[aa]
		else:
			charge = 0.0
		otop.write("RESI %4s    %5.3f\n" % (aa,charge))
		otop.write("GROU\n")
		otop.write("Atom CA %4s   %5.3f\n" % (aa,charge))
		otop.write("Bond CA +CA\n\n")

	otop.write("PRES DISU    0.0\n")
	otop.write("BOND 1CA 2CA\n\n")

	otop.write("END\n\n")
	otop.close()

def write_parm_header(filep,pairpot):
	"""write header of parameter file"""
	filep.write("* This CHARMM .param file describes a generic model\n")
	if pairpot == 'mj':
		filep.write("* using the miyazawa-jernigan pair potential\n")
	elif pairpot == 'kh':
		filep.write("* using the kim-hummer pair potential\n")
	elif pairpot == 'mjbln':
		filep.write("* using a MJ-BLN pair potential\n")
	elif pairpot == 'bln':
		filep.write("* using a BLN pair potential\n")
	filep.write("\n\nread param card\n")
	filep.write("* Parameters for Generic model\n\n")
	
def write_bonded(outp,defb,defa,eps_res,genangle,seqs):
	"""write bonded terms only to parameter file"""
	kbond = eps_res * 200.0
	kangle = eps_res * 40.0
	dihe_fac = eps_res * 0.4
	# for bonds, write all possible combinations 
	outp.write("BOND\n")
	naa = len(aalist)
	for i in range(naa):
		ai = aalist[i]
		for j in range(i+1):
			aj = aalist[j]
			outp.write("%-6s %-6s %12.6f %12.6f\n" % (ai,aj,kbond,defb))
	outp.write("\n")
	# restrict angles to only those that are necessary for this system
	alist = []
	for i in range(naa):
		alist.append(['CYS','CYS',aalist[i]])  # for disulfides
	for s in seqs:
		for i in range(len(s)-2):
			si = s[i]
			sj = s[i+1]
			sk = s[i+2]
			if [si,sj,sk] not in alist and [sk,sj,si] not in alist:
				alist.append([si,sj,sk])
	# for disulfides

	outp.write("ANGLE\n")
	if genangle: #generic angles, default
		for a in alist:
			outp.write("%-4s %-4s %-4s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" \
					% tuple(a+statangle))
		else:
			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" \
					% tuple(a+[kangle,defa]))
	outp.write("\n")
	dlist = []
	for s in seqs:
		for i in range(len(s)-3):
			si = s[i]
			sj = s[i+1]
			sk = s[i+2]
			sl = s[i+3]
			if [sl,sk,sj,si] in dlist:
				sys.stderr.write("########################################\n")
				sys.stderr.write("WARNING!!!!!\n")
				sys.stderr.write("both dihedrals:\n")
				sys.stderr.write("%4s %4s %4s %4s \n"%(sl,sk,sj,si))
				sys.stderr.write("%4s %4s %4s %4s \n"%(si,sj,sk,sl))
				sys.stderr.write("are present!\n")
				sys.stderr.write("########################################\n")
			if [si,sj,sk,sl] not in dlist:
				dlist.append([si,sj,sk,sl])
	outp.write("DIHEDRAL\n")
	dihe_parm = read_karanicolas_dihe()
	for d in dlist:
		ai = d[1]
		aj = d[2]
		tmpp =  dihe_parm[ai][aj]
		for t in tmpp:
			k,mult,phi0 = t
			k*=dihe_fac
			outp.write("%-6s %-6s %-6s %-6s %12.6f %2i %12.6f\n" % \
				tuple(d+[k,mult,phi0]))
	#for i in range(naa):
	#	ai = aalist[i]
	#	for j in range(naa):
	#		aj = aalist[j]
	#		tmpp =  dihe_parm[ai][aj]
	#		for t in tmpp:
	#			k,mult,phi0 = t
	#			k*=dihe_fac
	#			outp.write("%-6s %-6s %-6s %-6s %12.6f %2i %12.6f\n" % \
	#				('X',ai,aj,'X',k,mult,phi0))
	#for j in range(naa):
	#	aj = aalist[j]
	#	for i in range(naa):
	#		ai = aalist[i]
	#		for k in range(i+1):
	#			ak = aalist[k]
	#			outp.write("%-6s %-6s %-6s %12.6f %12.6f\n" \
	#					% (ai,aj,ak,kangle,defa))



def write_defnb_list(outp,defr,eps_res):
	outp.write("\n")
	#outp.write("""NONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -
	# note switch function!!!                   ******
	outp.write("""NONBONDED NBXMOD 3 ATOM CDIEL SWITCH VATOM VDISTANCE VSWITCH -
  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n""")
	eij_rep = -1.5e-3 * eps_res / 21.5 	# fudge
	for aa in aalist:
		outp.write("%4s  0.0  %12.6f %12.6f\n" % (aa,eij_rep,defr/2.0))
	outp.write('\n')

def write_nbfix(outp,defr,pairpot):
	outp.write("\nNBFIX\n")
	mjmap, ave_mj = read_miyazawa_jernigan()
	naa=len(aalist)
	if pairpot == 'mj':
		#for k in mjmap.keys():
		#	for j in mjmap[k].keys():
		for i in range(naa):
			ai = aalist[i]
			for j in range(i+1):
				aj = aalist[j]
				eij = mjmap[k][j]/ave_mj
				outp.write("%4s %4s %12.6e %12.6f\n" \
					% (ai,aj,eij,defr))
	elif pairpot == 'kh':
		LAMBDA = 0.159
		E0 = -2.27
		for i in range(naa):
			ai = aalist[i]
			sigi = res2khsig[ai]
			for j in range(i+1):
				aj = aalist[j]
				sigj = res2khsig[aj]
				sigij = (sigi+sigj)/2.
				epsij = LAMBDA*(mjmap[ai][aj]-E0)*0.6 # to get kcal/mol from kT
				outp.write("%4s %4s %12.6f %12.6f\n"%(ai,aj,epsij,sigij))

	elif pairpot == 'mjbln':
		mjbln = {}
		blnk = ['B','L','N']
		for k in blnk:
			mjbln[k] = {}
			for j in blnk:
				mjbln[k][j]  = [0.,0]
		for k in mjmap.keys():
			bln_k = blnmap[k]
			for j in mjmap[k].keys():
				bln_j = blnmap[j]
				mjbln[bln_k][bln_j][0] += mjmap[k][j]/ave_mj
				mjbln[bln_k][bln_j][1] += 1
		for kk in range(3):
			k = blnk[kk]
			for jj in range(kk+1):
				j = blnk[jj]
				gnum = mjbln[k][j][0]+mjbln[j][k][0]
				gden = mjbln[k][j][1]+mjbln[j][k][1]
				ave = gnum/float(gden)
				mjbln[k][j] = ave
				mjbln[j][k] = ave
		#for k in mjmap.keys():
		#	bln_k = blnmap[k]
		#	for j in mjmap[k].keys():
		#		bln_j = blnmap[j]
		for i in range(naa):
			ai = aalist[i]
			bln_ai = blnmap[ai]
			for j in range(i+1):
				aj = aalist[j]
				bln_aj = blnmap[aj]
				eij = mjbln[bln_ai][bln_aj] 
				outp.write("%4s %4s %12.6e %12.6f\n" \
					% (ai,aj,eij,defr))

	elif pairpot == 'bln':
		bln_pot = { ('B','B'): -2.0,  ('B','L'): -0.1, ('B','N'): -0.4, \
				('L','B'): -0.1,  ('L','L'): -0.1, ('L','N'): -0.4, \
				('N','B'): -0.4,  ('N','L'): -0.4, ('N','N'): -0.4 }
		#blnk = ['B','L','N']
		for i in range(naa):
			ai = aalist[i]
			bln_ai = blnmap[ai]
			for j in range(i+1):
				aj = aalist[j]
				bln_aj = blnmap[aj]
				eij = bln_pot[(bln_k,bln_j)] 
				outp.write("%4s %4s %12.6e %12.6f\n" \
					% (ai,aj,eij,defr))

def write_parms(file,defb,defa,defr,pairpot,genangle,residue_sequences):
	"""top level parameter writing function"""
	outp = open(file,"w")
	write_parm_header(outp,pairpot)
	eps_res = 1.89
	write_bonded(outp,defb,defa,eps_res,genangle,residue_sequences)	# common bit
	write_defnb_list(outp,defr,eps_res)
	write_nbfix(outp,defr,pairpot)
	outp.write("\nEND\n\n")
	outp.close()

def calc_dir(a,b):
	dx = b[0]-a[0]
	dy = b[1]-a[1]
	dz = b[2]-a[2]
	dr = math.sqrt(dx**2+dy**2+dz**2)
	return (dx/dr,dy/dr,dz/dr)

def sc_mult(bondlen,director):
	return (bondlen*director[0],bondlen*director[1],bondlen*director[2])

def rotate(angle,vector):
	x = vector[0]
	y = vector[1]
	cosa = math.cos(angle)
	sina = math.sin(angle)
	nu_x = cosa*x+sina*y
	nu_y = -sina*x+cosa*y
	return (nu_x,nu_y,0.0)

def addv(a,b):
	return (a[0]+b[0],a[1]+b[1],a[2]+b[2])

def make_sequ_coords(sequences,defb,defa):
	""" create linear coordinate data as for pdb file reading, but from sequence"""
	alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	crd_dat = { }
	xyz = {  }
#
	sign = 1.0
	for s in range(len(sequences)):
		chain = alpha[s]
		crd_dat[chain] = {}
		xyz[chain] = []
		sequence = sequences[s]
		for i in range(len(sequence)):
			residue = i+1
			#if sequence[i] not in aamap.keys():
			#	print "Don't know about amino acid type \"%s\"" % sequence[i]
			#	sys.exit(1)
			name = sequence[i]
			crd_dat[chain][residue] = {}
			crd_dat[chain][residue]['name'] = name
			if i==0:
				x,y,z = 0.0,0.0,s*15.
			elif i==1:
				x,y,z = defb,0.0,s*15.
			else:
				bondlen = defb
				tmpv = sc_mult(bondlen,director)
				tmpv = rotate(sign*(math.pi-defa*math.pi/180.0),tmpv)
				x,y,z = addv(xyz[chain][-1][1:],tmpv)
				sign *= -1.0
			xyz[chain].append((residue,x,y,z))
			crd_dat[chain][residue]["CA"] = (x,y,z)
			if i>0:
				director = calc_dir(xyz[chain][-2][1:],xyz[chain][-1][1:])
	return crd_dat,xyz
	#return xyz

def write_CA_pdb(filename, pdbdat, crd_dat, chains):
	"""write CA-only coordinates"""
	outp = open(filename,"w")
	at = 0
	chains.sort()
	for chain in chains:
		pdbc = pdbdat[chain]
		for p in pdbc:
			res = p[0]
			name = crd_dat[chain][res]['name']
			at+=1
			if res < 1000:
				outp.write("ATOM  %5i  CA  %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
					% (at,name,res,p[1],p[2],p[3],1.00,0.00,chain))
			else:
				outp.write("ATOM  %5i  CA  %3s   %4i   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n" \
					% (at,name,res,p[1],p[2],p[3],1.00,0.00,chain))
	outp.write("END\n")
	outp.close()


# write topology
######################################################################
def write_aa_top(filename):
	top_out = open(filename,"w")

	top_out.write("* Topology for generic ff\n*\n  20   1\n\n");

	aa = massmap.keys()
	naa = len(aa)

	for i in range(naa):
		a = aa[i]
		amass = massmap[a]
		top_out.write("MASS %3i %3s  %8.3f\n"%(i+1,a,amass))

	top_out.write("\nDECL +CA\n\nAUTOGENERATE ANGLES DIHEDRAL\n\n");

	for i in range(naa):
		a = aa[i]
		acharge = res2charge[a]
		top_out.write("RESI   %3s  %5.3f\n"%(a,acharge))
		top_out.write("GROUP\n")
		top_out.write("ATOM   CA   %s  %5.3f\n"%(a,acharge))
		top_out.write("BOND CA +CA\n\n")

	top_out.write("END\n\n")
	top_out.close()


