#!/usr/bin/env python

import sys,os,math

fddir = "fd_data"

# stolen from nab ===================================================
hxht = {}; hxrep = {}
helical_rise = {"arna": 2.81, "aprna": 3.00, "lbdna": 3.38,
		"abdna": 3.38, "sbdna": -3.38, "adna": 2.56 }
helical_twist = {"arna": 32.7, "aprna": 30.0, "lbdna": 36.0,
		"abdna": 36.0, "sbdna": 36.0, "adna": 32.7 }
# stolen from nab ===================================================

mass_map = { "P": 31., "O": 16., "C": 12., "N": 14. }

na_complement = { 'DNA': { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' },
		'RNA':{ 'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' } } 

def parse_fd(fdkey):
	"""parse fibre diffraction data"""
	aa_map = {}
	for l in open("%s/%s.dat"%(fddir,fdkey)).readlines():
		ls = l.split()
		res = ls[4]
		name = ls[0]
		r = float(ls[1])
		phi = float(ls[2])
		z = float(ls[3])
		if res not in aa_map.keys():
			aa_map[res] = []
		aa_map[res].append((name,r,phi,z))
	return aa_map

def cart2cyl(x,y):
	r = math.sqrt(x**2+y**2)
	phi = 180./math.pi*math.asin(y/r)
	if x < 0:
		phi = 180.-phi
	return r,phi

# once-off function to generate the "canonical" CG coordinates
def cg_reduce(aa_map):
	cg_map = {}
	for res in aa_map.keys():
		cg_map[res] = []
		px,py,pz,mp = 0., 0., 0., 0.
		rx,ry,rz,mr = 0., 0., 0., 0.
		bx,by,bz,mb = 0., 0., 0., 0.
		for atom in aa_map[res]:
			name,r,phi,z = atom
			x = r * math.cos(phi*math.pi/180.)
			y = r * math.sin(phi*math.pi/180.)
			mass = mass_map[name[0]]
			if name[-1] == "P":
				px += mass*x
				py += mass*y
				pz += mass*z
				mp += mass
			elif name[-1] == "'":
				rx += mass*x
				ry += mass*y
				rz += mass*z
				mr += mass
			else:
				bx += mass*x
				by += mass*y
				bz += mass*z
				mb += mass
		px /= mp; py /= mp; pz /= mp;
		rx /= mr; ry /= mr; rz /= mr;
		bx /= mb; by /= mb; bz /= mb;
		pr,pphi = cart2cyl(px,py)
		rr,rphi = cart2cyl(rx,ry)
		br,bphi = cart2cyl(bx,by)
		cg_map[res].append(('P',pr,pphi,pz))
		cg_map[res].append(('R',rr,rphi,rz))
		cg_map[res].append(('B',br,bphi,bz))
	return cg_map

def cg_base_shift(cg_map,off):
	"""move CG bases closer to helical axis"""
	new_cg_map = {}
	for res in cg_map.keys():
		new_cg_map[res] = []
		for cg in cg_map[res]:
			name,r,phi,z = cg
			if name == 'B':
				new_cg_map[res].append(('B',off,phi,z))
			else:
				new_cg_map[res].append(cg)
	return new_cg_map

def cg_base_scale(cg_map,sfac):
	"""move CG bases closer to helical axis"""
	new_cg_map = {}
	for res in cg_map.keys():
		new_cg_map[res] = []
		for cg in cg_map[res]:
			name,r,phi,z = cg
			off = r*sfac
			if name == 'B':
				new_cg_map[res].append(('B',off,phi,z))
			else:
				new_cg_map[res].append(cg)
	return new_cg_map

def write_fd(pmap,key):
	outp = open("%s/%s.dat"%(fddir,key),"w")
	for k in pmap.keys():
		for d in pmap[k]:
			name, r, phi, z = d
			outp.write("%5s %8.3f %8.3f %8.3f %4s\n"\
					%(name,r,phi,z,k))

def complement(seq,na_type = 'DNA'):
	comp = ""
	for s in seq:
		comp = na_complement[na_type][s] + comp
	return comp


def fd_build_duplex(seq,fd_type = 'abdna'):
	if fd_type.find('dna')>=0:
		na_type = "DNA"
	else:
		na_type = "RNA"
	if fd_type.find('_') > 0:
		fd_key = fd_type[:fd_type.find('_')]
	else:
		fd_key = fd_type
	canonical_rphiz = parse_fd(fd_type)
	# build 5'-3' strand
	strand53_crd = []
	chain = "A"
	resnum = 0
	nres = len(seq)
	for s in range(nres):
		resnum += 1
		total_rise = s*helical_rise[fd_key]
		total_twist = s*helical_twist[fd_key]
		for atom in canonical_rphiz[seq[s]]:
			name,r,phi,z = atom
			phi += total_twist
			z += total_rise
			x = r*math.cos(phi*math.pi/180.)
			y = r*math.sin(phi*math.pi/180.)
			strand53_crd += [(name,seq[s],chain,resnum,x,y,z)]
	# build complementary 3'-5' strand
	strand35_crd = []
	chain = "B"
	cseq = complement(seq,na_type)
	print seq
	print cseq
	for s in range(nres):
		resnum += 1
		total_rise = (nres-s-1)*helical_rise[fd_key]
		total_twist = (nres-s-1)*helical_twist[fd_key]
		for atom in canonical_rphiz[cseq[s]]:
			name,r,phi,z = atom
			phi = total_twist-phi
			z = total_rise - z
			x = r*math.cos(phi*math.pi/180.)
			y = r*math.sin(phi*math.pi/180.)
			strand35_crd += [(name,cseq[s],chain,resnum,x,y,z)]
	all =  strand53_crd + strand35_crd
	return all

def fd_build_duplex_nonmatch(seq1,seq2,fd_type = 'abdna'):
	if fd_type.find('dna')>=0:
		na_type = "DNA"
	else:
		na_type = "RNA"
	if fd_type.find('_') > 0:
		fd_key = fd_type[:fd_type.find('_')]
	else:
		fd_key = fd_type
	canonical_rphiz = parse_fd(fd_type)
	# build 5'-3' strand
	strand53_crd = []
	chain = "A"
	resnum = 0
	nres1 = len(seq1)
	for s in range(nres1):
		resnum += 1
		total_rise = s*helical_rise[fd_key]
		total_twist = s*helical_twist[fd_key]
		for atom in canonical_rphiz[seq1[s]]:
			name,r,phi,z = atom
			phi += total_twist
			z += total_rise
			x = r*math.cos(phi*math.pi/180.)
			y = r*math.sin(phi*math.pi/180.)
			strand53_crd += [(name,seq1[s],chain,resnum,x,y,z)]
	# build complementary 3'-5' strand
	strand35_crd = []
	chain = "B"
	#cseq = complement(seq,na_type)
	#print seq
	#print cseq
	nres2 = len(seq2)
	for s in range(nres2):
		resnum += 1
		total_rise = (nres1-s-1)*helical_rise[fd_key]
		total_twist = (nres1-s-1)*helical_twist[fd_key]
		for atom in canonical_rphiz[seq2[s]]:
			name,r,phi,z = atom
			phi = total_twist-phi
			z = total_rise - z
			x = r*math.cos(phi*math.pi/180.)
			y = r*math.sin(phi*math.pi/180.)
			strand35_crd += [(name,seq2[s],chain,resnum,x,y,z)]
	all =  strand53_crd + strand35_crd
	return all

def calc_bonded(coords):
	"""compute bonded (bonds,angles,torsions) terms from structure"""
	print "blah"

def calc_nonbonded(coords):
	"""compute non-bonded distances from structure"""
	print "blah"

def write_pdb(coords,filename):
	outp = open(filename, "w")
	outp.write("REMARK  File written by na_lib\n")
	for i in range(len(coords)):
		c = coords[i]
		#print c
		format = "%-6s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n"
		outp.write(format%("ATOM  ",i+1,c[0],c[1],c[2],c[3],c[4],
			c[5],c[6],1.0,0.0,""))







