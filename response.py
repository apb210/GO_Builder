#!/usr/bin/env python

#./enm.py xyz.pdb springs.dat output_directory
#
# xyz.pdb - minimized pdb
# springs.dat - list of atom pairs and spring constants (kcal/mol)
# output_directory - frequencies and pdb's for lowest 50 modes 
#		will go here

import sys,os,math,numarray,numarray.linear_algebra

def vecout(v,name):
	outp = open(name,"w")
	map(lambda x: outp.write("%12.5e\n"%(x)), v)
	outp.close()

def read_dVdr(file):
	inp = open(file,"r")
	raw = map(lambda x: x.split(), inp.readlines())
	natom = 0
	dVdr = []
	for r in raw:
		i,j,hij = int(r[0]),int(r[1]),float(r[2])
		dVdr.append((i-1,j-1,hij))
		natom = max(i,j,natom)
	return dVdr, natom

def read_cons(file):
	inp = open(file,"r")
	raw = map(lambda x: x.split(), inp.readlines())
         
	dVdr = []
	for r in raw:
		i,j,hij,rij0 = int(r[0]),int(r[1]),float(r[2]),float(r[3])
		dVdr.append((i-1,j-1,hij,rij0))
	return dVdr

def read_pdb(pdb_file):
	pdat = filter(lambda x: x.find('ATOM')==0, open(pdb_file).readlines())
	txtdat = map(lambda x: x[:30], pdat)
	X = map(lambda x: float(x[30:38]), pdat)
	Y = map(lambda x: float(x[38:46]), pdat)
	Z = map(lambda x: float(x[46:54]), pdat)
	return txtdat, X, Y, Z
	
def write_pdb(file,txtdat, X, Y, Z):
	natom = len(X)
	outp = open(file,"w")
	for i in range(natom):
		outp.write("%30s%8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n" \
				% ( txtdat[i], X[i],Y[i],Z[i],1.0,0.0,"PROT" ))
	outp.write("END\n")
	outp.close()

def calc_hessian(dVdr,X,Y,Z):
	natom = len(X)
	H = numarray.zeros([3*natom,3*natom],numarray.Float32)
	for d in dVdr:
		i,j,hij = d 
		dr = ( X[i]-X[j], Y[i]-Y[j], Z[i]-Z[j] )
		r2 = dr[0]**2+dr[1]**2+dr[2]**2
		for p in range(3):
			drp = dr[p]
			for q in range(3):	# select dxij, dyij etc.
				drq = dr[q]
				dVdpdq = hij*drp*drq/r2
				for ij in (i,j):
					pidx = 3*ij+p
					if ij == i:
						sij = 1
					else:
						sij = -1
					for kl in (i,j):
						qidx = 3*kl+q
						if kl == i:
							skl = 1
						else:
							skl = -1
						H[pidx][qidx] += dVdpdq*sij*skl
	return H
	
def calc_f(dVdr,X,Y,Z):
	natom = len(X)
	f = numarray.zeros([3*natom],numarray.Float32)
	for d in dVdr:
		i,j,hij,rij0 = d 
		drij = ( X[i]-X[j], Y[i]-Y[j], Z[i]-Z[j] )
		rij = math.sqrt(drij[0]**2+drij[1]**2+drij[2]**2)
		for p in range(3):	# single loop over x,y,z
			drp = drij[p]
			#dVdpdq = hij*drp*drq/r2
			dVdp = hij*drp/rij*(rij-rij0)
			for ij in (i,j):
				pidx = 3*ij+p
				if ij == i:
					sij = -1
				else:
					sij = +1
				f[pidx] += dVdp*sij
	return f

	
def response(H,f,nmode,T):
	R = 0.0019861376	# kcal/mol
	beta = 1./(R*T)
	evals, evecs = numarray.linear_algebra.Heigenvectors(H)
	dim = len(evecs[0])
	dr = numarray.zeros(dim,numarray.Float32)
	for m in range(6,6+nmode):
		dot = numarray.dot(evecs[m],f)
		print dot,evals[m]
		# TODO!!! CHECK THIS AMPLITUDE FACTOR
		dr += math.sqrt(beta/evals[m])*dot*evecs[m]
	return dr

#======================================================================
pdbfile = sys.argv[1]
dVdrfile = sys.argv[2]
consfile = sys.argv[3]
odir = sys.argv[4]
#======================================================================

if not os.path.exists(odir):
	os.mkdir(odir)

txtdat, X, Y, Z = read_pdb(pdbfile)
dVdr, natom = read_dVdr(dVdrfile)
cons = read_cons(consfile)
H = calc_hessian(dVdr,X,Y,Z)
f = calc_f(cons,X,Y,Z)
vecout(f,"f.dat")

#dr = response(H,f,50,300.0)
dr = response(H,f,500,300.0)
#notbonds = 3*natom-(natom-1)
#dr = response(H,f,notbonds,300.0)
vecout(dr,"dr.dat")
evals, evecs = numarray.linear_algebra.Heigenvectors(H)

ndump = 30
for i in range(ndump):
	Xtmp = numarray.array(map(lambda x: evecs[i][x], range(0,3*natom,3)),numarray.Float32)
	Ytmp = numarray.array(map(lambda x: evecs[i][x], range(1,3*natom,3)),numarray.Float32)
	Ztmp = numarray.array(map(lambda x: evecs[i][x], range(2,3*natom,3)),numarray.Float32)
	Xtmp += X
	Ytmp += Y
	Ztmp += Z
	opdb = "%s/%i.pdb" % (odir,i)
	write_pdb(opdb,txtdat, Xtmp, Ytmp, Ztmp)

Xtmp = numarray.array(map(lambda x: dr[x], range(0,3*natom,3)),numarray.Float32)
Ytmp = numarray.array(map(lambda x: dr[x], range(1,3*natom,3)),numarray.Float32)
Ztmp = numarray.array(map(lambda x: dr[x], range(2,3*natom,3)),numarray.Float32)
Xtmp += X
Ytmp += Y
Ztmp += Z

opdb = "%s/response.pdb" % odir

write_pdb(opdb,txtdat, Xtmp, Ytmp, Ztmp)

