#!/usr/bin/env python

#./enm.py xyz.pdb springs.dat output_directory
#
# xyz.pdb - minimized pdb
# springs.dat - list of atom pairs and spring constants (kcal/mol)
# output_directory - frequencies and pdb's for lowest 50 modes 
#		will go here

import sys,os,math,numarray,numarray.linear_algebra

def read_dVdr(file):
	inp = open(file,"r")
	raw = map(lambda x: x.split(), inp.readlines())
	natom = 0
	dVdr = []
	for r in raw:
		i,j,hij = int(r[0]),int(r[1]),float(r[2])
		dVdr.append((i-1,j-1,hij))
		natom = max(i,j,natom)
	#H = numarray.zeros([natom,natom],numarray.Float64)
	#for d in dVdr:
	#	i,j,hij = d 
	#	H[i-1][j-1] = hij
	#	H[j-1][i-1] = hij
	#return H
	return dVdr, natom

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
	H = numarray.zeros([3*natom,3*natom],numarray.Float64)
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
	
def normal_disp(X,Y,Z,evecs,amplitudes,T,usemodes=-1):
	XX, YY, ZZ = [], [], []
	natom = len(X)
	for i in range(natom):
		XX.append(X[i])
		YY.append(Y[i])
		ZZ.append(Z[i])
	if usemodes < 0:
		nmode = len(amplitudes)
	else:
		nmode = usemodes
	for m in range(nmode):
		ampl = amplitudes[m]
		for i in range(natom):
			XX[i] += evecs[m][3*i]*ampl
			YY[i] += evecs[m][3*i+1]*ampl
			ZZ[i] += evecs[m][3*i+1]*ampl
	return XX, YY, ZZ


