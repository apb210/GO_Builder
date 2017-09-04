#!/usr/bin/env python

import sys,os,math

# make stupid ENM

def read_pdb(pdb_file):
	pdat = filter(lambda x: x.find('ATOM')==0, open(pdb_file).readlines())
	txtdat = map(lambda x: x[:30], pdat)
	X = map(lambda x: float(x[30:38]), pdat)
	Y = map(lambda x: float(x[38:46]), pdat)
	Z = map(lambda x: float(x[46:54]), pdat)
	return txtdat, X, Y, Z
	

pdbfile = sys.argv[1]
cut = float(sys.argv[2])
k = float(sys.argv[3]) #kcal/mol/A**2

txt, X,Y,Z = read_pdb(pdbfile)

natom = len(X)

for i in range(natom):
	for j in range(i+1,natom):
		dx = X[i] - X[j]
		dy = Y[i] - Y[j]
		dz = Z[i] - Z[j]
		dr = math.sqrt(dx**2+dy**2+dz**2)
		if dr<cut:
			sys.stdout.write("%5i %5i %8.3f\n"%(i+1,j+1,k))



