#!/usr/bin/env python 

import sys,os,enm,numarray,numarray.linear_algebra

#======================================================================
pdbfile = sys.argv[1]
dVdrfile = sys.argv[2]
odir = sys.argv[3]
#======================================================================

if not os.path.exists(odir):
	os.mkdir(odir)

txtdat, X, Y, Z = enm.read_pdb(pdbfile)
dVdr, natom = enm.read_dVdr(dVdrfile)
H = enm.calc_hessian(dVdr,X,Y,Z)

ndump = 50
#print dVdr
#print H[1]
#sys.exit(0)
evals, evecs = numarray.linear_algebra.Heigenvectors(H)

evout = open("%s/evals.dat"%(odir),"w")
for e in evals:
	evout.write("%12.5e\n"%(e))
evout.close()

for i in range(ndump):
	Xtmp = 10*numarray.array(map(lambda x: evecs[i][x], range(0,3*natom,3)),numarray.Float64)
	Ytmp = 10*numarray.array(map(lambda x: evecs[i][x], range(1,3*natom,3)),numarray.Float64)
	Ztmp = 10*numarray.array(map(lambda x: evecs[i][x], range(2,3*natom,3)),numarray.Float64)
	Xtmp += X
	Ytmp += Y
	Ztmp += Z
	enm.write_pdb("%s/evec_%i.pdb"%(odir,i),txtdat, Xtmp, Ytmp, Ztmp)

