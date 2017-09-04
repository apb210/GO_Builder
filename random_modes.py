#!/usr/bin/env python 

# write superposition of random modes @ temperature T
import sys,os,enm,numarray,numarray.linear_algebra,numarray.random_array,math

#======================================================================
pdbfile = sys.argv[1]
dVdrfile = sys.argv[2]	# spring constants assumed to be kcal/mol/(A**2)
odir = sys.argv[3]
nmode = int(sys.argv[4])
nstruct = int(sys.argv[5])
T = float(sys.argv[6])
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

dim = len(evals)-6
means = numarray.zeros(dim,numarray.Float64)
kT = 8.31*T/4184.0
sdevs = numarray.array(map(lambda x: math.sqrt(kT/x), evals[6:]),numarray.Float64)

for s in range(nstruct):
	amplitudes = numarray.random_array.normal(means,sdevs)
	Xtmp,Ytmp,Ztmp = enm.normal_disp(X,Y,Z,evecs[6:],amplitudes,nmode)
	enm.write_pdb("%s/%i.pdb"%(odir,s),txtdat, Xtmp, Ytmp, Ztmp)
	

