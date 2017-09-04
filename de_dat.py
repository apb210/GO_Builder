#!/usr/bin/env python

# make mutation input file for quick_epert 

import  sys

def parse_parm(inp):
	#inp = open(file).read()
	sections = inp.read().split('\n\n')
	header = sections[0].strip()
	return map(lambda x: x.strip().split('\n'), sections)

def mutate_nbfix(nbfix,res,scale):
	mdat = []
	for nb in nbfix:
		nbs = nb.split()
		ii,jj,eij,rij = nbs[0],nbs[1],float(nbs[2]),float(nbs[3])
		i = int(ii[1:])
		j = int(jj[1:])
		if i==res or j==res:
			de_ij = eij*scale
			mdat.append((i,j,de_ij,rij))
	return mdat


mutres = int(sys.argv[1])
mutfrac = float(sys.argv[2])

parmdat = parse_parm(sys.stdin)

mut_dat = mutate_nbfix(parmdat[7][1:],mutres,mutfrac)

for nb in mut_dat:
	sys.stdout.write("%5i %5i %12.6f %12.6f\n" % (nb))

