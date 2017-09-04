#!/usr/bin/env python

import sys,os,math


def parse_golist(file):
	data = map(lambda x: x.split(), open(file).readlines())
	golist = map(lambda x: (int(x[0]),int(x[1]),float(x[2]),float(x[3])),\
			data)
	return golist

def avenb(golists,cut):
	nstr = len(golists)
	ave_nnb = 0
	nbtmp = {}
	for g in golists:
		ave_nnb += len(g)
		for contact in g:
			i,j,rij,eij = contact
			#int(nbs[0]),int(nbs[1]),\
			#		float(nbs[2]),float(nbs[3])
			if (i,j) not in nbtmp.keys():
				nbtmp[(i,j)] = [ rij, eij, 1 ]
			else:
				nbtmp[(i,j)][0] += rij
				nbtmp[(i,j)][1] += eij
				nbtmp[(i,j)][2] += 1
	ave_nnb = float(ave_nnb)/float(nstr)

	avenb = {}
	nbt_keys = nbtmp.keys()
	nbt_keys.sort()
	for k in nbt_keys:
		freq = nbtmp[k][2]
		if float(freq)/float(nstr) > cut:
			if k not in avenb.keys():
				rij = nbtmp[k][0]/float(freq)
				eij = nbtmp[k][1]/float(freq)
				avenb[k] = (eij,rij)
	return avenb, ave_nnb


golists = map(parse_golist, sys.argv[2:])
cut = float(sys.argv[1])

avenb, ave_nnb = avenb(golists,cut)

for k in avenb.keys():
	i,j = k
	eij, rij = avenb[k]
	sys.stdout.write("%5i %5i %12.6f %12.6f\n" % (i,j,rij,eij))

