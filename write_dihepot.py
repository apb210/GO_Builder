#!/usr/bin/env python

import sys, math

dihe_parm = {}

for line in open("karanicolas_dihe_parm.dat").readlines():
	ls = line.split()
	a1, a2, k, m, d = ls[0],ls[1],float(ls[2]),int(ls[3]),float(ls[4])
	if a1 not in dihe_parm.keys():
		dihe_parm[a1] = {}
	if a2 not in dihe_parm[a1].keys():
		dihe_parm[a1][a2] = []
	dihe_parm[a1][a2].append((k,m,d))

for a1 in dihe_parm.keys():
	for a2 in dihe_parm[a1].keys():
		print a1, a2
		plist = dihe_parm[a1][a2]
		pout = open("dihedral_potentials/%s_%s.dat" % (a1,a2),"w")
		for t in range(-180,181,1):
			tt = float(t)*math.pi/180.0
			E=0
			for p in plist:
				E+=p[0]*(1+math.cos(float(p[1])*tt-p[2]*math.pi/180.0))
			pout.write("%8.3f %8.3f\n" % ( float(t), E ))
		pout.close()

