#!/usr/bin/env python

import sys,os,go_lib

metac, xyz = go_lib.make_coords(sys.stdin.readlines())
chains = metac.keys()
#print metac
#c = chains[0]
#for i in range(1,10):
#	print i, metac[c][i]["name"]
#sys.exit(0)


go_lib.write_CA_SC_pdb("stdout", xyz, metac, chains,1)
