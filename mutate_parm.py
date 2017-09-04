#!/usr/bin/env python

# 'mutate' interactions in parameter file by altering 
# interactions with a residue

import  sys

def parse_parm(inp):
	#inp = open(file).read()
	sections = inp.read().split('\n\n')
	header = sections[0].strip()
	return map(lambda x: x.strip().split('\n'), sections)

def mutate_nbfix(nbfix,res,scale):
	new_nbfix = []
	for nb in nbfix:
		nbs = nb.split()
		ii,jj,eij,rij = nbs[0],nbs[1],float(nbs[2]),float(nbs[3])
		i = int(ii[1:])
		j = int(jj[1:])
		if i==res or j==res:
			new_eij = eij*(1.0-scale)
		else:
			new_eij = eij
		new_nbfix.append((ii,jj,new_eij,rij))
	return new_nbfix


mutres = int(sys.argv[1])
mutfrac = float(sys.argv[2])

parmdat = parse_parm(sys.stdin)

new_nbfix = mutate_nbfix(parmdat[7][1:],mutres,mutfrac)

for block in parmdat[:7]:
	for b in block:
		sys.stdout.write("%s\n" % (b))
	sys.stdout.write("\n")

sys.stdout.write("\nNBFIX\n")
for nb in new_nbfix:
	sys.stdout.write("%-5s %-5s %12.6f %12.6f\n" % (nb))

sys.stdout.write("\nEND\n\n")

