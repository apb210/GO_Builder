#!/usr/bin/env python

import  sys

def parse_parm(file):
	inp = open(file).read()
	sections = inp.split('\n\n')
	header = sections[0].strip()
	#readparam = sections[1].strip()
	#bond = sections[2].strip()
	#angle = sections[3].strip()
	#dihedral = sections[4].strip()
	#nonbongen = sections[5].strip()
	#nonbondef = sections[6].strip()
	#nonbonatt = sections[7].strip()
	return map(lambda x: x.strip().split('\n'), sections)

def ave_bond(bond_dat):
	nstr = len(bond_dat)
	nbon = len(bond_dat[0])
	bondsum = []
	for i in range(nbon):
		bondsum.append(0.0)
	for bd in bond_dat:
		for i in range(nbon):
			bond = bd[i]
			bondsum[i] += float(bond.split()[3])
	ave_bond = []
	for i in range(nbon):
		bs = bond_dat[0][i].split()
		p,q,e = bs[0],bs[1],float(bs[2])
		blen = bondsum[i]/float(nstr)
		ave_bond.append((p,q,e,blen))
	return ave_bond

def ave_angle(angle_dat):
	nstr = len(angle_dat)
	nang = len(angle_dat[0])
	angsum = []
	stata = 0
	for i in range(nang):
		angsum.append(0.0)
	for ad in angle_dat:
		for i in range(nang):
			ang = ad[i].split()
			anglen = len(ang)
			if anglen > 5:
				stata = 1
				break
			angsum[i] += float(ang[4])
	ave_ang = []
	if stata == 1:
		for i in range(nang):
			as = angle_dat[0][i].split()
			p,q,r = as[0],as[1],as[2]
			ave_ang.append((p,q,r,106.4,91.7,26.3,130.0,0.1,4.3))
	else:
		for i in range(nang):
			as = angle_dat[0][i].split()
			p,q,r,e = as[0],as[1],as[2],float(as[3])
			ang = angsum[i]/float(nstr)
			ave_ang.append((p,q,r,e,ang))
	return ave_ang

def minrep(defnbdat):
	nstr = len(defnbdat)
	nnb = len(defnbdat[0])
	minrep = []
	for i in range(nnb):
		tmp = defnbdat[0][i].split()
		p, a, b = tmp[0], float(tmp[1]), float(tmp[2])
		minrep.append([p,a,b,1000.0])
	for dnb in defnbdat:
		for i in range(nnb):
			drep = float(dnb[i].split()[3])
			if drep < minrep[i][3]:
				minrep[i][3] = drep
	return minrep

def avenb(nbfix,cut):
	nstr = len(nbfix)
	ave_nnb = 0
	nbtmp = {}
	nntmp = {}
	for nbf in nbfix:
		ave_nnb += len(nbf)
		for nb in nbf:
			nbs = nb.split()
			ii,jj,eij,rij = nbs[0],nbs[1],\
					float(nbs[2]),float(nbs[3])
			i = int(ii[1:])
			j = int(jj[1:])
			if rij == 5.5:
				if (i,j) not in nntmp.keys():
					nntmp[(i,j)] = ( eij, rij )
				continue	# non-native...ignore
			if (i,j) not in nbtmp.keys():
				nbtmp[(i,j)] = [ eij, rij, 1 ]
			else:
				nbtmp[(i,j)][0] += eij
				nbtmp[(i,j)][1] += rij
				nbtmp[(i,j)][2] += 1
	ave_nnb = float(ave_nnb)/float(nstr)

	avenb = {}
	nbt_keys = nbtmp.keys()
	nbt_keys.sort()
	for k in nbt_keys:
		freq = nbtmp[k][2]
		if float(freq)/float(nstr) > cut:
			if k not in avenb.keys():
				eij = nbtmp[k][0]/float(freq)
				rij = nbtmp[k][1]/float(freq)
				avenb[k] = (eij,rij)
	return avenb, ave_nnb, nntmp


cfrac = float(sys.argv[1])
parmfiles = sys.argv[2:]

parmdat = []

for f in parmfiles:
	parmdat.append(parse_parm(f))

bond_dat = ave_bond(map(lambda x: x[2][1:], parmdat))
angle_dat = ave_angle(map(lambda x: x[3][1:], parmdat))
dihe_dat = parmdat[0][4]
nbdef = minrep(map(lambda x: x[6], parmdat))
nbatt,avennb,nonnb = avenb(map(lambda x: x[7][1:], parmdat),cfrac)

for l in parmdat[0][0]:
	sys.stdout.write("%s\n" % (l))
sys.stdout.write("\n")
for l in parmdat[0][1]:
	sys.stdout.write("%s\n" % (l))
sys.stdout.write("\nBOND\n")
for b in bond_dat:
	sys.stdout.write("%-5s %-5s %12.5f %12.5f\n" % ( b ))
sys.stdout.write("\nANGLE\n")
if len(angle_dat[0])==5:	# regular angles
	for a in angle_dat:
		sys.stdout.write("%-5s %-5s %-5s %12.5f %12.5f\n" % ( a ))
else:
	for a in angle_dat:
		sys.stdout.write("%-5s %-5s %-5s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n" % ( a ))

sys.stdout.write("\n")
for d in dihe_dat:
	sys.stdout.write("%s\n" % (d))

sys.stdout.write("\n")
for l in parmdat[0][5]:
	sys.stdout.write("%s\n" % (l))

#print nbdef

sys.stdout.write("\n")
for nb in nbdef:
	sys.stdout.write("%-5s %5.3f %12.6f %12.6f\n" % tuple( nb ))

sys.stdout.write("\nNBFIX\n")
nres = len(bond_dat)+1
for i in range(1,nres+1):
	for j in range(i+1,nres+1):
		if (i,j) in nbatt.keys():
			eij,rij = nbatt[(i,j)]
			sys.stdout.write("G%i G%i %12.6e %12.6f\n" % \
						(i,j,eij,rij))
		elif (i,j) in nonnb.keys():
			eij,rij = nonnb[(i,j)]
			sys.stdout.write("G%i G%i %12.6e %12.6f\n" % \
						(i,j,eij,rij))
#print len(nbatt)
#print avennb
sys.stdout.write("\nEND\n\n")

