#!/usr/bin/env python

"""Chop a selection out of a pdb file"""

from PDBFile import *

usage = """\
Usage:
	pdbchop p=xyz.pdb o=out.pdb [ c=chain ] i=initial f=final 
		[ r=[0|1] ] [ m=model ]
where
	xyz.pdb = input pdbfile
	out.pdb = output pdbfile
	chain = (optional) chain ID - else set to ' '
	initial = first residue (in original numbering scheme)
	final   = last residue
	renumber = renumber residues from 1? (default: on)"""

if __name__ == "__main__":
	import sys
	args = sys.argv[1:]
	#if (len(args) < 4):
	#	print usage
	#	sys.exit
	p = ''
	c = ' '
	i = f = -1
	r = 1
	m =0
	for arg in args:
		if arg[0] == 'p':
			p = arg[2:]
		elif arg[0] == 'o':
			o = arg[2:]
		elif arg[0] == 'c':
			c = arg[2:]
		elif arg[0] == 'i':
			i = int(arg[2:])
		elif arg[0] == 'f':
			f = int(arg[2:])
		elif arg[0] == 'r':
			r = int(arg[2:])
		elif arg[0] == 'm':
			m = int(arg[2:])
		else:
			print "Incorrect command-line usage."
			print usage
			sys.exit(1)

	if p=='' or o == '' or i==-1 or f==-1:
		#not properly initialized
		print "Incorrect command-line usage."
		print usage
		sys.exit(1)

	pf = PDBFile(p)
	if pf.nmodel > 0:
		print "PDB file contains multiple models, but none specified"
		print "Using model # 1 ..."
		m = 1
	offset = pf.res_offset(c, i, m)
	if offset < 0:
		print "Could not find residue # ", i
		sys.exit(1)
	length = f - i + 1
	pd = PDBData(pf, c, offset, length, m, r)
	po = PDBOutput(o) 
	po.append(pd, 'A')

		
