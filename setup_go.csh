#!/bin/csh -f


#./go_builder.py -t 350 -k ga -n 1 -g 0.0 charmm_setup/ga_polarh_wt.pdb
#./go_builder_sc.py -t 350 -k ga -n 1 -g 0.0 charmm_setup/ga_polarh_wt.pdb
#./go_builder.py -t 350 -k ga_geom -n 1 -g 0.0 -G charmm_setup/ga_polarh_wt.pdb
#./go_builder_sc.py -t 350 -k ga_geom -n 1 -g 0.0 -G charmm_setup/ga_polarh_wt.pdb

./go_builder.py -t 350 -k ga_orland -n 1 -g 0.0 -O charmm_setup/ga_polarh_wt.pdb

#Usage:
#	go_builder.py [-h] [--help] [-t Tf] [--tf=Tf] [-k key] [--key=key]
#	[-g gamma] [-n nrep] [-d nongo_dist] [--nrepeat=nrep] file.pdb 
#
#	in which:
#		Tf is the desired folding temperature (default = 350 K)
#		key is a handy key for naming the output files (default = junk)
#		nrep is the number of concatenated repeats to make (default = 1)
#		gamma is the fraction of the native contact energy assigned to
#			non-native contacts (default = 0.0)
#	
#	Only the pdb file is a required argument. It must have all backbone heavy
#	atoms and amide hydrogens defined, though side-chain heavy atoms are
#	used too.





