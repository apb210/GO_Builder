#!/bin/csh

./go_builder.py -t 350 -k lys -D 6-128:30-116:65-81:77-95 lysozyme.pdb

charmm < test_lys.inp > test_lys.out
