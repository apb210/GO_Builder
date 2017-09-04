#!/usr/bin/env python

import na_lib,sys,os

#print na_lib.complement("ATGGC")

# make reduced cylindrical coordinats

#for s in os.listdir("fd_data"):
#	ss = s[:-4]
#	if ss in "CVS" or s[-1] == 'g':
#		continue
#	aamap = na_lib.parse_fd(ss)
#	cgmap = na_lib.cg_reduce(aamap)
#	na_lib.write_fd(cgmap,ss+"_cg")
#



testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna")
na_lib.write_pdb(testbuild,"test_atat_b.pdb")
##
testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna_cg")
na_lib.write_pdb(testbuild,"test_atat_b_cg.pdb")

testbuild = na_lib.fd_build_duplex("AUAUAUAUAUAU","arna")
na_lib.write_pdb(testbuild,"test_auau_rna.pdb")

testbuild = na_lib.fd_build_duplex_nonmatch("GGCGAAGUCGA","AAGAUGGCGCC","arna")
na_lib.write_pdb(testbuild,"test_f9l_rna.pdb")

sys.exit(0)

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","adna_cg")
na_lib.write_pdb(testbuild,"test_cgcg_a_cg.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","adna_cgr")
na_lib.write_pdb(testbuild,"test_cgcg_a_cgr.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","adna")
na_lib.write_pdb(testbuild,"test_cgcg_a.pdb")
sys.exit(0)
#
testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna_cg_s0.6")
na_lib.write_pdb(testbuild,"test_atat_b_cg_s0.6.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","abdna_cg_s0.6")
na_lib.write_pdb(testbuild,"test_cgcg_b_cg_s0.6.pdb")

testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna_cgr")
na_lib.write_pdb(testbuild,"test_atat_b_cgr.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","abdna_cgr")
na_lib.write_pdb(testbuild,"test_cgcg_b_cgr.pdb")

testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna_cgx")
na_lib.write_pdb(testbuild,"test_atat_b_cgx.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","abdna")
na_lib.write_pdb(testbuild,"test_cgcg_b.pdb")

testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","abdna_cgx")
na_lib.write_pdb(testbuild,"test_cgcg_b_cgx.pdb")
#
testbuild = na_lib.fd_build_duplex("GCGCGCGCGCGC","abdna_cg")
na_lib.write_pdb(testbuild,"test_cgcg_b_cg.pdb")
#
#testbuild = na_lib.fd_build_duplex("ATATATATATAT","adna")
#na_lib.write_pdb(testbuild,"test_atat_a.pdb")
#
#testbuild = na_lib.fd_build_duplex("ATATATATATAT","adna_cg")
#na_lib.write_pdb(testbuild,"test_atat_a_cg.pdb")

sys.exit(0)

genseq = "AAACAGATTTTATCTGCCCGCTCAGGGCGTGA"
testbuild = na_lib.fd_build_duplex(genseq,"abdna")
na_lib.write_pdb(testbuild,"gen_abdna.pdb")

testbuild = na_lib.fd_build_duplex(genseq,"abdna_cg")
na_lib.write_pdb(testbuild,"gen_abdna_cg.pdb")

testbuild = na_lib.fd_build_duplex(genseq,"abdna_cgx")
na_lib.write_pdb(testbuild,"gen_abdna_cgx.pdb")
sys.exit(0)

testbuild = na_lib.fd_build_duplex("ATATATATATAT","abdna_cg")
na_lib.write_pdb(testbuild,"test_atat_b_cg.pdb")

testbuild = na_lib.fd_build_duplex("ATATATATATAT","adna")
na_lib.write_pdb(testbuild,"test_atat_a.pdb")

testbuild = na_lib.fd_build_duplex("ATATATATATAT","adna_cg")
na_lib.write_pdb(testbuild,"test_atat_a_cg.pdb")
