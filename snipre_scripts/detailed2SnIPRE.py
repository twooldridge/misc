#! /usr/bin/env python

import os, sys, io, re, argparse

'''
This script takes the "detailed" output from vcf2MK (https://github.com/russcd/vcf2MK), and formats it for input into SnIPRE.

Output format:
"geneID"	"PR"	"FR"	"PS"	"FS"	"Tsil"	"Trepl"	"nout"	"npop"
geneID = transcriptID
PR = polymorphic replacement (nonsynonymous)
FR = fixed replacement
PS = polymorphic synonymous
FS = fixed synonymous
Tsil = # silent sites
Trepl = # replacement sites
nout = outgroup sample size (we assume 2, as estimated from 1 diploid individual)
npop = population sample size (average used to estimate for each gene, rounded to the nearest integer)

Inputs:
--detailed vcf2MK detailed output
--out output file path for results
--maf desired maf cutoff for inclusion as polymorphic site (default 0.1)
'''

def	create_mk_dict(detailed_file):
	#Create a dictionary of transcripts, with each transcript its own dictionary of variables to be summarized
	mk_dict = {}
	detailed = open(detailed_file,"r")
	
	for l in detailed:
		l = l.strip().split("\t")
		trans = l[11].split(",")[0:-1]
		for tran in trans:
			if tran not in mk_dict:
				mk_dict[tran] = {}
				mk_dict[tran]["PR"] = [int(l[7])]
				mk_dict[tran]["FR"] = [int(l[5])]
				mk_dict[tran]["PS"] = [int(l[8])]
				mk_dict[tran]["FS"] = [int(l[6])]
				mk_dict[tran]["freq"] = [l[10]]
				mk_dict[tran]["n"] = [int(l[9])]
				mk_dict[tran]["foldDegen"] = [int(l[12])]
			else:
				mk_dict[tran]["PR"].append(int(l[7]))
				mk_dict[tran]["FR"].append(int(l[5]))
				mk_dict[tran]["PS"].append(int(l[8]))
				mk_dict[tran]["FS"].append(int(l[6]))
				mk_dict[tran]["freq"].append(l[10])
				mk_dict[tran]["n"].append(int(l[9]))
				mk_dict[tran]["foldDegen"].append(int(l[12]))
	detailed.close()
	return(mk_dict)

def maf_correction(maf_cutoff, maf_list, poly_list):
	#This function takes a list of minor allele frequencies, polymorphism scores (0 or 1), and minor allele frequency cutoffs. If a site falls below the cutoff but is listed as polymorphic (1), it will be corrected to 0.
	new_poly_list = []
	maf_cutoff = float(maf_cutoff)
	for i in range(0,len(poly_list)):
		if poly_list[i] == 1:
			if ((float(maf_list[i]) < maf_cutoff) or (float(maf_list[i]) > (1-maf_cutoff))):
				new_poly_list.append(0)
			else:
				new_poly_list.append(1)
		else:
			new_poly_list.append(poly_list[i])
	return(new_poly_list)

def mk_summarize(mk_dict,maf_cutoff):
	mk_dict_sum = {}
	for tran in mk_dict:
		tran_sum = []
		tran_sum.append(str(sum(maf_correction(maf_cutoff,mk_dict[tran]["freq"],mk_dict[tran]["PR"]))))
		tran_sum.append(str(sum(mk_dict[tran]["FR"])))
		tran_sum.append(str(sum(maf_correction(maf_cutoff,mk_dict[tran]["freq"],mk_dict[tran]["PS"]))))
		tran_sum.append(str(sum(mk_dict[tran]["FS"])))
		
		#Calculate number of silent and replacement sites
		tran_sil = []
		tran_repl = []
		for site in mk_dict[tran]["foldDegen"]:
			if site == 0:
				tran_repl.append(1)
			elif site == 2:
				tran_repl.append(0.5)
				tran_sil.append(0.5)
			elif site == 3:
				tran_repl.append(0.5)
				tran_sil.append(0.5)				
			elif site == 4:
				tran_sil.append(1)
			else:
				print "Weird degeneracy score " +str(site)+" for "+tran
		tran_sum.append(str(sum(tran_sil)))
		tran_sum.append(str(sum(tran_repl)))
		
		#Assumes 1 diploid outgroup
		tran_sum.append("2")
		
		#Avg sample size for each transcript
		n = sum(mk_dict[tran]["n"])/len(mk_dict[tran]["n"])
		tran_sum.append(str(n))
		
		mk_dict_sum[tran] = tran_sum
	return(mk_dict_sum)
	

def main():
	parser = argparse.ArgumentParser()
        parser.add_argument("--detailed", help="vcf2MK detailed output")
        parser.add_argument("--out", help="output file path for results")
        parser.add_argument("--maf", help="maf cutoff for inclusion as polymorphic site, default 0.1")
        
        args = parser.parse_args()
        detailed_file = args.detailed
        out_file = args.out
        maf = args.maf
        
        if not maf:
        	maf = 0.1

	mk_dict = create_mk_dict(detailed_file)
	mk_dict_sum = mk_summarize(mk_dict,maf)
	
	trans = sorted(mk_dict_sum.keys())

	out = open(out_file,"w")
	out.write("geneID\tPR\tFR\tPS\tFS\tTsil\tTrepl\tnout\tnpop\n")
	
	for tran in trans:
		out.write("%s\t%s\n"%(tran,"\t".join(mk_dict_sum[tran])))

	out.close()

if __name__ == "__main__":
    main()
