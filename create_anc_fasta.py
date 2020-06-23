#! /usr/bin/env python

from __future__ import print_function
import re, sys, os, sets, getopt
from itertools import izip        
import gzip
import re
import sys
import subprocess
import os
import argparse
import time

'''
Some code borrowed and modified from Sonal Singhal here: https://github.com/singhal/postdoc/blob/master/simple_ancestral_alleles.py
Need to run for each scaffold (chr)
'''

def get_vcf(vcffile, type):
        file = gzip.open(vcffile, 'r')
        var = dict()
        for i,l in enumerate(file):
                if not re.match("#", l):
                        d = re.split('\t', l)
                        alleles = [d[3]] + re.split(",", d[4])
                        # Only realy interested in biallelic sites
			if len(alleles) >3:
				continue
			# do not want to include indels in this variant dictionary
                        # cannot use these to polarize
			indel = False
                        for allele in alleles:
                            if len(allele) > 1:
                            	#Added this additional if statement because with only 1 sample GATK cannot distinguish the possibility of another allele. For our purposes we just set that to the reference allele though.
                            	if allele != "<NON_REF>":
                                	indel = True
                        
                        if not indel:
                                # do not want any alleles included that are at zero frequency
                                # not useful for anything but to show ref, and we already have that
                                allele_counts = dict()
                                for i in range(len(alleles)):
                                        #eprint(alleles[i])
					count = len(re.findall('\s%s[\/|\:]' % i, l)) + len(re.findall('\/%s' % i, l))
					#eprint(count)
					allele_counts[alleles[i]] = count
                                tot_n = sum(allele_counts.values())
                                #eprint(d[1])
				#eprint(tot_n)

				if tot_n > 0:
				
					af_real = dict()
	
					for i in allele_counts:

                	                        af = allele_counts[i] / float(tot_n)
                        	                if af > 0:
                                	        	af_real[i] = '%.3f' % af
					# for the ingroup, do not want to include fixed or polyallelic sites, so deal with that here
					if type == 'in':
						if len(af_real) <= 2:
                                			var[int(d[1])] = af_real
					else:
						var[int(d[1])] = af_real
        file.close()
        return var



## 08.28.18
## There is no longer a distinction between get_history and get_history_farout. 
## Both functions now return 'N's if a genotype is missing, regardless of how
## distant the outgroup is. See Allison's original script for details

def get_history(pos, var_out, ref):
	if pos in var_out:
		if len(var_out[pos]) == 1:
			return var_out[pos].keys()[0]	
		else:
			return var_out[pos]
	else:
		return 'N'

def get_history_farout(pos, var_out, ref):
	if pos in var_out:
		if len(var_out[pos]) == 1:
			return var_out[pos].keys()[0]
		else:
                        return var_out[pos]
	else:
		return 'N'


def check_equal(iterator):
	return len(set(iterator)) <= 1

#Added 5th outgroup (we have 3 far outgroups)
def call_ancestral_allele(b_in, b_out1, b_out2):

	near_out = 'N'
	#near_out = []
	far_out1 = 'N'
	#far_out2 = []
	
	#If not a dictionary, add base to near_out. If it is a dictionary (polymorphic), check whether both bases are shared with the ingroup. If only 1, add that to near_out. If both (or missing data), nothing added to near_out.
	
	poly = 1
	
	# If the outgroup is not polymorphic, store allele
	if not isinstance(b_out1, dict):
		if b_out1 in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']:
			near_out = b_out1
	# If outgroup is polymorphic, append all unique alleles
	if isinstance(b_out1, dict):
		poly = poly + 1
		joint_poly = list(b_out1.viewkeys() & b_in.keys())
		# The step below accounts for 'false' polymorphisms, where a missing genotype
		# '.' is present
		joint_poly = [x for x in joint_poly if x != '.']
		if len(joint_poly) >= 1:
			near_out = joint_poly
			
#	if not isinstance(b_out2, dict):
#		if b_out2 in ['A', 'T', 'C', 'G']:
#			near_out.append(b_out2)
#	if isinstance(b_out2, dict):
#		joint_poly = list(b_out2.viewkeys() & b_in.keys())
#		print joint_poly
#		if joint_poly == 1:
#			print 'here'
#			near_out.append(joint_poly[0])
			
	if not isinstance(b_out2, dict):
		if b_out2 in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']:
			far_out1 = b_out2
	if isinstance(b_out2, dict):
		poly = poly + 1
		joint_poly = list(b_out2.viewkeys() & b_in.keys())
		# The step below accounts for 'false' polymorphisms, where a missing genotype
                # '.' is present
                joint_poly = [x for x in joint_poly if x != '.']
		if len(joint_poly) >= 1:
			far_out1 = joint_poly
					
#	if not isinstance(b_out3, dict):
#		if b_out3 in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']:
#			far_out2.append(b_out3)
#	if isinstance(b_out3, dict):
#		poly = poly + 1
#		joint_poly = list(b_out3.viewkeys() & b_in.keys())
#		if len(joint_poly) == 1:
#			far_out2.append(joint_poly[0])	

#	if not isinstance(b_out5, dict):
#		if b_out5 in ['A', 'T', 'C', 'G']:
#			far_out2.append(b_out5)
#	if isinstance(b_out5, dict):
#		joint_poly = list(b_out5.viewkeys() & b_in.keys())
#		if joint_poly == 1:
#			far_out2.append(joint_poly[0])	
	

##############################################################################################
### THIS SECTION CONTAINS THE ACTUAL PARSIMONY LOGIC AND DETERMINATION OF ANCESTRAL ALLELE ###
##############################################################################################

	#Calculate near out ancestor
	if len(near_out) >= 2:
		if check_equal(near_out):
			near_out_anc = near_out[0]
		else:
			near_out_anc = near_out
	elif len(near_out) == 1:
		near_out_anc = near_out[0]
	else:
		near_out_anc = 'N'
		

	#Calculate far_out ancestor
	if len(far_out1) >= 2:
		if check_equal(far_out1):
			far_out1_anc = far_out1[0]
		else:
			far_out1_anc = far_out1
	elif len(far_out1) == 1:
		far_out1_anc = far_out1[0]
	else:
		far_out1_anc = 'N'
	
	# near_out_anc is NOT polymorphic
	if not isinstance(near_out_anc, list):
		# Simplest case, where both outgroups have the same allele, and allele is in ingroup
		if near_out_anc == far_out1_anc and near_out_anc in b_in:
			allele = near_out_anc
		# Near_anc is different or uncalled (N) from far_anc, then... 
		elif near_out_anc != far_out1_anc:
			# If near_anc is found in the ingroup, case closed
			if near_out_anc in b_in:
				allele = near_out_anc
			# If near_anc is not in ingroup (so therefore 'N' or other) let's check the next outgroup...
			else:
				# If that outgroup is polymorphic, we can't safely determine the ancestral allele
				if isinstance(far_out1_anc, list):
					allele = 'N'
				# If the outgroup is NOT polymorphic, and it's in the ingroup, then...	
				else:
					if far_out1_anc in b_in:
						allele = far_out1_anc
					else:
						allele = 'N'
				# 
				# Worst case: neither share alleles with b_in
				# This can happen if both outgroups are 'N', or a situation like:
				# Ingroup: A|T, out1 = 'G', out2 = 'C'
		else:
			allele = 'N'
	
	# near_out_anc IS polymorphic
	elif isinstance(near_out_anc, list):
		# If far_out1_anc is also poly, then forget it
		if isinstance(far_out1_anc, list):
			allele = 'N'
		# If monomorphic, then proceed
		else:
			# If far_out_anc is one of the alleles in near_out_anc...
			if far_out1_anc in near_out_anc:
				# And if far_out_anc can also be found in the ingroup...
				if far_out1_anc in b_in:
					allele = far_out1_anc
				# But if far_out_anc is not in the ingroup ('N' or otherwise), the anc. cannot be determined
				else:
					allele = 'N'
			else:
				allele = 'N'

	else:
		allele = 'N'
	
	return allele, poly, near_out_anc, far_out1_anc


		
def get_chromosome(genome, chr, out_stem):
	outfile = out_stem + '/fastas/' + os.path.basename(genome) + '_' + chr
	subprocess.call('samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
	out_f = open(outfile, 'r')
	locus_name = out_f.next()
	return outfile, out_f



#Make sure to send this to the vcf reader
def get_chromosome_vcf(vcf_file, species, chr, out_stem):
	outfile = out_stem + '/vcfs/' + species + '_' + chr + '.vcf.gz'
	#subprocess.call('/n/home13/ashultz/sw/progs/vcftools-0.1.15/bin/vcftools --gzvcf %s --chr %s --recode --stdout | gzip -c > %s' %(vcf_file, chr, outfile), shell=True)
	#subprocess.call('/Users/allisonshultz/miniconda2/bin/vcftools --gzvcf %s --chr %s --recode --stdout | gzip -c > %s' %(vcf_file, chr, outfile), shell=True)
	cmd=('~/Software/tabix/tabix %s %s | ~/Software/tabix/bgzip -c > %s' %(vcf_file, chr, outfile))
	subprocess.call(cmd, shell=True)
	return outfile

def pp_hist(base_info):
	if isinstance(base_info, dict):
		return "|".join("%s:%s" % bp_pair for bp_pair in base_info.items())
	else:
		return base_info

def trawl_genome(out_file, anc_fasta_file, chr, f_ref, var_in, var_out1, var_out2, start_pos):
	out = open(out_file, 'w')
	anc_fasta = open(anc_fasta_file, "w")
	out.write('chr,position,ingroup,outgroup1,outgroup2,ancestral,poly\n')
	anc_fasta.write('>%s\n' %(chr))
	
	eprint("Beginning trawling...")

	for line_count, (l_ref) in enumerate(f_ref):
		l_ref = l_ref.strip()
		
		if len(l_ref) > 60:
#			print 'Fasta lines not 60 bases, oh no!'
			quit()

		for ix, (b_ref) in enumerate(l_ref):
			# get the current genome position, need to iterate through reference and get each base value.	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = (start_pos) + (line_count * 60) + (ix)
			# this is a variable position
			if pos in var_in:
				#eprint((str(pos) + ' found'))
				b_in = var_in[pos]
				## VERY IMPORTANT ##
				## If input vcfs inclued all sites, use get_history_farout. This function will rightfully return an 
				## 'N' if no alternate alleles are presented (nothing mapped). If the input is a typical VCF, use get_history. 
				## This function will return the reference fasta allele, as it assumes that no alternate allele at that position
				## is more likely the result of a real variant call, rather than mapping error.
				b_out1call = get_history_farout(pos, var_out1, b_ref)
				b_out2call = get_history_farout(pos, var_out2, b_ref)
				#b_out3call = get_history_farout(pos, var_out3, b_ref)
				#b_out4call = get_history_farout(pos, var_out4, b_ref)
				#b_out5call = get_history_farout(pos, var_out5, b_ref)
				
				eprint(pos,b_in,b_out1call,b_out2call)	
				anc_allele, poly, near_out_anc, far_out1_anc = call_ancestral_allele(b_in, b_out1call, b_out2call)

				out.write('%s,%s,%s,%s,%s,%s,%s\n' % \
					(chr, pos, pp_hist(b_in), pp_hist(b_out1call), pp_hist(b_out2call), anc_allele, poly))
				anc_fasta.write(anc_allele)
				

			else:
				anc_fasta.write(b_ref)
		anc_fasta.write("\n")

	out.close()
	anc_fasta.close()
	return


def clean_up(chr_out, f_out):
	f_out.close()
	os.remove(chr_out)
	
def clean_up_vcf(chr_out_vcf):
	os.remove(chr_out_vcf)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def main():
	eprint("Beginning...")
	parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        parser.add_argument("--ref", help="reference genome fasta for which to run analysis")
        parser.add_argument("--vcf_in", help="vcf for ingroup for which to run analysis")
        parser.add_argument("--vcf_out1", help="vcf for outgroup 1 (CC) for which to run analysis")
        parser.add_argument("--vcf_out2", help="vcf for outgroup 2 (CP) for which to run analysis")
        parser.add_argument("--vcf_out3", help="vcf for outgroup 3 (PinEnu) for which to run analysis")
        parser.add_argument("--vcf_out4", help="vcf for outgroup 4 (CarEry) for which to run analysis")
        parser.add_argument("--vcf_out5", help="vcf for outgroup 5 (UraSib) for which to run analysis")
        parser.add_argument("--out", help="output file path and stem for results")
        	
        args = parser.parse_args()
        chr = args.chr
        genome_ref = args.ref
        vcf_in = args.vcf_in
        vcf_out1 = args.vcf_out1
        vcf_out2 = args.vcf_out2
        vcf_out3 = args.vcf_out3
        vcf_out4 = args.vcf_out4
        vcf_out5 = args.vcf_out5
        out_stem = args.out
	
	subprocess.Popen('mkdir -p ' + out_stem,shell=True)
	subprocess.Popen('mkdir -p ' + out_stem + '/fastas',shell=True)
	subprocess.Popen('mkdir -p ' + out_stem + '/vcfs',shell=True)
	subprocess.Popen('mkdir -p ' + out_stem + '/anc_fastas',shell=True)
	subprocess.Popen('mkdir -p ' + out_stem + '/anc_csvs',shell=True)
	time.sleep(5)

	eprint("Arguments parsed...")

	#For real script will need to add a step to run it through the vcftools command to select the scaffold of interest.
	#HF
 	#vcf_in = 'All_HFs_testfile.vcf.gz'
 	sp_in = "MANIC"
 	#CC and CP
	#vcf_out1 = 'All_CCs_testfile.vcf.gz'
	sp1 = 'GO'
	#vcf_out2 = 'All_CPs_testfile.vcf.gz'
	sp2 = 'IS'
	#Distant finches
	#vcf_out3 = 'PinEnu_1.vcf.gz'
	#sp3 = 'IS'
	#vcf_out4 = 'CarEry_1.vcf.gz'
	#sp4 = 'CarEry'
	#vcf_out5 = 'UraSib_1.vcf.gz'
	#sp5 = 'UraSib'

	#Reference (HF) fasta and output files (output csv and new fasta with ancestral alleles)
 	#genome_ref = 'final.assembly.homo.fa'
 	out_file = '%s/anc_csvs/%s.csv' % (out_stem,chr)
 	anc_fasta_file = '%s/anc_fastas/%s.fasta' % (out_stem,chr)
	
	eprint(("Reading in vcf for %s" % sp_in))
	vcf_chr_in = get_chromosome_vcf(vcf_in, sp_in, chr, out_stem)
 	var_in = get_vcf(vcf_chr_in, 'in')

	eprint(("Reading in vcf for %s" % sp1))	
	vcf_chr_out1 = get_chromosome_vcf(vcf_out1, sp1, chr, out_stem) 	
 	var_out1 = get_vcf(vcf_chr_out1, 'out')	

	eprint(("Reading in vcf for %s" % sp2))
	vcf_chr_out2 = get_chromosome_vcf(vcf_out2, sp2, chr, out_stem) 	
 	var_out2 = get_vcf(vcf_chr_out2, 'out')

	#eprint(("Reading in vcf for %s" % sp3))
	#vcf_chr_out3 = get_chromosome_vcf(vcf_out3, sp3, chr, out_stem)
 	#vcf_chr_out3 = 'PinEnu_scaffold_0_test.vcf.gz'
 	#var_out3 = get_vcf(vcf_chr_out3, 'out')

 	
 	#vcf_chr_out4 = get_chromosome_vcf(vcf_out4, sp4, chr)
 	#vcf_chr_out4 = 'CarEry_scaffold_0_test.vcf.gz'
 	#var_out4 = get_vcf(vcf_chr_out4, 'out')

 	#vcf_chr_out5 = get_chromosome_vcf(vcf_out5, sp5, chr)
 	#vcf_chr_out5 = 'UraSib_scaffold_0_test.vcf.gz'
 	#var_out5 = get_vcf(vcf_chr_out5, 'out')

	
	start = chr.split(':')[1]
	start = int(start.split('-')[0])

 	chr_ref, f_ref = get_chromosome(genome_ref, chr, out_stem)

 	trawl_genome(out_file, anc_fasta_file, chr, f_ref, var_in, var_out1, var_out2, start_pos = start)
 	
##Need to add clean up for vcfs - remove scaffold files that have been created.

 	clean_up_vcf(vcf_chr_in)
 	clean_up_vcf(vcf_chr_out1)
 	clean_up_vcf(vcf_chr_out2)
 	#clean_up_vcf(vcf_chr_out3)
 	#clean_up_vcf(vcf_chr_out4)
  	#clean_up_vcf(vcf_chr_out5)
 	clean_up(chr_ref, f_ref)



if __name__ == "__main__":
	main()
