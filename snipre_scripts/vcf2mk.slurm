#!/bin/bash

#SBATCH --mem=10000
#SBATCH --time=10:00:00
#SBATCH -p general,hoekstra,shared,holyhoekstra
#SBATCH -N 1
#SBATCH -n 1

ingroup=$1
outgroup=$2
gff=$3

id=$(basename $ingroup | sed 's/dp.//;s/.allSites.*//')
mkdir -p mk_files

~/Software/vcf2MK/vcf2MK -i -v <(cut -f1-8 $ingroup) -o <(cut -f1-8 $outgroup) -c $gff > mk_files/${id}.mk.raw  
