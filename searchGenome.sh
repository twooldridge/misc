#!/bin/bash
#SBATCH -p serial_requeue # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node
#SBATCH -t 10:00:00         # Running time of 2 hours
#SBATCH --mem 10000        # Memory request
#SBATCH -o searchGenome.out # Standard output
#SBATCH -e searchGenome.err # Standard error

module load R
module load bedtools2

usage="$(basename "$0") [-h] [-q query] [-e interval extension] [-g genome] [-d database] [-c chromosomeList]  [-b blast table] -- program to find matches for amino acid sequence(s) in a genome. Uses tblastn and only retains hits with e-value less than or equal to 1e-10

where:
    -h  show this help text
    -q	sequence(s) to query genome, in fasta format
    -e  # basepairs to extend BED intervals by on each side. Default=1000
    -g  genome to search
    -d  path to local BLAST database of searched genome, in form of ./databasePrefix
    -c  genome file describing length of each chromosome or scaffold. See bedtools slop documentation for details
    -b  blast table, if provided, option q not needed. Will skip the database building and search, proceeding straight to formatting and sequence extraction. Allows for the input of custom tables based on specific user preferences. "

while getopts ":hq:e:g:d:c:b:" opt; do
 case $opt in
  h)
   echo "$usage"
   exit
   ;;
  q)
   sequences=$OPTARG
   ;;
  e)
   extension=$OPTARG
   ;;
  g)
   genome=$OPTARG
   ;;
  d)
   database=$OPTARG
   ;;
  c)
   chromosomeList=$OPTARG
   ;;
  b)
   blast_output=$OPTARG
   ;;
  \?)
   echo "Invalid option: -$OPTARG" >&2
   exit 1
   ;;
  :)
   echo "Option -$OPTARG requires an argument." >&2
   exit 1
   ;;
 esac
done


if [ -z $genome ]; then
 echo "Need genome for blast search and extracting  query sequences, provide with option -g"
 exit 1
fi
 

if [ -z $blast_output ]
then

 if [ -z $sequences ]; then
  echo "Missing query sequences, provide in fasta format with option q"
  exit  
 fi

 if [ -z $database ]; then
  echo "No blast database provided, making one now"
  makeblastdb -in $genome -out DBSE -dbtype nucl -parse_seqids
  database=DBSE

 fi
 echo "Beginning blast"
 tblastn -query $sequences -out blast_output -outfmt '6 qseqid sseqid sstart send evalue bitscore score length' -db $database
 #while [[ $(sacct | grep blast_9.fa | grep -e RUNNING -e PENDING) ]];do  echo "Blast still running";sleep 5m; done

else
 echo "Blast table already provided, skipping to formatting"
fi

echo "Beginning formatting...ordering intervals and removing hits for evalue less greater than 1e-10. See /n/home11/twooldridge/scripts/R/genomeBlast.R to change" 
Rscript --vanilla /n/home11/twooldridge/scripts/R/genomeBlast.R blast_output regions.bed
sort -k1,1 -k2,2n regions.bed | awk '{$4="";gsub(/ +/," ")}1' > sortedRegions.bed
perl -p -i -e 's/ /\t/g' sortedRegions.bed
bedtools merge -i sortedRegions.bed > mergedRegions.bed
if [ -z $extension ]; then
 echo "No interval extension parameter was supplied; extending blast hit intervals by 1000bp on each side"
 extension=1000
fi
if [ -z $chromosomeList ]; then
 echo "No chromosome list provided for the genome, cannot work with bed intervals. Please provide with option -c"
 exit 1
fi
bedtools slop -i mergedRegions.bed -g $chromosomeList -b $extension > extendedRegions.bed
bedtools merge -i extendedRegions.bed > finalRegions.bed
mkdir -p BEDfiles
mv *.bed BEDfiles/

echo "Formatting complete, extracting regions from genome..."
bedtools getfasta -fi $genome -bed BEDfiles/finalRegions.bed -fo extractedSequences.fasta
echo "Done!"

