#!/bin/bash
#SBATCH -p general,holyhoekstra # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node
#SBATCH -t 05:00:00	  # Running time of 2 hours
#SBATCH --mem 10000	  # Memory request
#SBATCH -o postProcessing.out # Standard output
#SBATCH -e postProcessing.err # Standard error

while getopts ":t:k:s:" opt; do
 case $opt in
  t)
   inputTable=$OPTARG
   echo "Your input data table is $OPTARG" >&2
   ;;
  k)
   keyword=$OPTARG
   echo "Your search term is $OPTARG" >&2
   ;;
  s)
   sequenceFile=$OPTARG
   echo "Will create sequence database from $OPTARG" >&2
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

module load R
cat $inputTable | grep $keyword > key_hits.txt
Rscript --vanilla /n/home11/twooldridge/scripts/R/command_blast.R key_hits.txt candidates.txt
mkdir -p DBSE
makeblastdb -in $sequenceFile -out DBSE -dbtype nucl -parse_seqids
mv DBSE.n* DBSE/
blastdbcmd -db DBSE/DBSE -entry_batch candidates.txt -out=my_"$keyword"_seqs.fasta -outfmt '%f'

rm key_hits.txt candidates.txt
