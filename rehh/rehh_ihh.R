#!/usr/bin/env Rscript
suppressMessages(library(optparse))


option_list = list(make_option(c("-t", "--thap"), type="character", default=NULL,help="input haplotype file"),
			make_option(c("-m", "--map"), type="character",default=NULL,help="input map (snp info) file"),
			make_option(c("-c", "--chr"), type = "character",default=NULL,help="chromosome name"),
			make_option(c("-l", "--lime"), type = "double", default=0.05,help = "cutoff at which to stop computing ihs"), 
			make_option(c("-o","--outfile"),type="character",default="results.txt",help="Path to output file"),
			make_option(c("-n","--nthreads"),type="integer",default=1,help="Number of cores for multithreading"))


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if(length(opt)==0){
        print_help(opt_parser)
        stop()
}



################# LIBS ###################
print("Loading libraries...")
suppressMessages(library(data.table))
suppressMessages(library(rehh))


############### PERFORM CHECKS ################
if (is.null(opt$thap) || is.null(opt$map) || is.null(opt$chr)){
	print("Missing one or more input files, and or chromosme name")
	stop()
}


if (is.null(opt$outfile)){
	print("Output file path needed")
	stop()
}

############## LOAD DATA & RUN #################
print("Loading data...")
haps = data2haplohh(min_perc_geno.mrk=50,min_maf=0.01,hap_file=opt$thap, map_file=opt$map, haplotype.in.columns=TRUE, allele_coding = 'map', chr.name=opt$chr, verbose=FALSE)
print("Running scan_hh...")
hh_stats = scan_hh(haps,threads = opt$nthreads, limehh = opt$lime, limehhs = opt$lime)



############# WRITE OUT RESULTS ################
print("Writing out results")
fwrite(hh_stats,file = opt$outfile,sep="\t",col.names=TRUE,row.names=FALSE,nThread = opt$nthreads)







