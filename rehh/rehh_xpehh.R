#!/usr/bin/env Rscript
suppressMessages(library(optparse))


option_list = list(make_option(c("-a", "--pop1"), type="character", default=NULL,help="input iES/iHS file for population 1"),
			make_option(c("-b", "--pop2"), type="character",default=NULL,help="input iES/iHS file for population 2"),
			make_option(c("-o","--outprefix"),type="character",default="results.txt",help="Prefix for output files. Output will be written as ${prefix}_rsb.txt and ${prefix}_xpehh.txt"),
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
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

############### PERFORM CHECKS ################
if (is.null(opt$pop1) || is.null(opt$pop2)){
	print("Missing one or more input files")
	stop()
}


if (is.null(opt$outprefix)){
	print("Output file prefix needed")
	stop()
}

############## LOAD AND ANALYZE DATA  #################
print("Loading data...")
pop1.res = fread(opt$pop1,sep="\t",fill=TRUE,header=TRUE)
pop2.res = fread(opt$pop2,sep="\t",fill=TRUE,header=TRUE)

#if (nrow(pop1.res) != nrow(pop2.res)){
#	print("Pop1 and Pop2 files do not have the same number of rows. Have the same positions been analyzed in each dataset?")
#	stop()
#}

print("Calculating cross population statistics...")

print("Calculating rsb...")
rsb = ines2rsb(pop1.res,pop2.res)
save(rsb,file=paste0(opt$outprefix,"_rsb.RData"))
#pdf(file=paste0(opt$outprefix,"_rsb.pdf"))
#rsbplot(rsb,plot.pval=TRUE)
#dev.off()

print("Calculating xpehh...")
xpehh = ies2xpehh(pop1.res,pop2.res)
save(xpehh,file=paste0(opt$outprefix,"_xpehh.RData"))
#pdf(file=paste0(opt$outprefix,"_xpehh.pdf"))
#rsbplot(xpehh,plot.pval=TRUE)
#dev.off()


## Transform to bed ##
xpehh %<>% mutate(start=POSITION-1,end=POSITION) %>% dplyr::select(CHR,start,end,POSITION,contains("XPEHH"),LOGPVALUE)
rsb %<>% mutate(start=POSITION-1,end=POSITION) %>% dplyr::select(CHR,start,end,POSITION,contains("Rsb"),LOGPVALUE)


############# WRITE OUT RESULTS ################
print("Writing out results")
options(scipen = 999)
write.table(rsb,na = "NA",file = paste0(opt$outprefix,"_rsb.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(rsb %>% filter(Rsb > quantile(na.rm=TRUE,Rsb,0.995)),na = "NA",file = paste0(opt$outprefix,"_rsb_995.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(rsb %>% filter(Rsb > quantile(na.rm=TRUE,Rsb,0.99)),na = "NA",file = paste0(opt$outprefix,"_rsb_99.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(rsb %>% filter(Rsb > quantile(na.rm=TRUE,Rsb,0.95)),na = "NA",file = paste0(opt$outprefix,"_rsb_95.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

write.table(xpehh,na = "NA",file = paste0(opt$outprefix,"_xpehh.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(xpehh %>% filter(XPEHH > quantile(na.rm=TRUE,XPEHH,0.995)),na = "NA",file = paste0(opt$outprefix,"_xpehh_995.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(xpehh %>% filter(XPEHH > quantile(na.rm=TRUE,XPEHH,0.99)),na = "NA",file = paste0(opt$outprefix,"_xpehh_99.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(xpehh %>% filter(XPEHH > quantile(na.rm=TRUE,XPEHH,0.95)),na = "NA",file = paste0(opt$outprefix,"_xpehh_95.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)






