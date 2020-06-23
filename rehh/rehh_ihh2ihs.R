#!/usr/bin/env Rscript
suppressMessages(library(optparse))

# Args should be list of IHS/IES files, by chromosme for a single population
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
        print("No files provided")
        stop()
}

################# LIBS ###################
print("Loading libraries...")
suppressMessages(library(data.table))
suppressMessages(library(rehh))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
############## LOAD DATA & RUN #################
print("Loading data...")
print(paste0("reading file ",args[1]))
pop.res = fread(args[1],header=TRUE)
ihs = ihh2ihs(pop.res,freqbin=0,min_maf=0)
save(ihs,file=paste0(args[2],".RData"))
chroms = unique(ihs[[1]]$CHR)
for (chrom in chroms){
	ihs_chrom = list(ihs[[1]] %>% filter(CHR==chrom),ihs[[2]])
	save(ihs_chrom,file=paste0(args[2],"_",chrom,".RData"))
}


## Save plot ##
#layout(matrix(1:2,2,1))
#pdf(file=paste0(args[2],".pdf"),width=12,height=7)
#ihsplot(ihs,plot.pval=TRUE,ylim.scan=4)
#dev.off()


## Format artificially as a bed-type file, for masking
ihs[[1]] %<>% mutate(start = POSITION - 1,end = POSITION)
ihs[[1]] %<>% dplyr::select(CHR,start,end,POSITION,IHS,LOGPVALUE)
print("Data combined, writing out total file")
options(scipen = 999)
write.table(ihs[[1]],file = paste0(args[2],".bed"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
print("Filtering datasets")
high=ihs[[1]] %>% filter(IHS > quantile(na.rm=TRUE,IHS,0.995))
med=ihs[[1]] %>% filter(IHS > quantile(na.rm=TRUE,IHS,0.99))
low=ihs[[1]] %>% filter(IHS > quantile(na.rm=TRUE,IHS,0.95))
print("Writing filtered datasets")
write.table(high,file = paste0(args[2],"_995.bed"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(med,file = paste0(args[2],"_99.bed"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(low,file = paste0(args[2],"_95.bed"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)



