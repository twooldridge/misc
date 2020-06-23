#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
# test if output directory for intersected bed files is supplied. If not, assuming present directory
if (length(args)==0) {
  # Default output directory #
  stop("No input file provided")
} else if (length(args)==1){
	print("No output prefix provided;using basename")
	outname=basename(args[1])
} else {
	outname=args[2]
}

input_snipre=args[1]
data = read.table(input_snipre,header=TRUE)
system(paste("mkdir -p ",outname,sep=""))
#system(paste("cd ",outname,sep=""))
system(paste("cp ~/Software/SnIPRE/SnIPRE.bug ",outname,"/",sep=""))
setwd(paste(getwd(),"/",outname,"/",sep=""))
getwd()

## Part (1)  Empirical Bayes Implementation  (lme4 package, SnIPRE_source.R)
## Part (2)  Bayesian Implementation (R2WinBUGS package, B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary

#################################################################
## Part (1)  Empirical Bayes Implementation  (lme4 package)
#################################################################

source("~/Software/SnIPRE/SnIPRE_source.R")
source("~/Software/SnIPRE/my.jags2.R")
library(lme4)
library(R2jags)
library(arm)


#SnIPRE <-function(mydata)
# mydata: name of data set;
# mydata must have a header with the following columns: PS, PR, FS, FR, npop, nout, Tsil, Trepl (no particular order)
# outputs 2 objects:  new.dataset & model
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
eb.res = SnIPRE(data)

res = eb.res$new.dataset
model = eb.res$model

ebfile=paste("eb_",outname,".csv",sep="")
write.table(res, file = ebfile, sep  = ",", row.names = FALSE)



#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################

source("~/Software/SnIPRE/B_SnIPRE_source.R")
source("~/Software/SnIPRE/my.jags2.R")
library(lme4)
library(R2jags)
library(arm)

#BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
  # path will be where the chains are stored, and must also be where the ".bug" model is located
  # burnin, thin, and iter (number iterations after burnin) are for MCMC samples
BSnIPRE.run(data, burnin = 10000, thin = 4, iter = 15000)

# check to make sure it finished correctly:
# if a "sample" file is in your working directory (getwd()), or the path you sepecified)
# is empty or not there, there is a problem


load("samples")

res.mcmc <- samples

#BSnIPRE <- function(data.mcmc,mydata){
# outputs 2 objects:  new.dataset & effects
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
# the "effects" may be useful if you are interested in estimation
# of population parameters (gamma, constraint) with other assumptions than the PRF

b.res <- BSnIPRE(res.mcmc, data)

bres = b.res$new.dataset

bresfile=paste("bayesian_",outname,".csv",sep="")
write.table(bres, file = bresfile, sep  = ",", row.names = FALSE)

