#!/bin/sh

module load gcc/6.2.0-fasrc01 JAGS/4.1.0-fasrc01 R/3.4.2-fasrc02
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
export PATH="/n/home11/twooldridge/Software/SnIPRE/:$PATH"
