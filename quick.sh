#!/bin/bash

usage(){ echo "                                             
Usage: quick.sh [-m|--memory] [-t|--time] [-n|--threads] [-i|--jobid] "command.here"
	
	-m: job memory
	-t: job time
	-n: threads
	-i: jobid, preprended to *err and *out files                                                            
"                                                           
1>&2;exit 1;}

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
verbose=0
###Defaults###
mem=10000
jobtime=05:00:00
cores=1
##############

while :; do
	case $1 in
		-h|-\?|--help)
			usage    # Display a usage synopsis.
            		exit1
			;;
		-m|--memory)
			if [ "$2" ];then
				mem=$2
				shift
			fi
			;;
		-t|--time)
			if [ "$2" ];then
				jobtime=$2
				shift
			fi
			;;
		-n|--threads)
			if [ "$2" ];then
				cores=$2
				shift
			fi
			;;
		-i|--jobid)
			if [ "$2" ];then
				id=$2
				shift
			fi
			;;
		--)
			shift
			break
			;;
		-?*)	
			printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
			;;
		*)
			break
			;;
	esac
	shift
done


if [[ -z $1 ]];then echo "No commmand provided, exiting";exit;fi

BASECMD="--mem=$mem -t $jobtime -N 1 -n $cores -p commons,serial_requeue,hoekstra,shared --wrap=\"$1\""

if [[ ! -z $id ]];then FULLCMD="sbatch -J ${id} -e ${id}.err -o ${id}.out $BASECMD";else FULLCMD="sbatch $BASECMD";fi

echo $FULLCMD
eval $FULLCMD
