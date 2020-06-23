#!/bin/bash

usage(){ echo "                                             
Takes an fsc-style SFS and list of models (contained in templates folder) to build archicture for 
fastsimcoal model testing, and execute runs. 
                                                            
Usage: multiple_models.sh [-s|--sfs] [-p|--popsizes] [-o|--outdir] [-e|--execute] <model_1> <model_2> ... <model_n>
	
	-s : Input sfs, already formatted as needed for fastsimcoal. REQUIRED 
	-p : Diploid pop sizes, comma formatted (20,30). REQUIRED (for multiple pops at least)
	-o : Directory for output files. REQUIRED   
	-e : FLAG. If found, fsc runs will submit after constructing directories. 
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
execute="FALSE"
remove_old_files="FALSE"
##############


while :; do
	case $1 in
		-h|-h\?|--help)
			usage
			exit 1
			;;
		-s|--sfs)
			if [ "$2" ];then
				sfs=$2
				shift
			fi
			;;
                -p|--popsizes)
                        if [ "$2" ];then
                                popsizes=$2
                                shift
                        fi
                        ;;

		-o|--outdir)
			if [ "$2" ];then
				outdir=$2
				shift
			fi
			;;
		-e|--execute)
			if [ "$1" ];then
				execute="TRUE"
			fi
			;;
		-r|--remove)
			if [ "$1" ];then
				remove_old_files="TRUE"
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

sizes=$(echo $popsizes | awk -F"," '{OFS=" "} {$2=$2;print $0}') # $2=$2 tells awk to modify the input field $0. Otherwise it prints with the original comma separator
pops=$(echo $sizes | awk '{print NF}')

modelList="$@"

## Now start doing the work ##

mkdir -p ${outdir}
cp $sfs ${outdir}/

for model in $modelList;do
	if [[ $remove_old_files != "FALSE" ]];
		then rm -rf ${outdir}/${model}
	fi
	##
	cp -r templates/${model} ${outdir}/
	if [[ $model =~ "1Pop" ]];then
		cp $sfs ${outdir}/${model}/${model}_DAFpop0.obs
		dim=$(($(awk '{print NF}' $sfs | tail -n 1)-1))
		sed -i "s/diploid_num/$dim/" ${outdir}/${model}/${model}.tpl
	elif [[ $model =~ "2Pop" ]] || [[ $model =~ "INV" ]];then
		cp $sfs ${outdir}/${model}/${model}_DSFS.obs
		dip1=$(echo $sizes | cut -d" " -f1)
		dip2=$(echo $sizes | cut -d" " -f2)
		sed -i "s/diploid_num_pop1/$dip1/" ${outdir}/${model}/${model}.tpl
		sed -i "s/diploid_num_pop2/$dip2/" ${outdir}/${model}/${model}.tpl
	elif [[ $model =~ "3Pop" ]];then
		cp $sfs ${outdir}/${model}/${model}_DSFS.obs
		dip1=$(echo $sizes | cut -d" " -f1)
                dip2=$(echo $sizes | cut -d" " -f2)
		dip3=$(echo $sizes | cut -d" " -f3)
                sed -i "s/diploid_num_pop1/$dip1/" ${outdir}/${model}/${model}.tpl
                sed -i "s/diploid_num_pop2/$dip2/" ${outdir}/${model}/${model}.tpl
                sed -i "s/diploid_num_pop3/$dip3/" ${outdir}/${model}/${model}.tpl
	else
		cp $sfs ${outdir}/${model}/${model}_DSFS.obs
	fi
	
done

## Now submit
if [[ $execute != "FALSE" ]];then
	for model in $modelList;do
		cd ${outdir}/${model}
		echo -e "Starting run for ${outdir}/${model}"
		quick.sh -m 1000 -t 00:15:00 "sh /n/home11/twooldridge/scripts/launch_fsc_runs.sh $model"
		cd - &>/dev/null
	done
fi

