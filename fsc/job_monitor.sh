#!/bin/bash
function checkStatus(){
	status=$(sacct | grep $1.ex | awk '{print $5}')
	echo $status
}

jobList=$1
modelName=$(basename $1 | sed 's/runIDs.txt//;s/_\+$//g')

# The script will only reach this point when all jobs have finished
#rm -f incomplete
rm -f $modelName.results
header=$(head -n 1 run1/$modelName/$modelName.bestlhoods)
echo -e "$header\tAIC"  > $modelName.results
#current_best="-1000000000000"
for folder in run*/;do
	runid=$(echo $folder | sed 's/\///g')
	results_file="$folder/$modelName/$modelName.bestlhoods"
	if [ ! -f $results_file ];
	then 
		echo $folder | sed 's/\///g' >> incomplete;
	else
		results=$(tail -n+2 $results_file)
		new_score=$(echo $results | awk '{print $(NF-1)}')
		num_param=$(cat $modelName.k)
		AIC=$(echo -e "-2*$new_score+2*$num_param" | bc)
		#if (( $(echo "$new_score >= $current_best" | bc -l) ));then current_best=$new_score;best_run=$(echo $folder | sed 's/\///g');fi
		echo -e "$runid\t$results\t$AIC" >> $modelName.results
		#tail -n+2 $folder/$modelName/$modelName.bestlhoods >> $modelName.results;
	fi
done

# Now that all the jobs have finished, extract the best ML and calculate AIC	
#num_param=$(cat $modelName.k)
#AIC=$(echo -e "-2*$current_best+2*$num_param" | bc)
#echo -e "$best_run\t$current_best\t$AIC" > $modelName.AIC
