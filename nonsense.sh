#!/bin/bash


usage(){ echo "                                             
This is a nonsense script that by default adds the word "polychaete" to the end of each line in a text file. You can specify "abalone" instead with the -a|--abalone flag, or some other word with the -o|--other argument. 
                                                            
Usage: nonsense.sh [-w|--word] [-p|--print] <input_file>

        -w|--word    CHAR   : Append other word to end of each line in file [polychaete]
        -p|--print          : Print final command to console before running
        -h|--help           : Display this help message 
"                                                           
1>&2;exit 1;}



die() {
    printf '%s\n' "$1" >&2
    exit 1
}


# Set defaults. If arguments have been provided, these will be overrriden.
endword="polychaete"

# Initialize all the option variables.
# This ensures we are not contaminated by variables from the environment.
verbose=0

while :; do
        case $1 in
                -h|-h\?|--help)
                        usage
                        exit1
                        ;;
                -w|--word)
                        if [ "$2" ];then
                                endword=$2
                                shift
                        fi
                        ;;
                -p|--print)
                        print="placeholder"
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

## Now check if input text file has been provided
if [ -z $1 ];then
  echo "Input file missing! You dullard";
  exit;
else
  infile=$1
fi

## Now construct the command
CMD="sed 's/$/ ${endword}/g' ${infile}"

## Check print argument. If "print" exists, print the command
if [ ! -z $print ];then
  echo "The command is:" >&2 ## The '>&2' redirects this message to stderr. Good for a quick solution
  echo $CMD >&2
fi

## Finally, run the command
echo "Running..." >&2
eval $CMD
echo "Done!" >&2

##
