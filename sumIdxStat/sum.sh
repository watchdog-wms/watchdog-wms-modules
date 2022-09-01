#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'inputFile' '' 'path to input file' 'i'
DEFINE_string 'excludeChrom' '' 'chromosome to exclude from sum' 'c'
DEFINE_string 'outputFile' '' 'path to output file' 'o'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "parameters before check" # print param values, if in debug mode

# check if mandatory arguments are there
if [[ -z $FLAGS_inputFile ]]; then
	echo "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_outputFile ]] ; then
	echo "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi


### start script

IFS=',' read -ra chromosomes <<< "$FLAGS_excludeChrom"
echo ${chromosomes[@]}

past_sample="XXX"
sample_sum=0

while read line ; do
	if [[ $line == *"contigt"* ]] ; then #header
		continue
	fi
	chr=$(echo $line | cut -f1 -d' ')
	mapped=$(echo $line | cut -f3 -d' ')
	sample=$(echo $line | cut -f5 -d' ')
	if [[ "$past_sample" != "$sample" ]] ; then
		if [[ "$past_sample" != "XXX" ]] ; then
			echo -e "$past_sample\t$sample_sum" >> $FLAGS_outputFile	
		fi
		past_sample=$sample
		if [[ ! " ${chromosomes[@]} " =~ " ${chr} " ]] ; then
			sample_sum=$mapped
		else
			sample_sum=0
		fi
	else
		if [[ ! " ${chromosomes[@]} " =~ " ${chr} " ]] ; then
			sample_sum=$(($sample_sum+$mapped))
		fi
	fi
done < $FLAGS_inputFile
echo -e "$past_sample\t$sample_sum" >> $FLAGS_outputFile



writeParam2File "$FLAGS_returnFilePath" "samplesSum" "$FLAGS_outputFile"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK
