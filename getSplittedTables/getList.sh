#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'table' '' 'line entry of bedgraph table' 'e'
DEFINE_string 'outputDir' '' 'path to output folder' 'o'
DEFINE_string 'factor' '' '[optional] factor to generate files for only that factor' 'f'
DEFINE_string 'for' '' 'gives type coverage or metagenes to split table into' 't'
DEFINE_boolean 'debug' 'true' '[optional] prints out debug messages.' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'
DEFINE_boolean 'version' 'false' '[optional] print version' 'v'

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "parameters before check" # print param values, if in debug mode

# check if mandatory arguments are there
if [[ -z $FLAGS_outputDir ]]; then
	echo "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_for ]] ; then
	echo "Parameter -t must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_table ]]; then
	echo "Parameter -e must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi


### start script
IFS=',' read -r -a array <<< $FLAGS_factor
echo "${array[@]}"

mkdir -p $FLAGS_outputDir



if [[ $FLAGS_for == "coverage" ]] ; then
	counter=1
	while read line ; do
        	if [[ "$line" == "" ]] ; then #empty line
                	continue
        	fi
        	linearray=($line)
        	if [[ "${linearray[0]}" == "BEDGRAPH_POS" ]] ; then #header
                	continue
        	fi
        	
        	if [[ "${linearray[5]}" == "YES" ]] ; then   #use=yes
                	if [[ ($FLAGS_factor && (" ${array[@]} " =~ " ${linearray[2]} ")) || -z $FLAGS_factor ]] ; then 
                        	echo $line | sed s/" "/"\t"/g  > $FLAGS_outputDir"$counter"_table.csv
                        	((counter++))
                        	echo "print file "
                	fi
        	fi
	done < $FLAGS_table
elif [[ $FLAGS_for == "metagenes" ]] ; then
	factors=$(cat $FLAGS_table | cut -f3 | sort -u)
	experiments=$(cat $FLAGS_table | cut -f7 | sort -u)
	counter=1
	for f in $factors ; do
		if [[ ($FLAGS_factor && (" ${array[@]} " =~ " $f ")) || -z $FLAGS_factor ]] ; then
        		for exp in $experiments ; do
                		subtable=$(cat $FLAGS_table | awk -v col3="$f" -v col7="$exp" '$3==col3 && $7==col7 && $6=="YES" {print}')
				if [[ -z $subtable  ]] ; then
					continue
				fi
				cat $FLAGS_table | head -n1  > $FLAGS_outputDir"$counter"_table.csv
				echo "$subtable" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"""}' >> $FLAGS_outputDir"$counter"_table.csv
				((counter++))
                		echo "print file "
        		done
        	fi
	done
fi



writeParam2File "$FLAGS_returnFilePath" "list" "$FLAGS_outputDir"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK
