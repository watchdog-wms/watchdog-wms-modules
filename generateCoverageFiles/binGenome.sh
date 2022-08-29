#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'outputDir' '' 'path to output folder' 'o'
DEFINE_string 'bedgraphTable' '' 'path to table with bedgprahp paths' 'g'
DEFINE_string 'bedfile' '' 'path to bed file' 'b'
DEFINE_integer 'bins' '' 'number of bins to divide region' 'n'
DEFINE_string 'factor' '' '[optional] factor to generate files for only that factor' 'f'
DEFINE_string 'fixedBinSizeUpstream' '' "[optional] can be used to create fixed bins upstream; format: 'binsize:binnumber'" ''
DEFINE_string 'fixedBinSizeDownstream' '' "[optional] can be used to create fixed bins downstream; format: 'binsize:binnumber'" ''
DEFINE_boolean 'debug' 'true' '[optional] prints out debug messages.' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "parameters before check" # print param values, if in debug mode

# check if mandatory arguments are there
if [[ -z $FLAGS_outputDir ]]; then
	echo "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_bedgraphTable ]]; then
	echo "Parameter -g must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_bedfile ]]; then
	echoError "Parameter -b must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_bins ]]; then
	echo "Parameter -n must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi


### start script
binner=$SCRIPT_FOLDER"/../../modules/binGenome/binGenome.sh"

bedname=$(basename "$FLAGS_bedfile" | cut -d. -f1)
echo $bedname

IFS=',' read -r -a array <<< $FLAGS_factor
echo "${array[@]}"

while read line; do
	if [[ "$line" == "" ]] ; then #empty line
		continue
	fi
	linearray=($line)
	if [[ "${linearray[0]}" == "BEDGRAPH_POS" ]] ; then #header
		continue
	fi
	if [[ "${linearray[5]}" == "YES" ]] ; then   #use=yes
		if [[ ($FLAGS_factor && (" ${array[@]} " =~ " ${linearray[2]} ")) || -z $FLAGS_factor ]] ; then
		#if [[ ($FLAGS_factor && ("${linearray[2]}" == "$FLAGS_factor")) || -z $FLAGS_factor ]] ; then   #factor given and equal to column or not given
			echo "processing "${linearray[0]}"..."
			newdir=$FLAGS_outputDir"exp"${linearray[6]}"-"${linearray[2]}"/"
			[ -d $newdir ] || mkdir -p $newdir   # create output folder
			filename1=$(basename "${linearray[0]}" .bedgraph)"_anti"
			filename2=$(basename "${linearray[1]}" .bedgraph)"_anti"
			echo $filename1
			if [[ "${linearray[1]}" == "NA" ]] ; then   #bedgraph_neg empty, no strand
				if [[ " ${linearray[3]} " =~  "anti" ]] ; then
					if [[ -z $FLAGS_fixedBinSizeUpstream ]] ; then
						$binner --bedgraph ${linearray[0]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --bedgraphNames $filename1 --bedgraphNames $filename2
					else
						$binner --bedgraph ${linearray[0]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --bedgraphNames $filename1 --bedgraphNames $filename2 --fixedBinSizeUpstream $FLAGS_fixedBinSizeUpstream --fixedBinSizeDownstream $FLAGS_fixedBinSizeDownstream
					fi
				else
					if [[ -z $FLAGS_fixedBinSizeUpstream ]] ; then
                                                $binner --bedgraph ${linearray[0]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir
                                        else
                                                $binner --bedgraph ${linearray[0]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --fixedBinSizeUpstream $FLAGS_fixedBinSizeUpstream --fixedBinSizeDownstream $FLAGS_fixedBinSizeDownstream
                                        fi
				fi
				if [[ -z "$(ls -A $newdir)" ]] ; then
					exit $EXIT_WRITING_FAILED
				fi
			else   #bedgraph_neg given, stranded
				if [[ " ${linearray[3]} " =~  "anti" ]] ; then  #2x das selbe aber mit anti, muss andre namen wählen
					if [[ -z $FLAGS_fixedBinSizeUpstream ]] ; then
                                                $binner --bedgraphPos ${linearray[0]} --bedgraphNeg ${linearray[1]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --bedgraphNames $filename1 --bedgraphNames $filename2
                                        else
                                                $binner --bedgraphPos ${linearray[0]} --bedgraphNeg ${linearray[1]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --bedgraphNames $filename1 --bedgraphNames $filename2 --fixedBinSizeUpstream $FLAGS_fixedBinSizeUpstream --fixedBinSizeDownstream $FLAGS_fixedBinSizeDownstream
                                        fi
				else
					if [[ -z $FLAGS_fixedBinSizeUpstream ]] ; then
                                                $binner --bedgraphPos ${linearray[0]} --bedgraphNeg ${linearray[1]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir
                                        else
                                                $binner --bedgraphPos ${linearray[0]} --bedgraphNeg ${linearray[1]} --annotation $FLAGS_bedfile --bins $FLAGS_bins --outputDir $newdir --fixedBinSizeUpstream $FLAGS_fixedBinSizeUpstream --fixedBinSizeDownstream $FLAGS_fixedBinSizeDownstream
                                        fi
				fi
				if [[ -z "$(ls -A $newdir)" ]] ; then
					exit $EXIT_WRITING_FAILED
				fi
			fi
		fi
	fi
done < ${FLAGS_bedgraphTable}



writeParam2File "$FLAGS_returnFilePath" "coverageFiles" "$FLAGS_outputDir"
writeParam2File "$FLAGS_returnFilePath" "bedname" "$bedname"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK
