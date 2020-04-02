#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='which:python:DaPars_main.py:echo:wc:rm:cat:grep:sed:tr'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'controlCondition' '' 'name of the control condition' ''
DEFINE_string 'testCondition' '' 'name of the test condition' ''
DEFINE_string 'sampleAnnotation' '' 'annotation file with sample names in the first colum and sample condition in the second condition (header: sample\tcondition)' 's'
DEFINE_string 'excludeSamples' '' "names of samples that should be excluded from the analysis; separated by ','" 'e'
DEFINE_string 'wigFolder' '' 'folder containing the wig files (format: folder/samplename.bedgraph)' 'w'
DEFINE_string 'wigEnding' 'bedgraph' '[otpional] ending of the wig files' ''
DEFINE_string 'annotated3UTR' '' 'path to annotated 3 UTR regions created with DaPars_Extract_Anno.py' 'a'
DEFINE_string 'outputFile' '' 'path to the output file' 'o'
DEFINE_integer 'coverageCutoff' '30' '[optional] coverage threshold' ''
DEFINE_float 'FDRCutoff' '0.01' '[optional] FDR cutoff' ''
DEFINE_float 'PDUICutoff' '0.5' '[optional] degree of difference in APA usage ([0,1])' ''
DEFINE_float 'FoldChangeCutoff' '0.5' '[optional] log2 foldchange between the two conditions' ''
DEFINE_integer 'numberOfCondASamplesReachingCutoff' '' '[optional] number of samples from condition A that must pass the coverage cutoff; default: all samples' ''
DEFINE_integer 'numberOfCondBSamplesReachingCutoff' '' '[optional] number of samples from condition B that must pass the coverage cutoff; default: all samples' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	PTOOL=$(which DaPars_main.py)
	MESSAGE=$(grep -A 1 "def get_version():" "$PTOOL" | tail -n 1 | sed 's/return "//' | sed 's/"//')
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_sampleAnnotation" ]; then
	echoError "Parameter --sampleAnnotation must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_controlCondition" ]; then
	echoError "Parameter --controlCondition must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_testCondition" ]; then
	echoError "Parameter --testCondition must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_wigFolder" ]; then
	echoError "Parameter --wigFolder must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_annotated3UTR" ]; then
	echoError "Parameter -a must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_outputFile" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# check, if the input files exist
verifyFileExistence "$FLAGS_annotated3UTR"

# check, if name of control and test condition are in sample file
NUM_CONTROL=$(cut -f 2 "$FLAGS_sampleAnnotation" | grep -c "$FLAGS_controlCondition")
NUM_TEST=$(cut -f 2 "$FLAGS_sampleAnnotation" | grep -c "$FLAGS_testCondition")
if [ "$NUM_CONTROL" -eq 0 ] || [ "$NUM_TEST" -eq 0 ]; then
	echoError "Both test conditions must contain at least one sample ($FLAGS_controlCondition: $NUM_CONTROL / $FLAGS_testCondition: $NUM_TEST).";
	exit $EXIT_INVALID_ARGUMENTS
fi

# check, if the samples, which should be excluded are part of the annotation file
TMP_FILE=$(getTmpFile "DaPars")
cp "$FLAGS_sampleAnnotation" "$TMP_FILE"
IFS=',' read -ra EXC <<< "$FLAGS_excludeSamples"
for I in "${!EXC[@]}"; do 
	E="${EXC[$I]}"
	COUNT=$(cut -f 1 "$TMP_FILE" | grep "$E" -c)
	if [ $COUNT -eq 1 ]; then
		# get the type of the sample
		T=$(grep -E "^$E" "$TMP_FILE" | cut -f 2)
		if [ "$T" == "$FLAGS_controlCondition" ]; then
			NUM_CONTROL=$((NUM_CONTROL-1))
		else
			if [ "$T" == "$FLAGS_testCondition" ]; then
				NUM_TEST=$((NUM_TEST-1))
			fi
		fi
		# delete that line
		LINE=$(grep -n -E "^$E" "$TMP_FILE" | cut -f 1 -d ":")
		sedinline "${LINE}d" "$TMP_FILE"
	else
		echoError "Sample '$E' can not be excluded because it is not part of the sample list.";
		exit $EXIT_INVALID_ARGUMENTS
	fi
done

# check, if after removal enough samples are there
if [ "$NUM_CONTROL" -le 0 ] || [ "$NUM_TEST" -le 0 ]; then
	echoError "After removal of excluded samples at least one test condition contains no samples anymore ($FLAGS_controlCondition: $NUM_CONTROL / $FLAGS_testCondition: $NUM_TEST).";
	exit $EXIT_INVALID_ARGUMENTS
fi

# get the samples that are left
CONTROL=$(grep -E "\t$FLAGS_controlCondition" "$TMP_FILE" | cut -f 1 | tr '\n' ',')
TEST=$(grep -E "\t$FLAGS_testCondition" "$TMP_FILE" | cut -f 1 | tr '\n' ',')
IFS=',' read -ra EXC <<< "$CONTROL"

for I in "${!EXC[@]}"; do 
	F="${EXC[$I]}"
	FNAME="$FLAGS_wigFolder/$F.$FLAGS_wigEnding"
	verifyFileExistence "$FNAME"
	wigConditionA="$FNAME,$wigConditionA"
done
wigConditionA=${wigConditionA%?}
if [ -z "$FLAGS_numberOfCondASamplesReachingCutoff" ]; then
	FLAGS_numberOfCondASamplesReachingCutoff=$((I+1))
fi
IFS=',' read -ra EXC <<< "$TEST"
for I in "${!EXC[@]}"; do 
	F="${EXC[$I]}"
	FNAME="$FLAGS_wigFolder/$F.$FLAGS_wigEnding"
	verifyFileExistence "$FNAME"
	wigConditionB="$FNAME,$wigConditionB"
done
wigConditionB=${wigConditionB%?}
if [ -z "$FLAGS_numberOfCondBSamplesReachingCutoff" ]; then
	FLAGS_numberOfCondBSamplesReachingCutoff=$((I+1))
fi
# delete tmp file
rm -f "$TMP_FILE" 2>&1 > /dev/null

# verify that tool can write there
OUT_BASE=$(createOutputFolder "$FLAGS_outputFile/dummy")
if [ ! -w "$OUT_BASE" ]; then
	echoError "No write permissions in folder '$OUT_BASE'.";
	exit $EXIT_WRITING_FAILED
fi

# ensure that optional parameters have the correct ranges
if [ "$FLAGS_coverageCutoff" -le 0 ]; then
	echoError "Parameter --coverageCutoff must be greater than 0."
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_numberOfCondASamplesReachingCutoff" -le 0 ]; then
	echoError "Parameter --numberOfCondASamplesReachingCutoff must be greater than 0."
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_numberOfCondBSamplesReachingCutoff" -le 0 ]; then
	echoError "Parameter --numberOfCondBSamplesReachingCutoff must be greater than 0."
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ $(echo "$FLAGS_FDRCutoff < 0" | bc -l) -eq 1 ] || [ $(echo "$FLAGS_FDRCutoff > 1" | bc -l) -eq 1 ]; then
	echoError "Parameter --FDRCutoff must be between zero and one!"
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ $(echo "$FLAGS_PDUICutoff < 0" | bc -l) -eq 1 ]; then
	echoError "Value of PDUICutoff must be greater than 0."
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ $(echo "$FLAGS_FoldChangeCutoff < 0" | bc -l) -eq 1 ]; then
	echoError "Value of FoldChangeCutoff must be greater than 0."
	exit $EXIT_INVALID_ARGUMENTS
fi

# all parameter seem to be ok --> write config file in output folder
CONFIG_FILE="$FLAGS_outputFile.dapars.config"
OUT_BASE=$(dirname "$FLAGS_outputFile")
BASENAME=$(basename "$FLAGS_outputFile")
cat >"$CONFIG_FILE" <<EOL
#DaPars config file documented here: http://lilab.research.bcm.edu/dldcc-web/lilab/zheng/DaPars_Documentation/html/DaPars.html
#The following file is the result of step 1.
Annotated_3UTR=${FLAGS_annotated3UTR}

#A comma-separated list of BedGraph files of samples from condition A and B
Group1_Tophat_aligned_Wig=${wigConditionA}
Group2_Tophat_aligned_Wig=${wigConditionB}

# output folder and file
Output_directory=${OUT_BASE}
Output_result_file=${BASENAME}

#At least how many samples passing the coverage threshold in two conditions
Num_least_in_group1=${FLAGS_numberOfCondASamplesReachingCutoff}
Num_least_in_group2=${FLAGS_numberOfCondBSamplesReachingCutoff}
Coverage_cutoff=${FLAGS_coverageCutoff}

#Cutoff for FDR of P-values from Fisher exact test.
FDR_cutoff=${FLAGS_FDRCutoff}
PDUI_cutoff=${FLAGS_PDUICutoff}
Fold_change_cutoff=${FLAGS_FoldChangeCutoff}
EOL

# run it
cd /tmp
PTOOL=$(which DaPars_main.py)
MESSAGE=$(python "$PTOOL" "$CONFIG_FILE")
RET=$?

# check exit code
FAIL=0
if [ $RET -eq 0 ]; then
	# check if result file is there and contains some lines
	if [ -f "${FLAGS_outputFile}_All_Prediction_Results.txt" ]; then
		mv "${FLAGS_outputFile}_All_Prediction_Results.txt" "${FLAGS_outputFile}"
		COUNT=$(head -l 10 "$FLAGS_outputFile" | wc -c)
		if [ $COUNT -eq 0 ]; then
			FAIL=1
			echoError "Return file '$FLAGS_outputFile' is empty!"
		fi
	else
		echoError "Return file '$FLAGS_outputFile' was not found."
		FAIL=1
	fi
else
	FAIL=1
fi
exit
# output the error
if [ $FAIL -eq 1 ]; then
	echoError "DaPars run failed. Return file was deleted. See error of DaPars below"
	echoAError "error code: '$RET'"
	 # output the original message
	printf "$MESSAGE\n"
	rm -f "$FLAGS_output" 2>&1 > /dev/null
	rm -f "$CONFIG_FILE" 2>&1 > /dev/null
	exit $EXIT_FAILED
else
	exit $EXIT_OK
fi

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
