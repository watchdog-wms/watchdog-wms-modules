#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='bedtools:echo:wc:rm'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'bam' '' 'path to the position-based-sorted bam file' 'b'
DEFINE_string 'output' '' 'path to bedgraph file' 'o'
DEFINE_string 'contigSizes' '' '[optional] file containing the sizes of the contigs used in the bam file if ranges should be extended (format: <chrName><TAB><SIZE>)' 'c'
DEFINE_string 'returnFilePath' '' '[optional] path to the return variables file' ''
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(bedtools --version 2>&1 | sed -E 's/^bedtools //')
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_bam" ]; then
	echoError "Parameter -b must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# check, if the input files exist
verifyFileExistence "$FLAGS_bam"

if [ ! -z "$FLAGS_contigSizes" ]; then
	verifyFileExistence "$FLAGS_contigSizes"
fi

# verify that tool can write there
OUT_BASE=$(createOutputFolder "$FLAGS_output")
if [ ! -w "$OUT_BASE" ]; then
	echoError "No write permissions in folder '$OUT_BASE'.";
	exit $EXIT_WRITING_FAILED
fi

# run it
if [ ! -z "$FLAGS_contigSizes" ]; then
	MESSAGE=$(bedtools genomecov -split -bg -max 99999999 -ibam "$FLAGS_bam" -g "$FLAGS_contigSizes" > "$FLAGS_output")
else
	MESSAGE=$(bedtools genomecov -split -bg -max 99999999 -ibam "$FLAGS_bam" > "$FLAGS_output")
fi
RET=$?

# check exit code
FAIL=0
if [ $RET -eq 0 ]; then
	# check if bam file is there and contains some lines
	if [ -f "$FLAGS_output" ]; then
		COUNT=$(head -c 512 "$FLAGS_output" | wc -c)
		if [ $COUNT -eq 0 ]; then
			FAIL=1
			echoError "Wig file '$FLAGS_output' is empty!"
		fi
	else
		echoError "Wig file '$FLAGS_output' was not found."
		FAIL=1
	fi
else
	FAIL=1
fi

# output the error
if [ $FAIL -eq 1 ]; then
	echoError "bedtools run failed. Wig file was deleted. See error of bedtools below"
	echoAError "error code: '$RET'"
	 # output the original message
	printf "$MESSAGE\n"
	rm -f "$FLAGS_output" 2>&1 > /dev/null
	exit $EXIT_FAILED
else
	# write return value, if variable is set
	if [ ! -z "$FLAGS_returnFilePath" ]; then
		writeParam2File "$FLAGS_returnFilePath" "wiggleFile" "$FLAGS_output"
		blockUntilFileIsWritten "$FLAGS_returnFilePath"
	fi
	exit $EXIT_OK
fi

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
