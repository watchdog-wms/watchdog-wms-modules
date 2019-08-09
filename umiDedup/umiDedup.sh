#!/bin/bash
# todo: add support for new illumina header if needed
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='umi_tools:rm:file:tail:grep'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'bamFile' '' "path to the BAM file; UMI must be a suffix of the fastq id separated with '_'" 'i'
DEFINE_string 'outputFile' '' 'path to the de-duplicated BAM file' 'o'
DEFINE_boolean 'deleteOnSuccess' '1' '[optional] deletes the BAM file when deduplication was successfull' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_bamFile" ]; then
	echoError "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
else
	# check, if the input files exist
	verifyFileExistence "$FLAGS_bamFile"

	# check if it might be a bam file
	if [ $(file "$FLAGS_bamFile" | grep -c "gzip compressed") -ne 1 ]; then
		echoError "Input file '$FLAGS_bamFile' seems to be not a file in BAM format.";
		exit $EXIT_INVALID_ARGUMENTS
	fi
fi
if [ -z "$FLAGS_outputFile" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
# ensure that the parent folder is there
OUT_FOLDER=$(createOutputFolder "$FLAGS_outputFile")

# verify that tool can write there
if [ ! -w "$OUT_FOLDER" ]; then
	echoError "No write permissions in folder '$OUT_FOLDER'.";
	exit $EXIT_WRITING_FAILED
fi

# run it
MESSAGE=$(umi_tools dedup -I "$FLAGS_bamFile"  -S "$FLAGS_outputFile" 2>&1)
RET=$?

# check exit code and message
FAIL=1
if [ $RET -eq 0 ]; then
	# ensure log file also says that all is correct
	if [ $(printf "$MESSAGE\n" |tail -n 1 | grep -c "^# job finished in ") -eq 1 ]; then
		FAIL=0
	else
		echoError "umi_tools dedup run failed. Log file does not end with expected last line."
	fi
fi

# all was ok
if [ $RET -eq 0 ] && [ $FAIL -eq 0 ]; then
	if [ "$FLAGS_deleteOnSuccess" -eq 0 ]; then
		rm -f "$FLAGS_bamFile" 2>&1 > /dev/null
	fi

	# check, if return parameters must be set
	if [ ! -z "$FLAGS_returnFilePath" ]; then
		writeParam2File "$FLAGS_returnFilePath" "deduplicatedFile" "$FLAGS_outputFile"
		blockUntilFileIsWritten "$FLAGS_returnFilePath"
	fi

	# output the original message
	printf "$MESSAGE\n"
	exit $EXIT_OK
else
	echoError "umi_tools dedup run failed. Output file was deleted. See error of umi_tools below"
	echoAError "error code: '$RET'"
	 # output the original message
	printf "$MESSAGE\n"
	rm -f "$FLAGS_outputFile" 2>&1 > /dev/null
fi

# all seems to be ok
exit $EXIT_OK
