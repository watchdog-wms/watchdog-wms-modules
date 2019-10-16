#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='samtools'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'infile' '' 'input bam file(s)' 'i'
DEFINE_string 'outfile' '' 'output bam file' 'o'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(samtools 2>&1 | head -n 3 | tail -n 1)
	echo $MESSAGE
	exit $EXIT_OK
fi


# check if mandatory arguments are there
if [ -z "$FLAGS_outfile" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_infile" ]; then
	echoError "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #



COMMAND="samtools merge $FLAGS_outfile"


IFS=',' read -ra FILES <<< "$FLAGS_infile"
for I in "${!FILES[@]}"; do 
        FILE="${FILES[$I]}"  
        verifyFileExistence "$FILE"
        COMMAND="$COMMAND $FILE"
done

MESSAGE=$(eval $COMMAND 2>&1)
RET=$?

# check exit code
if [ $RET -eq 0 ]; then
        # output the original message
        printf "$MESSAGE\n" 1>&2
        exit $EXIT_OK
else
        echoError "mergeBAM run failed. Output file was deleted. See log of mergeBAM below"
        echoAError "error code: '$RET'"
         # output the original message
        printf "$MESSAGE\n"
        rm -f "$FLAGS_out" 2>&1 > /dev/null
        exit $EXIT_FAILED
fi



# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
