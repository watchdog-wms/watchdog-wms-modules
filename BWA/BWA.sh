#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='bwa:grep'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoNew "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'index' '' 'genome index' 'i'
DEFINE_string 'in1' '' 'in1.fastq' ''
DEFINE_string 'in2' '' '[optional] in2.fastq' ''
DEFINE_string 'out' '' 'outfile' 'o'
DEFINE_boolean 'all' '1' 'output all alignments for SE or unpaired PE [default false]' 'a'
DEFINE_boolean 'ignoreIndexExistence' '1' 'do not throw an error if index does not exist [default false]' 'e'
DEFINE_integer 'numberOfThreads' '1' '[optional] number of threads [default 1]' 'n'
DEFINE_integer 'minimumSeedLength' '19' '[optional] minimum seed length [default 19]' 'm'
DEFINE_boolean 'version' '1' "print version information and quit" ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
        
        MESSAGE=$(bwa 2>&1 | grep 'Version' | cut -d " " -f 2);

        echo $MESSAGE
        exit $EXIT_OK
fi


INDEX_FILE="$FLAGS_index.amb"
if [ $FLAGS_ignoreIndexExistence -eq 0 ]; then
    
    if [ ! -f "$INDEX_FILE" ]; then
        printf "$FLAGS_index does not exist. Ignore and exit!\n"
        exit $EXIT_OK
    fi

fi

verifyFileExistence "$INDEX_FILE"


# check if mandatory arguments are there
if [ -z "$FLAGS_index" ]; then
	echoError "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_in1" ]; then
	echoError "Parameter in1 must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_out" ]; then
	echoError "Parameter out must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_numberOfThreads" -lt 1 ] ; then
	echoError "Parameter -n must be greater than 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi


# create output directory if not already exists
dirname_out=$(dirname "$FLAGS_out")
if [ ! -d "$dirname_out" ]; then
    mkdir -p "$dirname_out"
fi


printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

COMMAND="bwa mem -t $FLAGS_numberOfThreads -k $FLAGS_minimumSeedLength";

if [ "$FLAGS_all" -eq 0 ]; then
    COMMAND="$COMMAND -a"
fi

COMMAND="$COMMAND $FLAGS_index $FLAGS_in1";

if [ -n "$FLAGS_in2" ]; then
    COMMAND="$COMMAND $FLAGS_in2"
fi

COMMAND="$COMMAND | grep -v 'AS:i:0'> $FLAGS_out";

MESSAGE=$(eval $COMMAND 2>&1)
RET=$?

echo $MESSAGE;
exit $RET

# check exit code
if [ $RET -eq 0 ]; then
        # output the original message
        printf "$MESSAGE\n"
        exit $EXIT_OK
else
        echoError "BWA run failed. Output file was deleted. See log of BWA below"
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
