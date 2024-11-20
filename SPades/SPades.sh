#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='python3'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'forward' '' 'Path to FastQ or FastQ.gz file containing the forward reads' 'f'
DEFINE_string 'reverse' '' 'Path to FastQ or FastQ.gz file containing the reverse reads' 'r'
DEFINE_string 'cons_path' '' 'Path to file containing consensus sequences (from svCaller)' 'c'
DEFINE_string 'outFolder' '' 'Path to output folder, where SPades stores all its resulting files.' 'o'
DEFINE_integer 'memory' '40' '[optional] RAM limit.' 'm'
DEFINE_boolean 'ignoreConsensusExistence' '1' 'do not throw an error if file containing consensus sequences does not exist' 'e'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(rnaspades.py -v 2>&1)
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_outFolder" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_cons_path" ]; then
	echoError "Parameter -c must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_reverse" ]; then
	echoError "Parameter -r must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_forward" ]; then
	echoError "Parameter -f must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_memory" -lt 1 ] ; then
	echoError "Parameter -m must be greater than 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi

# TODO: add further parameter checks here which refer to input ranges or valid parameters in general here

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

if [ $FLAGS_ignoreConsensusExistence -eq 0 ]; then
    
    if [ ! -f "$FLAGS_cons_path" ]; then
        printf "$FLAGS_cons_path does not exist. Ignore and exit!\n"
        exit $EXIT_OK
    fi

fi

mkdir -p "$FLAGS_outFolder"

command_to_call_out="rnaspades.py -1 $FLAGS_forward -2 $FLAGS_reverse -o $FLAGS_outFolder -m $FLAGS_memory"
MESSAGE=$(eval $command_to_call_out 2>&1)
RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile SPades)
touch "$TMP_FILE"
trap "rm -f \"$TMP_FILE\" 2>&1 > /dev/null" EXIT
printf "$MESSAGE" > "$TMP_FILE"
MESSAGE_ERROR=$($SCRIPT_FOLDER/../../core_lib/errorChecker.sh "$TMP_FILE" "truncated file:fail to open file" 2>&1)
CODE_ERROR=$?
rm -f "$TMP_FILE" 2>&1 > /dev/null

if [ $CODE_ERROR -ne 0 ]; then
	echoError "Error checker found some errors, see found errors below"
	echo -e "$MESSAGE_ERROR"
	exit $EXIT_FAILED
else
	# check exit code
	if [ $FAIL -eq 0 ] && [ $RET -eq 0 ]; then
		# output the original message
		printf "$MESSAGE\n"
		
		exit $EXIT_OK
	else
		FAIL=1
	fi
	if [ $FAIL -eq 1 ]; then
		echoError "Run failed. See error below"
		echoAError "error code: '$RET'"
		 # output the original message
		printf "$MESSAGE\n"
		exit $EXIT_FAILED
	fi
fi

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT

