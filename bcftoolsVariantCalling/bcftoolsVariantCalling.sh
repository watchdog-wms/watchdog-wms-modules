#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='bcftools'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'reference' '' 'Path to the file containing the reference genome.' 'r'
DEFINE_string 'bamfile' '' 'Path to the input bam file.' 'b'
DEFINE_string 'vcf' '' 'Path of the output vcf file.' 'v'
DEFINE_integer 'maxdepth' '100000' 'Maximum number of reads per position.' 'm'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'i'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(bcftools -v 2>&1)
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_vcf" ]; then
	echoError "Parameter -v must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_bamfile" ]; then
	echoError "Parameter -b must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_reference" ]; then
	echoError "Parameter -r must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_maxdepth" -gt 1000000000 ] || [ "$FLAGS_maxdepth" -lt 100 ]; then
	echoError "Parameter -m must be between [100, 1000000000]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi

# TODO: add further parameter checks here which refer to input ranges or valid parameters in general here

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

command_to_call="bcftools mpileup -Ou -f $FLAGS_reference -d $FLAGS_maxdepth"

dirname=$(dirname "$FLAGS_vcf")
if [ ! -d "$dirname" ]; then
    mkdir -p "$dirname"
fi

command_to_call="$command_to_call $FLAGS_bamfile | bcftools call -mv -Ob -o $FLAGS_vcf"

MESSAGE=$(eval $command_to_call 2>&1)
RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile bcftoolsVariantCalling)
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
