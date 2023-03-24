#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi
# define parameters
DEFINE_string 'input' '' 'input file' 'i'
DEFINE_string 'output' '' 'output file' 'o'
DEFINE_float 'probability' '' 'probability of keeping a read (pair)' 'p'
DEFINE_string 'pathToPicard' '' 'path to picard.jar' 'j'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
    echo "1.0"
    exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_input" ]; then
	echoError "Parameter --input must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter --output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_probability" ]; then
	echoError "Parameter --probability must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT #####################################################

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# TODO: write functional part of the script here
# check, if the input files exist
if test -f "$FLAGS_input" ; then
	echo "file exists"
else
	echo "file doesnt exist"
	exit $EXIT_OK
fi
verifyFileExistence "$FLAGS_input"

echo "java -jar picard.jar DownsampleSam I=$FLAGS_input O=$FLAGS_output P=$FLAGS_probability VALIDATION_STRINGENCY=SILENT"
java -jar $FLAGS_pathToPicard/picard.jar DownsampleSam I=$FLAGS_input O=$FLAGS_output P=$FLAGS_probability VALIDATION_STRINGENCY=SILENT
#/home/proj/software/picard/picard-tools-2.18.14
exit $EXIT_OK


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
