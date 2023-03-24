#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi
# define parameters
DEFINE_string 'idxstats' '' 'idxstats file' 'i'
echo "idx"
DEFINE_string 'output' '' 'output table file' 'o'
DEFINE_string 'exclude' '' '[optional] chromosomes to be excluded, comma separated' 'e'
DEFINE_string 'samples' '' '[optional] samples to be used, comma separated' 's'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'

echo "parse param"
# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
    echo "1.0"
    exit $EXIT_OK
fi

# check if mandatory arguments are there
echo "check param"
if [ -z "$FLAGS_idxstats" ]; then
	echoError "Parameter --idxstats must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
else
	echo "set correct"
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter --output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT #####################################################

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# TODO: write functional part of the script here
# check, if the input files exist
echo "start script"
echo $FLAGS_idxstats
if test -f "$FLAGS_idxstats" ; then
	echo "file exists"
fi
verifyFileExistence "$FLAGS_idxstats"

if [ -n "$FLAGS_exclude" ]; then
	echo "exlusion set"
	excl=$FLAGS_exclude
else
	excl=''
	echo "exclusion not set"
	#check all chromosomes
fi

if [ -n "$FLAGS_samples" ]; then
	echo "samples set"
	samp=$FLAGS_samples
	#only check those samples
else
	echo "samples not set"
	samp='all'
fi
 echo "python3 $SCRIPT_FOLDER/calcDownsample.py $FLAGS_idxstats $excl $samp $FLAGS_output"
python3 $SCRIPT_FOLDER/calcDownsample.py $FLAGS_idxstats $excl $samp $FLAGS_output


echo "done"
echo $FLAGS_idxstats
exit $EXIT_OK
echo "exit ok"

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
