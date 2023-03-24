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
DEFINE_string 'output' '' 'output directory' 'o'
DEFINE_string 'name' '' 'name' 'n'
DEFINE_integer 'd1' '10000' '[optional] maximum distance to gene end' 'd'
DEFINE_integer 'd2' '5000' '[optional] maximum distance to last added dOCR' 'p'
DEFINE_string 'annotation' '' 'genome annotation file' 'a'
DEFINE_string 'gene' '' '[optional] get in_gene_length and fraction coverage, default false' 'g'

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_input" ]; then
	echoError "Parameter -input must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter -output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_name" ]; then
	echoError "Parameter -name must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_annotation" ]; then
	echoError "Parameter -annotation must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT #####################################################

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# check, if the input files exist
if test -f "$FLAGS_input" ; then
	echo "file exists"
fi
verifyFileExistence "$FLAGS_input"


if [ ! -d "$FLAGS_output" ]; then
  mkdir $FLAGS_output;
	echo "made dir"
fi

echo "$FLAGS_gene"
if [ -z "$FLAGS_gene" ]; then
	echo "gene is not set";
	java -cp $SCRIPT_FOLDER/dist/dOCRTools.jar de.lmu.ifi.bio.docrtools.ATACSummarizer -b $FLAGS_input -o $FLAGS_output -name $FLAGS_name  -d1 $FLAGS_d1 -d2 $FLAGS_d2 -a $FLAGS_annotation	
else
	echo "gene is set";
	java -cp $SCRIPT_FOLDER/dist/dOCRTools.jar de.lmu.ifi.bio.docrtools.ATACSummarizer -b $FLAGS_input -o $FLAGS_output -name $FLAGS_name  -d1 $FLAGS_d1 -d2 $FLAGS_d2 -a $FLAGS_annotation -gene


fi



exit $EXIT_OK


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
