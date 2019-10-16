#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='tar'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'infile' '' 'input file, must be *tar, *tar.gz or *tar.bz2' 'i'
DEFINE_string 'outputdir' '' '[optional] output directory for extracting archive' 'o'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(tar --version | head -n 1 | awk '{print $4}')
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_infile" ]; then
	echoError "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

if [[ ! $FLAGS_infile == *.tar ]] && [[ ! $FLAGS_infile == *.tar.bz2 ]] && [[ ! $FLAGS_infile == *.tar.gz ]]; then
	 echo "Input file must be either *.tar, *.tar.gz or *tar.bz2"
	 exit $EXIT_INVALID_ARGUMENTS
fi


printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# check, if the input files exist
verifyFileExistence "$FLAGS_infile"

COMMAND="tar";
if [[ $FLAGS_infile == *.tar ]] ; then
	 COMMAND="$COMMAND xfv";
fi
if [[ $FLAGS_infile == *.tar.gz ]] ; then
	 COMMAND="$COMMAND xfvz";
fi
if [[ $FLAGS_infile == *.tar.bz2 ]] ; then
	 COMMAND="$COMMAND xfvj";
fi

COMMAND="$COMMAND '$FLAGS_infile'";

if [ ! -z "$FLAGS_outputdir" ]; then
	if [ ! -d "$FLAGS_outputdir" ]; then
		MESSAGE=$(mkdir "$FLAGS_outputdir")
        CODE=$?
        if [ $CODE -ne 0 ]; then
			echoError "Creation of output directory failed. See error below"
			echoAError "$MESSAGE, error code: '$CODE'"
			exit $EXIT_FAILED
		fi
	fi
	COMMAND="$COMMAND -C '$FLAGS_outputdir'";
fi


MESSAGE=$(eval $COMMAND)
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "Extraction of '$FLAGS_infile' failed."
	echoError "$MESSAGE, error code: '$CODE'"
	exit $EXIT_FAILED
else
	exit $EXIT_OK
fi


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
