#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='cat:sed:touch:rm:ls:sort:printf:find'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoAError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi














 
# define parameters
DEFINE_string 'folder' '' 'one ore more input folders; one for each pattern' 'f'
DEFINE_string 'output' '' 'write results to a file; one line per found file' 'o'
DEFINE_string 'sep' ',' 'separator between entries; default: ,' 's'
DEFINE_string 'pattern' '' 'one ore more unix file pattern (e.g. *.txt) that are used to find files matching that pattern; one pattern corresponds to one input folder path' 'p'
DEFINE_integer 'maxdepth' '' 'descend at most n levels of folders; default: 0' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_folder" ]; then
	echoError "Parameter -f must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_pattern" ]; then
	echoError "Parameter -p must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
FLAGS_maxdepth=$(ensureLowerBound "$FLAGS_maxdepth" 0)

# split the stuff
IFS=',' read -ra IN <<< "$FLAGS_folder"
IFS=',' read -ra PATTERNS <<< "$FLAGS_pattern"
COUNT=0

if [ "${#PATTERNS[@]}" -eq 1 ]; then
	SINGLE_PATTERN=1
else
	if [ "${#PATTERNS[@]}" -ne "${#IN[@]}" ]; then
		echoError "Only one or the same number of patterns than input folders must be given!"
		exit $EXIT_MISFORMATED_INPUT
	fi
	SINGLE_PATTERN=0
fi

IFS=$'\n'
# process all pattern
for I in "${!IN[@]}"; do 
	P="${IN[$I]}"

	# verify that folder is there
	verifyFolderExistence "$P"
	if [ "$SINGLE_PATTERN" -eq 1 ]; then
		PAT="${PATTERNS[0]}"
	else
		PAT="${PATTERNS[$I]}"
	fi

	# get all files which match that pattern
	RET=$(find "$P" -maxdepth ${FLAGS_maxdepth} -name "$PAT" 2> /dev/null)
	if [ $? -eq 0 ]; then
		for FILE in $RET; do 
			FILE=$(readlink -e "$FILE")
			verifyFileExistence "$FILE" # test if file is there and readable
			FILES[$COUNT]="$FILE"
			COUNT=$((COUNT+1))
		done
	else
		echoError "Pattern '$PAT' in '$P' did not match any files."
		exit $EXIT_MISFORMATED_INPUT
	fi
done

# sort the results
FILES_SORTED=($(printf "%s\n" ${FILES[@]} | sort ))

# create output folder
if [ ! -z "$FLAGS_output" ]; then
	OUT_BASE_FOLDER=$(createOutputFolder "$FLAGS_output")
	printf "%s\n" ${FILES_SORTED[@]} > "$FLAGS_output"
fi

# joins results
RES=$(printf "%s${FLAGS_sep}" ${FILES_SORTED[@]})

writeParam2File "$FLAGS_returnFilePath" "foundFiles" "$RES"
blockUntilFileIsWritten "$FLAGS_returnFilePath"
exit $EXIT_OK
