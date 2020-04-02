#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='wget:du:cut:sed'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'uri' '' "URI pointing to the resource to download; multiple URIs muste be separated by ','" 'u'
DEFINE_string 'output' '' 'path to a folder in which the downloaded files should be stored; filename remains untouched' 'o'
DEFINE_string 'rename' '' "renames the file to that name; multiple names must be provided in case of multiple URIs (separated by ',')" 'r'
DEFINE_boolean 'disableSizeCheck' '1' '[optional] flag that can be used to disable the size check that checks if a file is greater than 1KB' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE=$(wget --version 2>&1 | head -n 1 | sed -E 's/built on.*//')
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_uri" ]; then
	echoError "Parameter -u must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
# ensure that the parent folder is there
OUT_FOLDER=$(createOutputFolder "$FLAGS_output/x")

# download the files
IFS=',' read -ra URIS <<< "$FLAGS_uri"
IFS=',' read -ra RENAMES <<< "$FLAGS_rename"

# ensure that rename is either empty or has the same number of names as URIs are given
if [ ! -z "${FLAGS_rename}" ] && [ ${#URIS[@]} -ne ${#RENAMES[@]} ]; then
	echoError "The number of --uri parameters (${#URIS[@]}) must be the same as the number of given --rename names (${#RENAMES[@]})."
	exit $EXIT_INVALID_ARGUMENTS	
fi

for I in "${!URIS[@]}"; do 
	URI="${URIS[$I]}"
	BASE=$(echo "$URI" | sed -E 's|.+/||g')
	if [ ! -z "${FLAGS_rename}" ]; then
		RENAME="${RENAMES[$I]}"
		BASE=$RENAME
	fi
	executeCommand "wget --no-verbose -O \"$FLAGS_output/$BASE\" \"$URI\"" "/dev/null" "wget"
	# ensure that the new file is there
	verifyFileExistence "$FLAGS_output/$BASE"

	# ensure that the file has a file size greater than 1KB
	if [ $FLAGS_disableSizeCheck -eq 0 ]; then
		SIZE=$(du -B 1 "$FLAGS_output/$BASE" | cut -f 1 )
		if [ $SIZE -lt 1024 ]; then
			echoError "Size of file '$FLAGS_output/$BASE' seems to be too small ($SIZE bytes). Disable this check with --disableSizeCheck"
			exit $EXIT_FAILED
		fi
	fi
	FILES_all="${FILES_all},${FLAGS_output}/${BASE}"
done
FILES_all=$(echo "${FILES_all}" | sed -E 's/^,//') # remove trailing ','

# check, if return parameters must be set
if [ ! -z "$FLAGS_returnFilePath" ]; then
	writeParam2File "$FLAGS_returnFilePath" "downloadedFiles" "$FILES_all" # multiple files are separated by ','
	writeParam2File "$FLAGS_returnFilePath" "downloadedFolder" "$FLAGS_output"
	writeParam2File "$FLAGS_returnFilePath" "numberOfFiles" "${#URIS[@]}"
	blockUntilFileIsWritten "$FLAGS_returnFilePath"
fi

# all seems to be ok
exit $EXIT_OK
