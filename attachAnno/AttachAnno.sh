#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='Rscript:tr:wc'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

INSTALL_DIR="$SCRIPT_FOLDER/R/"
SCRIPT_PATH="AttachAnno.R"

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi
USED_TOOLS='^R|getopt|'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'targetFile' '' 'path to char-separated table file with header' 't'
DEFINE_string 'outputFile' '' 'path to the annotated output file' 'o'
DEFINE_string 'targetIDcolumn' '' 'name of the column of the target file that should be used to merge the table with the annotation file(s)' ''
DEFINE_string 'annotationFile' '' '1:n path(s) to annotation table file(s) that should be attached; n values expected; separated by ,' ''
DEFINE_string 'annotationIDcolumn' '' 'name(s) of the column(s) of the annotation file(s) that should be used to merge the table with the annotation file(s); 1 or n values expected; separated by ,' ''
DEFINE_string 'annotationSep' '\t' '[optional] separating char in the target file; 1 or n values expected; separated by ,; default: TAB' ''
DEFINE_string 'targetSep' '\t' '[optional] separating char in the annotation file(s); default: TAB' ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_targetFile" ]; then
	echoError "Parameter --targetFile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_targetIDcolumn" ]; then
	echoError "Parameter --targetIDcolumn must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_annotationFile" ]; then
	echoError "Parameter --annotationFile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_annotationIDcolumn" ]; then
	echoError "Parameter --annotationIDcolumn must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_outputFile" ]; then
	echoError "Parameter --outputFile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 


# test if files are there and readable
verifyFileExistence "$FLAGS_targetFile"

IFS=',' read -ra EXC <<< "$FLAGS_annotationFile"
FLAGS_annotationFile=""
N=0
for I in "${!EXC[@]}"; do 
	F="${EXC[$I]}"
	verifyFileExistence "$F"
	F=$(abspath $F)	
	FLAGS_annotationFile="${FLAGS_testFiles}${F},"
	N=$((N+1))
done
FLAGS_annotationFile=$(echo $FLAGS_annotationFile | sed 's#,$##')

# check, if annotation file is set
if [ $N -eq 0 ]; then
	echoError "Parameter --annotationFile must contain at least one absolute file path!"
	exit $EXIT_INVALID_ARGUMENTS
fi
# ensure that a correct number of additional arguments is given
N1=$(echo "$FLAGS_annotationIDcolumn" | tr -cd ',' | wc -c)
N2=$(echo "$FLAGS_annotationSep" | tr -cd ',' | wc -c)
C1=$(echo "$FLAGS_annotationIDcolumn" | wc -c)
C2=$(echo "$FLAGS_annotationSep" | wc -c)

if [ $N -ne $((N1+1)) ] && [ $N1 -ne 0 ] && [ $C1 -ne 0 ]; then
	echoError "Parameter --annotationIDcolumn is expected either one or n times!"
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ $N -ne $((N2+1)) ] && [ $N2 -ne 0 ] && [ $C2 -ne 0 ]; then
	echoError "Parameter --annotationSep is expected either one or n times!"
	exit $EXIT_INVALID_ARGUMENTS
fi

# all parameters seem to be ok, run the script!
RUN2END=$(getTmpFile "AttachAnnotation")
COMMAND="Rscript '$INSTALL_DIR/$SCRIPT_PATH' --confirmRun2EndFile '$RUN2END'"
# build the command
for PARAM in $__flags_longNames; do
	# ignore that parameter since it is only for the module
	if [ "$PARAM" == "debug" ] || [ "$PARAM" == "help" ]; then
		continue
	fi

	V_NAME='FLAGS_'"$PARAM"
	VALUE=${!V_NAME}
	V_NAME='__flags_'"$PARAM"'_type'
	TYPE=${!V_NAME}
	
	# test if value for that parameter was set
	if [ ! -z "$VALUE" ]; then
		# test if flag
		if [ $TYPE -eq 1 ]; then
			# test if flag is enabled
			if  [ $VALUE -eq 0 ]; then
				COMMAND="$COMMAND --$PARAM 1"
			else
				COMMAND="$COMMAND --$PARAM 0"
			fi
		else 
			COMMAND="$COMMAND --$PARAM '$VALUE'"
		fi
	fi
done
# execute the command
BASEDIR=$(dirname "$FLAGS_output")
mkdir -p "$BASEDIR"
LOG="$FLAGS_output.log"
if [ -e "$LOG" ]; then
	rm "$LOG"
fi
executeCommand "$COMMAND" "$LOG" "Rscript DETest"
# test, if the script run to its end
verifyRun2End "$RUN2END"

# if, we come till here, all should be ok.
exit $EXIT_OK
