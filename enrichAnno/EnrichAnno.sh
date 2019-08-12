#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# check, if used tools are installed
USED_TOOLS='Rscript'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

INSTALL_DIR="$SCRIPT_FOLDER/R/"
SCRIPT_PATH="GO_KEGG.R"

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi
USED_TOOLS='^R|getopt|clusterProfiler|pathview'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'backgroundFile' '' 'path to file with header containing a list of ENSEMBL or GENDCODE identifiers that should be used as backgroud' 'b'
DEFINE_string 'testFiles' '' 'path to file(s) with header containing a list of ENSEMBL or GENDCODE identifiers that should be used for enrichment testing' 't'
DEFINE_string 'orgDB' '' 'name of the orgDB that should be used as GO annotation; if package is missing it is installed via biocLite' 'd'
DEFINE_string 'keggDBName' '' 'organism code for KEGG (e.g. mmu / hsa); http://www.genome.jp/kegg/catalog/org_list.html; if not supported by KEGGREST parameter will be ignored' 'k'
DEFINE_float 'pValueCutoff' '0.01' 'p-Value cutoff for significant results' 'p'
DEFINE_boolean 'plotKegg' 'true' 'if enabled, plots are created for KEGG pathways' ''
DEFINE_string 'output' '' 'path to output basename; folder is created if not existent' 'o'
DEFINE_string 'suffix' '' 'suffix is inserted before basename of output if set; if a absolute path basename is applied' 's'
DEFINE_string 'foldchangeCol' '' 'name of the colum that contains the log2FC' ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_testFiles" ]; then
	echoError "Parameter --testFiles must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_backgroundFile" ]; then
	echoError "Parameter --backgroundFile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_keggDBName" ]; then
	echoError "Parameter --keggDBName must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_orgDB" ]; then
	echoError "Parameter --orgDB must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter --output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
else
	if [ ! -z "$FLAGS_suffix" ]; then
		suffix=$(basename "$FLAGS_suffix")
		BASE=$(basename "$FLAGS_output")
		DIR=$(dirname "$FLAGS_output")
		FLAGS_output="${DIR}/${suffix}/$BASE"
	fi
fi
printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 


# test if files are there and readable
verifyFileExistence "$FLAGS_backgroundFile"
FLAGS_backgroundFile=$(abspath $FLAGS_backgroundFile)

IFS=',' read -ra EXC <<< "$FLAGS_testFiles"
FLAGS_testFiles=""
for I in "${!EXC[@]}"; do 
	F="${EXC[$I]}"
	verifyFileExistence "$F"
	F=$(abspath $F)	
	FLAGS_testFiles="${FLAGS_testFiles}${F},"
done
FLAGS_testFiles=$(echo $FLAGS_testFiles | sed 's#,$##')

# check, if cutoffs have the correct value
if [ $(echo "$FLAGS_pValueCutoff < 0" | bc -l) -eq 1 ] || [ $(echo "$FLAGS_pValueCutoff > 1" | bc -l) -eq 1 ]; then
	echoError "Parameter --pValueCutoff must be between zero and one!"
	exit $EXIT_INVALID_ARGUMENTS
fi

# all parameters seem to be ok, run the script!
RUN2END=$(getTmpFile "DETest")
COMMAND="Rscript '$INSTALL_DIR/$SCRIPT_PATH' --confirmRun2EndFile '$RUN2END'"
# build the command
for PARAM in $__flags_longNames; do
	# ignore that parameter since it is only for the module
	if [ "$PARAM" == "debug" ] || [ "$PARAM" == "suffix" ] || [ "$PARAM" == "help" ]; then
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
echo $COMMAND
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
