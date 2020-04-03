#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='java:echo:printf:rm:grep:head:wc:nproc'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'jarPath' '' 'path to GEM jar file' ''
DEFINE_string 'readDistribution' '' "read spatial distribution file" ''
DEFINE_string 'expt' '' "aligned read file" ''
DEFINE_boolean 'gpsOnly' 'true' "run in GPS only mode" ''
DEFINE_integer 'k' '8' "length of the k-mer for motif finding, use --k or (--kmin & --kmax); GEM parameter" ''
DEFINE_integer 'kMin' '6' "min value of k, e.g. 6; GEM parameter" ''
DEFINE_integer 'kMax' '13' "max value of k, e.g. 13; GEM parameter" ''
DEFINE_string 'seed' '' "exact k-mer string to jump start k-mer set motif discovery; GEM parameter" ''
DEFINE_string 'genome' '' "the path to the genome sequence directory, for motif finding; GEM parameter" ''
DEFINE_string 'outputPrefix' '' "output folder name and file name prefix" ''
DEFINE_string 'control' '' "[optional] aligned reads file for control" ''
DEFINE_string 'chrSize' '' "[optional] genome chrom.sizes file with chr name/length pairs" ''
DEFINE_string 'format' '' "[optional] read file format: BED/SAM/BOWTIE/ELAND/NOVO; default: BED" ''
DEFINE_integer 'sizeInBp' ' ' "[optional] size of mappable genome in bp (default is estimated from genome chrom sizes)" ''
DEFINE_float 'alphaValue' ' ' "[optional] minimum alpha value for sparse prior (default is esitmated from the whole dataset coverage)" ''
DEFINE_float 'qValue' '2' "[optional] significance level for q-value, specify as -log10(q-value) (default=2, q-value=0.01)" ''
DEFINE_integer 'threads' '0' "[optional] maximum number of threads to run GEM in paralell; default: #CPU" ''
DEFINE_integer 'kSeqs' '5000' "number of binding events to use for motif discovery; GEM parameter" ''
DEFINE_boolean 'useFixedAlpha' 'false' "use a fixed user-specified alpha value for all the regions" ''
DEFINE_boolean 'JASPAROutput' 'true' "output motif PFM in JASPAR format; GEM parameter" ''
DEFINE_boolean 'MEMEOutput' 'true' "output motif PFM in MEME format; GEM parameter" ''
DEFINE_boolean 'HOMEROutput' 'true' "output motif PFM in HOMER format; GEM parameter" ''
DEFINE_boolean 'BEDOutput' 'true' "output binding events in BED format for UCSC Genome Browser" ''
DEFINE_boolean 'NarrowPeakOutput' 'true' "output binding events in ENCODE NarrowPeak format" ''
DEFINE_integer 'memoryPerThread' '2048' 'total memory per thread in MB if running on local host; otherwise memory limit of executor might be set; default: 2048' ''
DEFINE_string 'workingDir' '/usr/local/storage/' 'path to working directory' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'version' 'false' "display the version of GEM" ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ] && [ ! -z "$FLAGS_jarPath" ]; then
	VER=$(java -jar "$FLAGS_jarPath" --help 2>&1 | head -n 2 | grep -oE "([0-9]+\.([0-9]\.)*[0-9]+\.?)" | head -n 1)
	echo $VER
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_jarPath" ]; then 
	echoError "Parameter --jarPath must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_readDistribution" ]; then 
	echoError "Parameter --readDistribution must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_expt" ]; then 
	echoError "Parameter --expt must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi
# not GPS only mode
if [ "$FLAGS_gpsOnly" -eq 1 ]; then
	if [ -z "$FLAGS_k" ]; then 
		echoError "Parameter --k must be set. (see --help for details)"; 
		exit $EXIT_MISSING_ARGUMENTS
	fi
	if [ -z "$FLAGS_kMin" ]; then 
		echoError "Parameter --kMin must be set. (see --help for details)"; 
		exit $EXIT_MISSING_ARGUMENTS
	fi
	if [ -z "$FLAGS_kMax" ]; then 
		echoError "Parameter --kMax must be set. (see --help for details)"; 
		exit $EXIT_MISSING_ARGUMENTS
	fi
	if [ -z "$FLAGS_seed" ]; then 
		echoError "Parameter --seed must be set. (see --help for details)"; 
		exit $EXIT_MISSING_ARGUMENTS
	fi
	if [ -z "$FLAGS_genome" ]; then 
		echoError "Parameter --genome must be set. (see --help for details)"; 
		exit $EXIT_MISSING_ARGUMENTS
	fi
fi
if [ -z "$FLAGS_outputPrefix" ]; then 
	echoError "Parameter --outputPrefix must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi

# check, if the input files exist
verifyFileExistence "$FLAGS_jarPath"
verifyFileExistence "$FLAGS_readDistribution"

if [ ! -z "$FLAGS_control" ]; then 
	verifyFileExistence "$FLAGS_control"
fi


# not GPS only mode
if [ "$FLAGS_gpsOnly" -eq 1 ]; then
	if [ ! -z "$FLAGS_chrSize" ]; then 
		verifyFileExistence "$FLAGS_chrSize"
	fi
	verifyFolderExistence "$FLAGS_genome"

	# verify that not valid integers are not set
	FLAGS_kMin=$(ensureLowerBound "$FLAGS_kMin" 1)
	FLAGS_kMax=$(ensureLowerBound "$FLAGS_kMax" 1)
	FLAGS_k=$(ensureLowerBound "$FLAGS_k" 1)
fi

# use all cores
if [ "$FLAGS_threads" -eq 0 ]; then
	FLAGS_threads=$(nproc --all)
fi

# create output folder
OUT_BASE=$(createOutputFolder "$FLAGS_outputPrefix/.dummy")
verifyFolderExistence "$OUT_BASE"

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
MEMORY=$(getMemoryForJava "$FLAGS_threads" "$FLAGS_memoryPerThread")
COMMAND="java $MEMORY -XX:+UseConcMarkSweepGC -XX:NewSize=300M -XX:MaxNewSize=300M -jar '$FLAGS_jarPath'"

# build the command
for PARAM in $__flags_longNames; do
	# ignore that parameter since it is only for the module
	if [ "$PARAM" == "jarPath" ] || [ "$PARAM" == "returnFilePath" ] || [ "$PARAM" == "memoryPerThread" ] || [ "$PARAM" == "gpsOnly" ] || [ "$PARAM" == "returnFilePath" ]; then
		continue
	fi

	# in GPS only mode
	if [ "$FLAGS_gpsOnly" -eq 0 ]; then
		if [ "$PARAM" == "k" ] || [ "$PARAM" == "kMin" ] || [ "$PARAM" == "kMax" ] || [ "$PARAM" == "seed" ] || [ "$PARAM" == "genome" ]  || [ "$PARAM" == "JASPAROutput" ] || [ "$PARAM" == "MEMEOutput" ] || [ "$PARAM" == "HOMEROutput" ]; then
			continue
		fi
	fi

	V_NAME='FLAGS_'"$PARAM"
	VALUE=${!V_NAME}
	V_NAME='__flags_'"$PARAM"'_type'
	TYPE=${!V_NAME}
	MOD=0

	# leave untouched
	if [ "$PARAM" == "expt" ] || [ "$PARAM" == "seed" ] || [ "$PARAM" == "genome" ]; then
		MOD=1
	fi

	# change the name back
	if [ "$PARAM" == "kMin" ] || [ "$PARAM" == "kMax" ] || [ "$PARAM" == "kSeqs" ]; then
		PARAM=$(echo $PARAM | sed 's/\([A-Z]\)/_\L\1/')
		MOD=1
	fi

	# special cases
	if [ "$PARAM" == "JASPAROutput" ] || [ "$PARAM" == "MEMEOutput" ] || [ "$PARAM" == "HOMEROutput" ] || [ "$PARAM" == "BEDOutput" ] || [ "$PARAM" == "NarrowPeakOutput" ]; then
		PARAM="out"$(echo $PARAM | sed 's/Output//')

		if [ "$PARAM" == "outNarrowPeak" ]; then
			PARAM="outNP"
		fi
		MOD=1
	fi

	if [ "$PARAM" == "readDistribution" ]; then
		PARAM="d"
		MOD=1
	fi
	if [ "$PARAM" == "outputPrefix" ]; then
		PARAM="out"
		MOD=1
	fi
	if [ "$PARAM" == "control" ]; then
		PARAM="ctrl"
		MOD=1
	fi
	if [ "$PARAM" == "chrSize" ]; then
		PARAM="g"
		MOD=1
	fi
	if [ "$PARAM" == "useFixedAlpha" ]; then
		PARAM="fa"
		MOD=1
	fi

	# get only first char
	if [ $MOD -eq 0 ]; then
		PARAM=$(echo $PARAM | cut -b 1)
	fi
	
	# test if value for that parameter was set
	if [ ! -z "$VALUE" ] && [ "$VALUE" != " " ]; then
		# test if flag
		if [ $TYPE -eq 1 ]; then
			# test if flag is enabled
			if  [ $VALUE -eq 0 ]; then
				COMMAND="$COMMAND --$PARAM"
			fi
		else 
			COMMAND="$COMMAND --$PARAM '$VALUE'"
		fi
	fi
done

# get tmp folder to work in
FLAGS_workingDir=${FLAGS_workingDir%/}
TMPF=$(getTmpFile GEM "$FLAGS_workingDir")
mkdir -p "$TMPF"
deleteFolderOnExit "$TMPF"

# run it
FAIL=0
MESSAGE=$(eval "cd "$TMPF" && $COMMAND" 2>&1 | tee "${FLAGS_outputPrefix}.log")
CODE=$?

# check, if tool run to end.
COUNT=$(tail -n 3 "${FLAGS_outputPrefix}.log" | grep -Ec "^Total running time:")

if [ $COUNT -ne 1 ]; then
	echoError "GEM was terminated before it finished."
	FAIL=1
else
	# check exit code
	if [ $CODE -eq 0 ]; then
		# remove SAM validation warnings
		TMP_FILE1=$(getTmpFile)
		grep -v "Ignoring SAM validation error:" "${FLAGS_outputPrefix}.log" > "$TMP_FILE1"

		# check for exception
		MESSAGE_ERROR=$($SCRIPT_FOLDER/../../core_lib/errorChecker.sh "$TMP_FILE1" 2>&1)
		CODE_ERROR=$?
		if [ $CODE_ERROR -ne 0 ]; then
			echoError "Error checker found some errors, see found errors below"
			FAIL=2
			echoError "$MESSAGE_ERROR"
		else
			# move the result files to the correct folder	
			BASE=$(basename "${FLAGS_outputPrefix}")
			mv -n "${FLAGS_outputPrefix}/$BASE/" "${FLAGS_outputPrefix}/" > /dev/null 2>&1 

			# check if GPS file is there
			if [ ! -e "${FLAGS_outputPrefix}/${BASE}.GPS_events.txt" ]; then
				echoError "GPS event file '$TMP_FOLDER/mapping.sam' was not found."
				FAIL=1
			fi

		fi
	else
		FAIL=1
	fi
fi

# check how the status is
if [ $FAIL -eq 0 ]; then
	writeParam2File "$FLAGS_returnFilePath" "GPS_events" "${FLAGS_outputPrefix}/${BASE}.GPS_events.txt"
	writeParam2File "$FLAGS_returnFilePath" "outputPrefix" "${FLAGS_outputPrefix}"
	writeParam2File "$FLAGS_returnFilePath" "outputFolder" "${FLAGS_outputPrefix}/${BASE}_outputs"
	blockUntilFileIsWritten "$FLAGS_returnFilePath"
	echoInfo "GEM run finished!"
	exit $EXIT_OK
else
	if [ $FAIL -eq 1 ]; then
		echoError "GEM run failed. See errors above and log of GEM below:"
		echoAError "error code: '$CODE'"
		# output the original message
		echoAInfo "$MESSAGE"
	else
		echoError "GEM run failed. See errors above:"
		echoAError "error code: '$CODE'"
	fi
	exit $EXIT_FAILED
fi

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
