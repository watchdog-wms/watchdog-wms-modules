#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='spring:sed:echo:tail:grep:du'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'fastq' '' 'path to one or two (PE datasets) fastq files; possible endings: .fastq, .fq, .fastq.gz or .fq.gz file' 'f'
DEFINE_string 'spring' '' 'path to compressed spring file; possible endings: .spring or .tar' 's'
DEFINE_boolean 'compress' '0' '[optional] if true the fastq files are compressed; otherwise the spring file is decompressed' 'c'
DEFINE_boolean 'preserveOrder' '0' '[optional] preserve read order; default: true' 'p'
DEFINE_boolean 'quality' '0' '[optional] retain quality values during compression; default: true' ''
DEFINE_boolean 'ids' '0' '[optional] retain read identifiers during compression; default: 0' ''
DEFINE_string 'qualityMode' 'lossless' "[optional] possible values: 'lossless', 'qvz qv_ratio', 'ill_bin' or 'binary thr high low'; default: lossless" 'q'
DEFINE_boolean 'long' '1' '[optional] use for compression of arbitrarily long reads; default: false' 'l'
DEFINE_string 'decompressRange' '' "[optional] decompress only reads (or read pairs for PE datasets) from start to end (both inclusive); e.g. '1 100'" 'r'
DEFINE_integer 'threads' '1' '[optional] number of cores to use' 't'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_string 'workingDir' '/usr/local/storage/' 'path to working directory' 'w'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# mode parameters (1 == set)
GZIP_MODE=0
DECOMPRESS_MODE=$FLAGS_compress

# check if mandatory arguments are there
if [ -z "$FLAGS_fastq" ]; then
	echoError "Parameter -f must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_spring" ]; then
	echoError "Parameter -s must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

# check of optional parameters
FLAGS_spring=$(abspath "${FLAGS_spring}")
FLAGS_fastq=$(echo "$FLAGS_fastq" | sed 's/,\W*$//')
IFS=',' read -ra FASTQ <<< "$FLAGS_fastq"
IFS=$'\n'

# get absolute path
for I in "${!FASTQ[@]}"; do 
	FASTQ[$I]=$(abspath "${FASTQ[$I]}")
done

# test if one or two fastq files are given
if [ "${#FASTQ[@]}" -eq 2 ]; then
	FASTQ_FILES="'${FASTQ[0]}' '${FASTQ[1]}'"
	PAIRED="true"
else
	if [ "${#FASTQ[@]}" -ne 1 ]; then
		echoError "Parameter -f does only accept one or two fastq files (see --help for details)";
		exit $EXIT_MISFORMATED_INPUT
	fi
	FASTQ_FILES="'${FASTQ[0]}'"
	PAIRED="false"
fi

# test fastq file endings
if [ $DECOMPRESS_MODE -eq 0 ]; then
	echoInfo "file size of input files in MB:"
fi
for I in "${!FASTQ[@]}"; do 
	F="${FASTQ[$I]}"
	if [[ "$F" == *.fastq ]] || {[ "$F" == *.fq ]] || [{ "$F" == *.fastq.gz ]] || [[ "$F" == *.fq.gz ]]; then
		if [[ "$F" == *.fastq.gz ]] || [[ "$F" == *.fq.gz ]]; then
			GZIP_MODE=1
		fi
	else
		echoError "Parameter -f must be have .fastq, .fq, .fastq.gz or .fq.gz ending. (see --help for details)";
		exit $EXIT_INVALID_ARGUMENTS
	fi

	# ensure that the file is there
	if [ $DECOMPRESS_MODE -eq 0 ]; then
		verifyFileExistence "$F"
		# determine size of input fastq files in MB
		du -B 1M "$F"
	else
		if [ -e "$F" ]; then
			echoError "Output file '$F' already exists!";
			exit $EXIT_INVALID_ARGUMENTS
		fi
	fi
done

# ensure that spring file ending is provided
if [[ "$FLAGS_spring" != *.spring ]] && [[ "$FLAGS_spring" != *.tar ]]; then
	echoError "Parameter -s must be have .spring or .tar ending. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
else
	# ensure that the file is there
	if [ $DECOMPRESS_MODE -eq 1 ]; then
		verifyFileExistence "$FLAGS_spring"
	else
		if [ -e "$FLAGS_spring" ]; then
			echoError "Output file '$FLAGS_spring' already exists!";
			exit $EXIT_INVALID_ARGUMENTS
		fi
	fi
	
fi
if [ "$FLAGS_threads" -gt 32 ] || [ "$FLAGS_threads" -lt 1 ]; then
	echoError "Parameter -t must be between [1, 32]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ ! -z "$FLAGS_decompressRange" ]; then
	IFS=' ' read -ra RANGE <<< "$FLAGS_decompressRange"
	IFS=$'\n'
	if [ "${#RANGE[@]}" -ne 2 ] || [ "${RANGE[1]}" -lt "${RANGE[0]}" ]; then
		echoError "Parameter -r must contain two integers with fistN <= lastN. (see --help for details)";
		exit $EXIT_INVALID_ARGUMENTS
	fi
fi
if [ "$FLAGS_qualityMode" != "lossless" ]; then
	echoWarn "Parameter -q is not validated. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# get base output name
FLAGS_workingDir=${FLAGS_workingDir%/}
B=$(basename "${FLAGS_spring}")
BASENAME=${B%.*}

# get temporary folder
TMP=$(getTmpFile spring "$FLAGS_workingDir")
mkdir -p "$TMP"
deleteFolderOnExit "$TMP"
cd "$TMP"

# build the command
COMMAND="spring"

if [ "$DECOMPRESS_MODE" -eq 1 ]; then
	COMMAND="$COMMAND --decompress --input-file '$FLAGS_spring' --output-file $FASTQ_FILES"
	PREFIX_MODE="de"
else
	COMMAND="$COMMAND --compress --input-file $FASTQ_FILES --output-file '$FLAGS_spring'"
	PREFIX_MODE=""
fi
if [ "$FLAGS_preserveOrder" -eq 1 ]; then
	COMMAND="$COMMAND --allow-read-reordering"
fi
if [ "$FLAGS_quality" -eq 1 ]; then
	COMMAND="$COMMAND --no-quality"
fi
if [ "$FLAGS_ids" -eq 1 ]; then
	COMMAND="$COMMAND --no-ids"
fi
if [ ! -z "$FLAGS_qualityMode" ]; then
	COMMAND="$COMMAND --quality-opts '$FLAGS_qualityMode'"
fi
if [ "$FLAGS_long" -eq 1 ]; then
	COMMAND="$COMMAND --long"
fi
if [ "$GZIP_MODE" -eq 1 ]; then
	COMMAND="$COMMAND --gzipped_fastq"
fi
if [ ! -z "$FLAGS_decompressRange" ]; then
	COMMAND="$COMMAND --decompress-range '$FLAGS_decompressRange'"
fi

COMMAND="$COMMAND --num-threads '$FLAGS_threads' --working-dir '$FLAGS_workingDir'"
echo $COMMAND

# execute the command
LOG="$(abspath "$(dirname "${FLAGS_spring}")")/${BASENAME}.${PREFIX_MODE}compress.spring.log"
if [ -e "$LOG" ]; then
	rm "$LOG"
fi
executeCommand "$COMMAND" "$LOG" "SPRING"


# some additional tests
CHECKS=0
COUNT=$(tail -n 5 "$LOG" | grep -c -e "^Total time for ")
if [ $COUNT -ne 1 ]; then
	echoError "Spring was terminated before it finished."
	CHECKS=1
fi

# test if files are there
if [ $DECOMPRESS_MODE -eq 1 ]; then
	for I in "${!FASTQ[@]}"; do 
		F="${FASTQ[$I]}"
		if [ ! -e "$F" ]; then 
			echoError "Output file '$F' is missing!"
			CHECKS=1
		fi
	done
	OUTFILE="${FASTQ[0]},${FASTQ[1]}"
else
	if [ ! -e "$FLAGS_spring" ]; then 
		echoError "Output file '$FLAGS_spring' is missing!"
		CHECKS=1
	fi
	OUTFILE="$FLAGS_spring"
fi

# check status
if [ $CHECKS -eq 0 ]; then
	if [ ! -z "$FLAGS_returnFilePath" ]; then
		writeParam2File "$FLAGS_returnFilePath" "createdFile" "$OUTFILE"
		writeParam2File "$FLAGS_returnFilePath" "isPairedEnd" "$PAIRED"
		blockUntilFileIsWritten "$FLAGS_returnFilePath"
	fi
	exit $EXIT_OK
else
	echoError "SPRING run failed!"
	if [ $DECOMPRESS_MODE -eq 1 ]; then
		for I in "${!FASTQ[@]}"; do 
			F="${FASTQ[$I]}"
			rm -f "$F"
		done
	else
		rm -f "$FLAGS_spring"
	fi
	exit $EXIT_FAILED
fi

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
