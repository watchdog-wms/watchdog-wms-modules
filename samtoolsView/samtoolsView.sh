#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='echo:printf:rm:grep:head:wc:samtools'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'bamoutput' '1' "output BAM" 'b'
DEFINE_boolean 'cramoutput' '1' "output CRAM (requires reference sequence)" 'C'
DEFINE_boolean 'fastCompression' '1' "use fast BAM compression (implies bamoutput)" '1'
DEFINE_boolean 'uncompressedBam' '1' "uncompressed BAM output (implies bamoutput)" 'u'
DEFINE_boolean 'includeHeader' '1' "include header in SAM output" 'h1'
DEFINE_boolean 'printOnlyHeader' '1' "print SAM header only (no alignments)" 'H'
DEFINE_boolean 'printCounts' '1' "print only the count of matching records" 'c'
DEFINE_string 'output' '' "output file name" 'o'
DEFINE_string 'outputReadsNotSelected' '' "output reads not selected by filters to FILE" 'U'
DEFINE_string 'referenceLengths' '' "FILE listing reference names and lengths (see long help)" 't'
DEFINE_string 'bedfile' '' "only include reads overlapping this BED FILE" 'L'
DEFINE_string 'readgroup' '' "only include reads in read group STR" 'r'
DEFINE_string 'readgroupFile' '' "only include reads with read group listed in FILE" 'R'
DEFINE_integer 'mappingquality' '0' "only include reads with mapping quality at least INT" 'q'
DEFINE_string 'library' '' "only include reads in library STR" 'l'
DEFINE_integer 'minquerylength' '' "only include reads with number of CIGAR operations consuming  query sequence at least INT" 'm'
DEFINE_integer 'bitsset' '0' "only include reads with all bits set in INT set in FLAG" 'f'
DEFINE_integer 'bitsnotset' '0' "only include reads with none of the bits set in INT set in FLAG" 'F'
DEFINE_string 'readTagToStrip' '' "read tag to strip (repeatable)" 'x'
DEFINE_string 'collapseCIGAROperation' '' "collapse the backward CIGAR operation" 'B'
DEFINE_float 'seed' '0' "integer part sets seed of random number generator, rest sets fraction of templates to subsample" 's'
DEFINE_string 'threads' '' "number of BAM/CRAM compression threads" '@'
DEFINE_string 'printLongHelp' '' "print long help, including note about region specification" 'h2'
DEFINE_string 'inputfmtoption' '' "Specify a single input file format option in the form  of OPTION or OPTION=VALUE" ''
DEFINE_string 'outputfmt' '' "Specify output format (SAM, BAM, CRAM)" 'O'
DEFINE_string 'outputfmtoption' '' "Specify a single output file format option in the form  of OPTION or OPTION=VALUE" ''
DEFINE_string 'reference' '' "Reference sequence FASTA FILE" 'T'
DEFINE_string 'inbam' '' "input bam file" ''
DEFINE_string 'insam' '' "input sam file" ''
DEFINE_string 'incram' '' "input cram file" ''
DEFINE_string 'region' '' "region selected" ''
DEFINE_boolean 'version' 'false' 'prints the version' 'v'
DEFINE_boolean 'debug' 'false' 'prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
        MESSAGE=$(samtools 2>&1 | head -n 3 | tail -n 1)
        echo $MESSAGE
        exit $EXIT_OK
fi


# check if mandatory arguments are there
if [ -z "$FLAGS_inbam" ] && [ -z "$FLAGS_insam" ] && [ -z "$FLAGS_incram" ]; then 
	echoError "One of inbam, insam or incram must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_output" ] ; then 
	echoError "Output file must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi


# further parameter checks  which refer to input ranges or valid parameters in general


# check, if the input files exist



printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# functional part of the script
flagsAsString=""
if  [ "$FLAGS_bamoutput" -eq 0 ]; then
	flagsAsString="$flagsAsString -b"
fi
if  [ "$FLAGS_cramoutput" -eq 0 ]; then
	flagsAsString="$flagsAsString -C"
fi
if  [ "$FLAGS_fastCompression" -eq 0 ]; then
	flagsAsString="$flagsAsString -1"
fi
if  [ "$FLAGS_uncompressedBam" -eq 0 ]; then
	flagsAsString="$flagsAsString -u"
fi
if  [ "$FLAGS_includeHeader" -eq 0 ]; then
	flagsAsString="$flagsAsString -h"
fi
if  [ "$FLAGS_printOnlyHeader" -eq 0 ]; then
	flagsAsString="$flagsAsString -H"
fi
if  [ "$FLAGS_printCounts" -eq 0 ]; then
	flagsAsString="$flagsAsString -c"
fi
if  [ ! -z "$FLAGS_output" ]; then
	flagsAsString="$flagsAsString -o  $FLAGS_output"
fi
if  [ ! -z "$FLAGS_outputReadsNotSelected" ]; then
	flagsAsString="$flagsAsString -U  $FLAGS_outputReadsNotSelected"
fi
if  [ ! -z "$FLAGS_referenceLengths" ]; then
	flagsAsString="$flagsAsString -t  $FLAGS_referenceLengths"
fi
if  [ ! -z "$FLAGS_bedfile" ]; then
	flagsAsString="$flagsAsString -L  $FLAGS_bedfile"
fi
if  [ ! -z "$FLAGS_readgroup" ]; then
	flagsAsString="$flagsAsString -r  $FLAGS_readgroup"
fi
if  [ ! -z "$FLAGS_readgroupFile" ]; then
	flagsAsString="$flagsAsString -R  $FLAGS_readgroupFile"
fi
if  [ ! -z "$FLAGS_mappingquality" ]; then
	flagsAsString="$flagsAsString -q  $FLAGS_mappingquality"
fi
if  [ ! -z "$FLAGS_library" ]; then
	flagsAsString="$flagsAsString -l  $FLAGS_library"
fi
if  [ ! -z "$FLAGS_minquerylength" ]; then
	flagsAsString="$flagsAsString -m  $FLAGS_minquerylength"
fi
if  [ ! -z "$FLAGS_bitsset" ]; then
	flagsAsString="$flagsAsString -f  $FLAGS_bitsset"
fi
if  [ ! -z "$FLAGS_bitsnotset" ]; then
	flagsAsString="$flagsAsString -F  $FLAGS_bitsnotset"
fi
if  [ ! -z "$FLAGS_readTagToStrip" ]; then
	flagsAsString="$flagsAsString -x  $FLAGS_readTagToStrip"
fi
if  [ ! -z "$FLAGS_collapseCIGAROperation" ]; then
	flagsAsString="$flagsAsString -B  $FLAGS_collapseCIGAROperation"
fi
if  [ ! -z "$FLAGS_seed" ]; then
	flagsAsString="$flagsAsString -s  $FLAGS_seed"
fi
if  [ ! -z "$FLAGS_threads" ]; then
	flagsAsString="$flagsAsString -@  $FLAGS_threads"
fi
if  [ ! -z "$FLAGS_printLongHelp" ]; then
	flagsAsString="$flagsAsString -?  $FLAGS_printLongHelp"
fi
if  [ ! -z "$FLAGS_inputfmtoption" ]; then
	flagsAsString="$flagsAsString --input-fmt-option  $FLAGS_inputfmtoption"
fi
if  [ ! -z "$FLAGS_outputfmt" ]; then
	flagsAsString="$flagsAsString -O  $FLAGS_outputfmt"
fi
if  [ ! -z "$FLAGS_outputfmtoption" ]; then
	flagsAsString="$flagsAsString --output-fmt-option  $FLAGS_outputfmtoption"
fi
if  [ ! -z "$FLAGS_reference" ]; then
	flagsAsString="$flagsAsString -T  $FLAGS_reference"
fi
if  [ ! -z "$FLAGS_inbam" ]; then
	flagsAsString="$flagsAsString  $FLAGS_inbam"
fi
if  [ ! -z "$FLAGS_insam" ]; then
	flagsAsString="$flagsAsString  $FLAGS_insam"
fi
if  [ ! -z "$FLAGS_incram" ]; then
	flagsAsString="$flagsAsString  $FLAGS_incram"
fi
if  [ ! -z "$FLAGS_region" ]; then
	flagsAsString="$flagsAsString  $FLAGS_region"
fi

# run it
MESSAGE=$(samtools view $flagsAsString)


RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile samtoolsView)
touch "$TMP_FILE"
trap "rm -f \"$TMP_FILE\" 2>&1 > /dev/null" EXIT
printf "$MESSAGE" > "$TMP_FILE"
MESSAGE_ERROR=$($SCRIPT_FOLDER/../../core_lib/errorChecker.sh "$TMP_FILE" "truncated file:fail to open file" 2>&1)
CODE_ERROR=$?
rm -f "$TMP_FILE" 2>&1 > /dev/null

if [ $CODE_ERROR -ne 0 ]; then
	echoError "Error checker found some errors, see found errors below"
	echo -e "$MESSAGE_ERROR"
	exit $EXIT_FAILED
else
	# check exit code
	if [ $FAIL -eq 0 ] && [ $RET -eq 0 ]; then
		# output the original message
		printf "$MESSAGE\n"
		writeParam2File "$FLAGS_returnFilePath" "outputFile" "$FLAGS_output"
		blockUntilFileIsWritten "$FLAGS_returnFilePath"
		exit $EXIT_OK
	else
		FAIL=1
	fi
	if [ $FAIL -eq 1 ]; then
		echoError "Run failed. See error below"
		echoAError "error code: '$RET'"
		 # output the original message
		printf "$MESSAGE\n"
		exit $EXIT_FAILED
	fi
fi


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
