#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='echo:printf:rm:grep:head:wc:picard:java'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'returnFilePath' '' "path to the return variables file" ''
DEFINE_string 'Histogram_FILE' '' "File to write insert size Histogram chart to.  Required." 'H'
DEFINE_string 'INPUT' '' "Input SAM/BAM/CRAM file.  Required." 'I'
DEFINE_string 'OUTPUT' '' "The file to write the output to.  Required." 'O'
DEFINE_string 'arguments_file' '' "[optional] read one or more arguments files and add them to the command line  This argument may be specified 0 or more times. Default value: null." ''
DEFINE_boolean 'ASSUME_SORTED' '0' "[optional] If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}" 'AS'
DEFINE_integer 'COMPRESSION_LEVEL' '' "[optional] Compression level for all compressed files created (e.g. BAM and VCF).  Default value: 5." ''
DEFINE_boolean 'CREATE_INDEX' '1' "[optional] Whether to create an index when writing VCF or coordinate sorted BAM output.  Default value: false. Possible values: {true, false}" ''
DEFINE_boolean 'CREATE_MD5_FILE' '1' "[optional] Whether to create an MD5 digest for any BAM or FASTQ files created. Default value:false. Possible values: {true, false}" ''
DEFINE_float 'DEVIATIONS' '' "[optional] Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0." ''
DEFINE_string 'GA4GH_CLIENT_SECRETS' '' "[optional] Google Genomics API client_secrets.json file path.  Default value: client_secrets.json." ''
DEFINE_integer 'HISTOGRAM_WIDTH' '' "[optional] Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null." 'W'
DEFINE_boolean 'INCLUDE_DUPLICATES' '1' "[optional] If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false}" ''
DEFINE_integer 'MAX_RECORDS_IN_RAM' '' "[optional] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000." ''
DEFINE_string 'METRIC_ACCUMULATION_LEVEL' '' "[optional] The level(s) at which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP}" 'LEVEL'
DEFINE_integer 'MIN_HISTOGRAM_WIDTH' '' "[optional] Minimum width of histogram plots. In the case when the histogram would otherwise betruncated to a shorter range of sizes, the MIN_HISTOGRAM_WIDTH will enforce a minimum range.  Default value: null. " 'MW'
DEFINE_float 'MINIMUM_PCT' '' "[optional] When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05." 'M'
DEFINE_boolean 'QUIET' '1' "[optional] Whether to suppress job-summary info on System.err.  Default value: false. Possible values: {true, false}" ''
DEFINE_string 'REFERENCE_SEQUENCE' '' "[optional] Reference sequence file.  Default value: null." 'R'
DEFINE_integer 'STOP_AFTER' '' "[optional] Stop after processing N reads, mainly for debugging.  Default value: 0." ''
DEFINE_string 'TMP_DIR' '' "[optional] One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null." ''
DEFINE_boolean 'USE_JDK_DEFLATER' '1' "[optional] Use the JDK Deflater instead of the Intel Deflater for writing compressed output. Default value: false. Possible values: {true, false}" 'use_jdk_deflater'
DEFINE_boolean 'USE_JDK_INFLATER' '1' "[optional] Use the JDK Inflater instead of the Intel Inflater for reading compressed input. Default value: false. Possible values: {true, false}" 'use_jdk_inflater'
DEFINE_string 'VALIDATION_STRINGENCY' 'STRICT' "[optional] Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT} " ''
DEFINE_string 'VERBOSITY' '' "[optional] Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}" ''
DEFINE_boolean 'version' '1' "[optional] display the version number for this tool" 'v'
DEFINE_boolean 'debug' '1' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

tool="picard CollectInsertSizeMetrics";


if [ "$FLAGS_version" -eq 0 ]; then
	# TODO obtain version of important external software or remove --version flag
	MESSAGE=$($tool --version true)
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_Histogram_FILE" ]; then 
	echoError "Parameter Histogram_FILE must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_INPUT" ]; then 
	echoError "Parameter INPUT must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_OUTPUT" ]; then 
	echoError "Parameter OUTPUT must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi


# further parameter checks  which refer to input ranges or valid parameters in general


# check, if the input files exist
verifyFileExistence "$FLAGS_INPUT"



printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# functional part of the script
flagsAsString=""
if [ "$FLAGS_ASSUME_SORTED" == 0 ]; then
 flagsAsString="$flagsAsString --AS true"
fi
if [ "$FLAGS_CREATE_INDEX" == 0 ]; then
 flagsAsString="$flagsAsString --CREATE_INDEX true"
fi
if [ "$FLAGS_CREATE_MD5_FILE" == 0 ]; then
 flagsAsString="$flagsAsString --CREATE_MD5_FILE true"
fi
if [ "$FLAGS_INCLUDE_DUPLICATES" == 0 ]; then
 flagsAsString="$flagsAsString --INCLUDE_DUPLICATES true"
fi
if [ "$FLAGS_QUIET" == 0 ]; then
 flagsAsString="$flagsAsString --QUIET true"
fi
if [ "$FLAGS_USE_JDK_DEFLATER" == 0 ]; then
 flagsAsString="$flagsAsString --use_jdk_deflater true"
fi
if [ "$FLAGS_USE_JDK_INFLATER" == 0 ]; then
 flagsAsString="$flagsAsString --use_jdk_inflater true"
fi
if  [ ! -z "$FLAGS_Histogram_FILE" ]; then
	flagsAsString="$flagsAsString --H $FLAGS_Histogram_FILE"
fi
if  [ ! -z "$FLAGS_INPUT" ]; then
	flagsAsString="$flagsAsString --I $FLAGS_INPUT"
fi
if  [ ! -z "$FLAGS_OUTPUT" ]; then
	flagsAsString="$flagsAsString --O $FLAGS_OUTPUT"
fi
if  [ ! -z "$FLAGS_arguments_file" ]; then
	flagsAsString="$flagsAsString --arguments_file $FLAGS_arguments_file"
fi
if  [ ! -z "$FLAGS_COMPRESSION_LEVEL" ]; then
	flagsAsString="$flagsAsString --COMPRESSION_LEVEL $FLAGS_COMPRESSION_LEVEL"
fi
if  [ ! -z "$FLAGS_DEVIATIONS" ]; then
	flagsAsString="$flagsAsString --DEVIATIONS $FLAGS_DEVIATIONS"
fi
if  [ ! -z "$FLAGS_GA4GH_CLIENT_SECRETS" ]; then
	flagsAsString="$flagsAsString --GA4GH_CLIENT_SECRETS $FLAGS_GA4GH_CLIENT_SECRETS"
fi
if  [ ! -z "$FLAGS_HISTOGRAM_WIDTHW" ]; then
	flagsAsString="$flagsAsString --W $FLAGS_HISTOGRAM_WIDTHW"
fi
if  [ ! -z "$FLAGS_MAX_RECORDS_IN_RAM" ]; then
	flagsAsString="$flagsAsString --MAX_RECORDS_IN_RAM $FLAGS_MAX_RECORDS_IN_RAM"
fi
if  [ ! -z "$FLAGS_METRIC_ACCUMULATION_LEVEL" ]; then
	flagsAsString="$flagsAsString --LEVEL $FLAGS_METRIC_ACCUMULATION_LEVEL"
fi
if  [ ! -z "$FLAGS_MIN_HISTOGRAM_WIDTH" ]; then
	flagsAsString="$flagsAsString --MW $FLAGS_MIN_HISTOGRAM_WIDTH"
fi
if  [ ! -z "$FLAGS_MINIMUM_PCT" ]; then
	flagsAsString="$flagsAsString --M $FLAGS_MINIMUM_PCT"
fi
if  [ ! -z "$FLAGS_REFERENCE_SEQUENCE" ]; then
	flagsAsString="$flagsAsString --R $FLAGS_REFERENCE_SEQUENCE"
fi
if  [ ! -z "$FLAGS_STOP_AFTER" ]; then
	flagsAsString="$flagsAsString --STOP_AFTER $FLAGS_STOP_AFTER"
fi
if  [ ! -z "$FLAGS_TMP_DIR" ]; then
	flagsAsString="$flagsAsString --TMP_DIR $FLAGS_TMP_DIR"
fi
if  [ ! -z "$FLAGS_VALIDATION_STRINGENCY" ]; then
	flagsAsString="$flagsAsString --VALIDATION_STRINGENCY $FLAGS_VALIDATION_STRINGENCY"
fi
if  [ ! -z "$FLAGS_VERBOSITY" ]; then
	flagsAsString="$flagsAsString --VERBOSITY $FLAGS_VERBOSITY"
fi
# run it
MESSAGE=$($tool $flagsAsString)


RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile insertSizeMetrics)
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
		writeParam2File "$FLAGS_returnFilePath" "outputBamFile" "TODO: add value"
		writeParam2File "$FLAGS_returnFilePath" "outputHistogramFile" "TODO: add value"
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
