#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='samtools:java'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'infile' '' 'input file' 'i'
DEFINE_string 'method' 'mpileup2snp' 'method: mpileup2snp, mpileup2indel or mpileup2cns' 'd'
DEFINE_string 'jar' '' 'jar file' 'j'
DEFINE_string 'reference' '' 'reference sequence' 'g'
DEFINE_integer 'minCoverage' '' 'paramIntegerRange_varScanMpileup' 'm'
DEFINE_integer 'minReads2' '' 'Minimum supporting reads at a position to call variants' 'r'
DEFINE_integer 'minAvgQual' '' 'Minimum base quality at a position to count a read ' 'a'
DEFINE_float 'minVarFreq' '' 'Minimum variant allele frequency threshold' 'v'
DEFINE_float 'minFreqForHom' '' 'Minimum frequency to call homozygote' 'f'
DEFINE_float 'pValue' '' 'Default p-value threshold for calling variants' 'p'
DEFINE_integer 'strandFilter' '' 'Ignore variants with >90% support on one strand' 's'
DEFINE_integer 'outputVcf' '' 'If set to 1, outputs in VCF format' 'o'
DEFINE_string 'vcfSampleList' '' 'For VCF output, a list of sample names in order, one per line' 'l'
DEFINE_integer 'variants' '' 'Report only variant (SNP/indel) positions' 'V'
DEFINE_boolean 'version' '1' "print version information and quit" ''
DEFINE_string 'output' '' 'output file' 'O'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
        
        MESSAGE=$(java -jar $FLAGS_jar 2>&1 | grep '^VarScan' | cut -d " " -f 2);

        echo $MESSAGE
        exit $EXIT_OK
fi


# check if mandatory arguments are there
if [ -z "$FLAGS_infile" ]; then
	echoError "Parameter --infile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_jar" ]; then
	echoError "Parameter --jar must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

if [ -z "$FLAGS_reference" ]; then
	echoError "Parameter --reference must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter --output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi


if [ -z "$FLAGS_output" ]; then
	echoError "Parameter --output must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi


printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# create output directory if not already exists
dirname_out=$(dirname "$FLAGS_output")
if [ ! -d "$dirname_out" ]; then
    mkdir -p "$dirname_out"
fi

# check, if the input files exist
verifyFileExistence "$FLAGS_reference"

SAMTOOLS_CMD="samtools mpileup -f $FLAGS_reference";

IFS=',' read -ra FILES <<< "$FLAGS_infile"
for I in "${!FILES[@]}"; do 
        FILE="${FILES[$I]}"
        verifyFileExistence "$FILE"
        SAMTOOLS_CMD="$SAMTOOLS_CMD $FILE";
done



flagsAsString=""
if  [ ! -z "$FLAGS_minCoverage" ]; then
        flagsAsString="$flagsAsString --min-coverage $FLAGS_minCoverage"
fi
if  [ ! -z "$FLAGS_minReads2" ]; then
        flagsAsString="$flagsAsString --min-reads2 $FLAGS_minReads2"
fi
if  [ ! -z "$FLAGS_minAvgQual" ]; then
        flagsAsString="$flagsAsString --min-avg-qual $FLAGS_minAvgQual"
fi
if  [ ! -z "$FLAGS_minVarFreq" ]; then
        flagsAsString="$flagsAsString --min-var-freq $FLAGS_minVarFreq"
fi
if  [ ! -z "$FLAGS_minFreqForHom" ]; then
        flagsAsString="$flagsAsString --min-freq-for-hom $FLAGS_minFreqForHom"
fi
if  [ ! -z "$FLAGS_pValue" ]; then
        flagsAsString="$flagsAsString --p-value $FLAGS_pValue"
fi
if  [ ! -z "$FLAGS_strandFilter" ]; then
        flagsAsString="$flagsAsString --strand-filter $FLAGS_strandFilter"
fi
if  [ ! -z "$FLAGS_outputVcf" ]; then
        flagsAsString="$flagsAsString --output-vcf $FLAGS_outputVcf"
fi
if  [ ! -z "$FLAGS_vcfSampleList" ]; then
        flagsAsString="$flagsAsString --vcf-sample-list $FLAGS_vcfSampleList"
fi
if  [ ! -z "$FLAGS_variants" ]; then
        flagsAsString="$flagsAsString --variants $FLAGS_variants"
fi


# run it
MESSAGE=$($SAMTOOLS_CMD | java -jar $FLAGS_jar $FLAGS_method $flagsAsString > $FLAGS_output);


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
