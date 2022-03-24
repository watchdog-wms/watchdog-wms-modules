#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='STAR'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_integer 'runThreadN' '' '[optional] int: number of threads to run STAR' 'r'
DEFINE_string 'genomeDir' '' 'string: path to the directory where genome files will be generated' 'g'
DEFINE_string 'genomeFastaFiles' '' 'string(s): path(s) to the fasta files with the genome sequences, separated by spaces. These files should be plain text FASTA files, they *cannot* be zipped.' 'f'
DEFINE_string 'sjdbGTFfile' '' '[optional] string: path to the GTF file with annotations' 's'
DEFINE_integer 'sjdbOverhang' '100' '[optional] int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)' 'j'
DEFINE_string 'sjdbGTFtagExonParentTranscript' 'transcript_id' '[optional] string: GTF attribute name for parent transcript ID (default "transcript_id" works for GTF files)' 'd'
DEFINE_string 'sjdbFileChrStartEnd' '' '[optional] string(s): path to the files with genomic coordinates (chr <tab> start <tab> end <tab> strand) for the splice junction introns.' 'b'
DEFINE_integer 'genomeSAindexNbases' '' '[optional] int: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).' 'S'
DEFINE_integer 'genomeChrBinNbits' '' '[optional] int: =log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome will occupy an integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).' 'C'


DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	# TODO obtain version of important external software or remove --version flag
	MESSAGE=$(STAR --version)
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_genomeFastaFiles" ]; then
	echoError "Parameter -f must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_genomeDir" ]; then
	echoError "Parameter -g must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -n "$FLAGS_runThreadN" ] && [ "$FLAGS_runThreadN" -lt 1 ] ; then
	echoError "Parameter -r must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ -n "$FLAGS_sjdbOverhang" ] && [ "$FLAGS_sjdbOverhang" -lt 1 ] ; then
	echoError "Parameter -s must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ -n "$FLAGS_genomeSAindexNbases" ] && [ "$FLAGS_genomeSAindexNbases" -lt 1 ] ; then
	echoError "Parameter -S must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ -n "$FLAGS_genomeChrBinNbits" ] && [ "$FLAGS_genomeChrBinNbits" -lt 1 ] ; then
	echoError "Parameter -C must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi


printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# TODO: write functional part of the script here


if [ -n "$FLAGS_sjdbGTFfile" ]; then
	verifyFileExistence "$FLAGS_sjdbGTFfile";
	exit $EXIT_MISSING_ARGUMENTS
fi

# split the input array for fasta files and check if all files are there
FASTA_FILES="--genomeFastaFiles";
IFS=',' read -ra IN <<< "$FLAGS_genomeFastaFiles"
for I in "${!IN[@]}"; do 
        FILE="${IN[$I]}"
        # test if file is there and readable
        verifyFileExistence "$FILE"
        FASTA_FILES="$FASTA_FILES $FILE";
done

COMMAND="STAR --runMode genomeGenerate $FASTA_FILES --genomeDir $FLAGS_genomeDir";

if [ -n "$FLAGS_runThreadN" ]; then
	COMMAND="$COMMAND --runThreadN $FLAGS_runThreadN";
fi

if [ -n "$FLAGS_sjdbGTFfile" ]; then
	COMMAND="$COMMAND --sjdbGTFfile $FLAGS_sjdbGTFfile";
fi

if [ -n "$FLAGS_sjdbFileChrStartEnd" ]; then
	
	FILES="--sjdbFileChrStartEnd"
	IFS=',' read -ra IN <<< "$FLAGS_sjdbFileChrStartEnd"
	for I in "${!IN[@]}"; do 
			FILE="${IN[$I]}"
			# test if file is there and readable
			verifyFileExistence "$FILE"
			FILES="$FILES $FILE";
	done
	
	COMMAND="$COMMAND $FILES";
fi

if [ -n "$FLAGS_genomeChrBinNbits" ]; then
	COMMAND="$COMMAND --genomeChrBinNbits $FLAGS_genomeChrBinNbits";
fi

if [ -n "$FLAGS_genomeSAindexNbases" ]; then
	COMMAND="$COMMAND --genomeSAindexNbases $FLAGS_genomeSAindexNbases";
fi


echo $COMMAND;

MESSAGE=$($COMMAND)


RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile STARgenomeGenerate)
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
