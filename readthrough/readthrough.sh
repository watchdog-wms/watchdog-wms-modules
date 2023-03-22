#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='java'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'annotation' '' 'annotation file path' 'a'
DEFINE_string 'genecounts' '' 'gene read count file' 'g'
DEFINE_string 'input' '' 'input file' 'i'
DEFINE_integer 'readinLength' '5000' '[optional] length of upstream window in which read-in is calculated' 'n'
DEFINE_string 'output' '' 'output file' 'o'
DEFINE_integer 'strandedness' '0' 'strandedness: 0=not strandspecific, 1=first read indicates strand, 2=second read indicates strand' 's'
DEFINE_integer 'readthroughLength' '5000' '[optional] length of downstream window in which read-through is calculated' 't'
DEFINE_integer 'overlap' '25' '[optional] minimum overlap of read to be counted for read-through/in window' 'v'
DEFINE_string 'idxstats' '' '[optional] idxstats file with numbers of mapped reads per chromosome, necessary for calculating downstream FPKM and transcription in dOCR regions' 'x'
DEFINE_string 'normFactor' '' '[optional] factor for normalizing by mapped reads and gene length for downstream FPKM calculation' 'm'
DEFINE_string 'exclude' '' '[optional] chromosomes to exclude from calculating total mapped reads, separated by ,' 'e'
DEFINE_string 'excludeType' '' '[optional] gene types to exclude when determining genes with no other genes up- or down-stream, separated by ,' 'z'
DEFINE_string 'dOCRFile' '' '[optional] file containing dOCR lengths' 'd'
DEFINE_integer 'windowLength' '1000' '[optional] number of steps for evaluating transcription on dOCRs' 'w'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_annotation" ]; then
	echoError "Parameter -a must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_genecounts" ]; then
	echoError "Parameter -g must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_input" ]; then
	echoError "Parameter -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_readinLength" -gt 1000000000 ] || [ "$FLAGS_readinLength" -lt 1 ]; then
	echoError "Parameter -n must be between [1, 1000000000]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ -z "$FLAGS_output" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_strandedness" -gt 2 ] || [ "$FLAGS_strandedness" -lt 0 ]; then
	echoError "Parameter -s must be between [0, 2]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_readthroughLength" -gt 1000000000 ] || [ "$FLAGS_readthroughLength" -lt 1 ]; then
	echoError "Parameter -t must be between [1, 1000000000]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_overlap" -gt 200 ] || [ "$FLAGS_overlap" -lt 1 ]; then
	echoError "Parameter -v must be between [1, 200]. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi



printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# check, if the input files exist
verifyFileExistence "$FLAGS_annotation"
verifyFileExistence "$FLAGS_genecounts"
verifyFileExistence "$FLAGS_input"


cd $SCRIPT_FOLDER/dist;

MEMORY=$(getMemoryForJava 1 4096)
COMMAND="java $MEMORY -cp dOCRTools.jar de.lmu.ifi.bio.docrtools.ReadthroughCalculatorForWatchdog -a '$FLAGS_annotation' -o '$FLAGS_output' -i '$FLAGS_input' -g '$FLAGS_genecounts' -s '$FLAGS_strandedness' -v '$FLAGS_overlap' -n '$FLAGS_readinLength' -t '$FLAGS_readthroughLength'";

echo $COMMAND;

if [ -n "$FLAGS_idxstats" ] ; then
    COMMAND="$COMMAND -x '$FLAGS_idxstats'"
fi

if [ -n "$FLAGS_normFactor" ] ; then
    COMMAND="$COMMAND -m '$FLAGS_normFactor'"
fi

if [ -n "$FLAGS_exclude" ] ; then
    COMMAND="$COMMAND -e '$FLAGS_exclude'"
fi

if [ -n "$FLAGS_excludeType" ] ; then
    COMMAND="$COMMAND -z '$FLAGS_excludeType'"
fi

if [ -n "$FLAGS_dOCRFile" ] ; then
    COMMAND="$COMMAND -d '$FLAGS_dOCRFile' -w '$FLAGS_windowLength'"
fi


echo $COMMAND;

MESSAGE=$(eval $COMMAND 2>&1)
RET=$?


# check exit code
if [ $RET -eq 0 ]; then
        # output the original message
        printf "$MESSAGE\n"

        exit $EXIT_OK
else
        echoError "readthrough calculation run failed. Output file was deleted."
        echoAError "error code: '$RET'"
         # output the original message
        printf "$MESSAGE\n" 1>&2
        rm -f "$FLAGS_output*" 2>&1 > /dev/null
        exit $EXIT_FAILED
fi


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
