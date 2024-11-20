#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='samtools:awk'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'bam' '' 'Path to bam file, where clipped reads should be extracted' 'b'
DEFINE_string 'out' '' 'Path to output bam file, were only clipped reads are stored.' 'o'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE="clippedReadsExtractor 1.0"
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_out" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_bam" ]; then
	echoError "Parameter -b must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

# test if file is there and readable
verifyFileExistence "$FLAGS_bam"

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

dirname=$(dirname "$FLAGS_out")
if [ ! -d "$dirname" ]; then
    mkdir -p "$dirname"
fi

# Generate unique file names
header="${dirname}/header_tmp_$(basename "$FLAGS_bam" .bam).sam"
clipped="${dirname}/clipped_tmp_$(basename "$FLAGS_bam" .bam).sam"


# Get BAM header
samtools view -H "$FLAGS_bam" > "$header"
# Extract clipped reads from BAM file
samtools view "$FLAGS_bam" | awk '$6 ~ /S|H/' > "$clipped"
# Combine header and clipped reads to new BAM file
cat "$header" "$clipped" | samtools view -Sb - > "$FLAGS_out"
# Index resulting BAM file
samtools index "$FLAGS_out"
# Remove header and clipped reads files (they are just temporary)
rm "$header"
rm "$clipped"

# If no error occurred previously it should be fine
exit $EXIT_OK

# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
