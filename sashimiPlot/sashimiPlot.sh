#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='echo:printf:rm:grep:head:wc:samtools:python:R:sashimi-plot.py'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

USED_TOOLS='^R|ggplot2|data.table|gridExtra'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS")
CODE=$?

if [ $CODE -ne 0 ]; then
        echoError "$MESSAGE"
        exit $EXIT_TOOLS_MISSING
fi


# define parameters
DEFINE_string 'help' '' "show this help message and exit" 'h1'
DEFINE_string 'bam' '' "Individual bam file or file with a list of bam files.  In the case of a list of files the format is tsv:  1col: id for bam file, 2col: path of bam file, 3+col:  additional columns" 'b'
DEFINE_string 'coordinates' '' "Genomic region. Format: chr:start-end (1-based)" 'c'
DEFINE_string 'outprefix' 'sashimi' "Prefix for plot file name" 'o'
DEFINE_string 'outstrand' 'both' "Only for --strand other than 'NONE'. Choose which  signal strand to plot: <both> <plus> <minus>" 'S'
DEFINE_integer 'mincoverage' '1' "Minimum number of reads supporting a junction to be  drawn" 'M'
DEFINE_string 'junctionsbed' '' "Junction BED file name" 'j'
DEFINE_string 'gtf' '' "Gtf file with annotation (only exons is enough)" 'g'
DEFINE_string 'strand' 'NONE' "Strand specificity: <NONE> <SENSE> <ANTISENSE>  <MATE1_SENSE> <MATE2_SENSE>" 's'
DEFINE_boolean 'shrink' '1' "Shrink the junctions by a factor for nicer display" ''
DEFINE_integer 'overlay' '' "Index of column with overlay levels (1-based)" 'O'
DEFINE_string 'aggr' '' "Aggregate function for overlay: <mean> <median>  <mean_j> <median_j>. Use mean_j | median_j to keep  density overlay but aggregate junction counts" 'A'
DEFINE_integer 'colorfactor' '' "Index of column with color levels (1-based)" 'C'
DEFINE_float 'alpha' '0.5' "Transparency level for density histogram" ''
DEFINE_string 'palette' '' "Color palette file. tsv file with >=1 columns, where  the color is the first column" 'P'
DEFINE_integer 'labels' '' "Index of column with labels (1-based)" 'L'
DEFINE_float 'height' '2' "Height of the individual signal plot in inches" ''
DEFINE_float 'annheight' '1.5' "Height of annotation plot in inches" ''
DEFINE_float 'width' '10' "Width of the plot in inches" ''
DEFINE_integer 'basesize' '14' "Base font size of the plot in pch" ''
DEFINE_string 'outformat' 'pdf' "Output file format: <pdf> <svg> <png> <jpeg> <tiff>" 'F'
DEFINE_integer 'outresolution' '300' "Output file resolution in PPI (pixels per inch).  Applies only to raster output formats" 'R'

DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode


printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# functional part of the script
flagsAsString=""
if [ "$FLAGS_shrink" == 0 ]; then
 flagsAsString="$flagsAsString --shrink  "
fi
if  [ ! -z "$FLAGS_help" ]; then
	flagsAsString="$flagsAsString --help  $FLAGS_help"
fi
if  [ ! -z "$FLAGS_bam" ]; then
	verifyFileExistence "$FLAGS_bam"
	flagsAsString="$flagsAsString -b  $FLAGS_bam"
fi
if  [ ! -z "$FLAGS_coordinates" ]; then
	flagsAsString="$flagsAsString -c  $FLAGS_coordinates"
fi
if  [ ! -z "$FLAGS_outprefix" ]; then
	flagsAsString="$flagsAsString -o  $FLAGS_outprefix"
fi
if  [ ! -z "$FLAGS_outstrand" ]; then
	flagsAsString="$flagsAsString -S  $FLAGS_outstrand"
fi
if  [ ! -z "$FLAGS_strand" ]; then
	flagsAsString="$flagsAsString -s  $FLAGS_strand"
fi
if  [ ! -z "$FLAGS_mincoverage" ]; then
	flagsAsString="$flagsAsString -M  $FLAGS_mincoverage"
fi
if  [ ! -z "$FLAGS_junctionsbed" ]; then
	flagsAsString="$flagsAsString -j  $FLAGS_junctionsbed"
fi
if  [ ! -z "$FLAGS_gtf" ]; then
	flagsAsString="$flagsAsString -g  $FLAGS_gtf"
fi
if  [ ! -z "$FLAGS_overlay" ]; then
	flagsAsString="$flagsAsString -O  $FLAGS_overlay"
fi
if  [ ! -z "$FLAGS_aggr" ]; then
	flagsAsString="$flagsAsString -A  $FLAGS_aggr"
fi
if  [ ! -z "$FLAGS_colorfactor" ]; then
	flagsAsString="$flagsAsString -C  $FLAGS_colorfactor"
fi
if  [ ! -z "$FLAGS_alpha" ]; then
	flagsAsString="$flagsAsString --alpha  $FLAGS_alpha"
fi
if  [ ! -z "$FLAGS_palette" ]; then
	flagsAsString="$flagsAsString -P  $FLAGS_palette"
fi
if  [ ! -z "$FLAGS_labels" ]; then
	flagsAsString="$flagsAsString -L  $FLAGS_labels"
fi
if  [ ! -z "$FLAGS_height" ]; then
	flagsAsString="$flagsAsString --height  $FLAGS_height"
fi
if  [ ! -z "$FLAGS_annheight" ]; then
	flagsAsString="$flagsAsString --ann-height  $FLAGS_annheight"
fi
if  [ ! -z "$FLAGS_width" ]; then
	flagsAsString="$flagsAsString --width  $FLAGS_width"
fi
if  [ ! -z "$FLAGS_basesize" ]; then
	flagsAsString="$flagsAsString --base-size  $FLAGS_basesize"
fi
if  [ ! -z "$FLAGS_outformat" ]; then
	flagsAsString="$flagsAsString -F  $FLAGS_outformat"
fi
if  [ ! -z "$FLAGS_outresolution" ]; then
	flagsAsString="$flagsAsString -R  $FLAGS_outresolution"
fi
# run it
MESSAGE=$(sashimi-plot.py $flagsAsString)

echo $flagsAsString;

RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile sashimiPlot)
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
