#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'out' '' 'path to output folder' 'o'
DEFINE_string 'inputregs' '' 'single region to compute' 'i'
DEFINE_string 'everyPos' '' 'count every Position of read' 'e'
DEFINE_string 'bams' '' 'path to bam files' 'b'
DEFINE_string 'strandness' '' '0 if unstranded, 1 if forward' 's'
DEFINE_string 'pattern' '' 'pattern for bams' 'p'
DEFINE_string 'pseudocount' '' 'pseudocount to subtract' 'z'
DEFINE_string 'numrandomizations' '' 'num of randomizations' 'q'
DEFINE_string 'sampleAnnotation' '' 'samples annotated to names' 'a'
DEFINE_boolean 'debug' 'true' '[optional] prints debug messages.' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'
DEFINE_boolean 'version' 'false' '[optional] print version' 'v'



# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"

# check if mandatory arguments are there
if [[ -z $FLAGS_out ]]; then
        echo "Parameter -o must be set. (see --help for details)";
        exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_inputregs ]] ; then
        echo "Parameter -i must be set. (see --help for details)";
        exit $EXIT_MISSING_ARGUMENTS
fi
#extend param checks!!!


if [[ -z $FLAGS_pseudocount ]] ; then
	FLAGS_pseudocount=1
fi
if [[ -z $FLAGS_numrandomizations ]] ; then
        FLAGS_numrandomizations=1000
fi
if [[ -z $FLAGS_everyPos ]] ; then
        FLAGS_everyPos="false"
fi

printParamValues "parameters"





IN=$FLAGS_inputregs
arrIN=(${IN//_/ })


chr=${arrIN[0]}
s=${arrIN[1]}
e=${arrIN[2]}
strand=${arrIN[3]}
if [[ -z "$strand" ]] ; then
	strand="+"
fi
echo "calling window: $chr $s $e $strand"
plotfile=$FLAGS_out"/"$chr"-"$s"-"$e"/amss_plots.pdf"
echo $plotfile
if [[ -f "$plotfile" ]] ; then
	echo "window already exists"
	exit $EXIT_OK
fi
	
start=`date +%s`
for b in $FLAGS_bams*; do
	f=$(basename "$b")
	if [[ $f =~ $FLAGS_pattern ]] ; then
		echo $f
		#fetch counts of bam file in the specified genomic region
		python3 $SCRIPT_FOLDER"/quantify_curves_difference.py" --chr $chr --start $s --end $e --givenstrand $strand --strandness $FLAGS_strandness --bam $FLAGS_bams/$f --out $FLAGS_out$chr"-"$s"-"$e"/" --everyPos $FLAGS_everyPos	#fetches specified region and counts from bam
	fi
done
runtime=$((`date +%s`-start))
echo "TIME select window counts "$runtime

#computes AMSS
#Rscript $SCRIPT_FOLDER"/quantify_curves_efficient.R" $chr $s $e $FLAGS_out $FLAGS_sampleAnnotation $FLAGS_pseudocount $FLAGS_numrandomizations
Rscript $SCRIPT_FOLDER"/quantify_curves_efficient.R" -c $chr -s $s -t $e -o $FLAGS_out -a $FLAGS_sampleAnnotation -p $FLAGS_pseudocount -n $FLAGS_numrandomizations
#rm -r $FLAGS_out$chr"-"$s"-"$e"/counts/"

#writes amss regions into annotation files for dexseq per single dir
Rscript $SCRIPT_FOLDER"/createSAF_quantCurves_efficient.R" $FLAGS_out"/"$chr"-"$s"-"$e"/" $FLAGS_sampleAnnotation $strand





writeParam2File "$FLAGS_returnFilePath" "out" "$FLAGS_out"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK


