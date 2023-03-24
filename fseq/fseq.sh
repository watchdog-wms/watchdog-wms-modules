#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi
# define parameters
DEFINE_string 'dir' '' 'working directory' 'd'
DEFINE_string 'bam' '' 'bam file' 'b'
DEFINE_string 'name' '' 'name' 'n'
DEFINE_string 'pathToFseq' '' 'path to Fseq jar' 'p'
DEFINE_string 'heapSize' '-Xmx32000M' '[optional] choose JAVA OPTS heap size' 'j'
DEFINE_integer 'mergeDist' '0' '[optional] distance for merging' 'm'
DEFINE_integer 'minLen' '0' '[optional] min length for fragments' 'l'


JAVAOPTS="-Xmx32000M"
#DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
#DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
#DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' 'd'
# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

# check if mandatory arguments are there
if [ -z "$FLAGS_dir" ]; then
	echoError "Parameter -dir must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_bam" ]; then
	echoError "Parameter -bam must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_name" ]; then
	echoError "Parameter -name must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_pathToFseq" ]; then
	echoError "Parameter -pathToFseq must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT #####################################################

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# TODO: write functional part of the script here
echo $FLAGS_bam
# check, if the input files exist
if test -f "$FLAGS_bam" ; then
	echo "file exists"
fi
verifyFileExistence "$FLAGS_bam"

if [ ! -d "$FLAGS_dir/BED" ]; then
  mkdir $FLAGS_dir/BED;
	echo "made dir" $FLAGS_dir/BED
fi

if [ ! -d "$FLAGS_dir/BED/$FLAGS_name" ]; then
  mkdir $FLAGS_dir/BED/$FLAGS_name;
	echo "made dir"$FLAGS_dir/BED/$FLAGS_name
fi
#/home/proj/software/FSeq/F-seq/dist~/fseq/lib/fseq.jar
#/home/proj/projekte/sequencing/Illumina/ATAC_WT_HSV/
bedtools bamtobed -i $FLAGS_bam > $FLAGS_dir/BED/$FLAGS_name.bed

if [ -L "$FLAGS_pathToFseq" ]; then
	echo "symlink"
else
	REALDIR="`dirname "$FLAGS_pathToFseq"`"
	CLASSPATH="$REALDIR/commons-cli-1.1.jar:$FLAGS_pathToFseq"
	echo $CLASSPATH
fi

for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do
	echo $i
	ESCAPE="_"
	grep -P "$i\t" $FLAGS_dir/BED/$FLAGS_name.bed> $FLAGS_dir/BED/$FLAGS_name\_$i.bed;
	echo 'running : 'java $FLAGS_heapSize -cp $CLASSPATH edu.duke.igsp.gkde.Main -of bed -o $FLAGS_dir/BED/$FLAGS_name $FLAGS_dir/BED/$FLAGS_name$ESCAPE$i.bed
	java $FLAGS_heapSize -cp $CLASSPATH edu.duke.igsp.gkde.Main -of bed -o $FLAGS_dir/BED/$FLAGS_name $FLAGS_dir/BED/$FLAGS_name$ESCAPE$i.bed;
done
$SCRIPT_FOLDER/mergeBEDFiles.pl $FLAGS_dir/BED/$FLAGS_name/ $FLAGS_mergeDist $FLAGS_minLen;

rm -r $FLAGS_dir/BED/$FLAGS_name/
rm -r $FLAGS_dir/BED/$FLAGS_name.bed
rm -r $FLAGS_dir/BED/$FLAGS_name'_chr'*.bed

exit $EXIT_OK

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
