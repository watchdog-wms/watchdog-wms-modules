#!/usr/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'annot' '' 'file name of annotation file' 'a'
DEFINE_string 'annot_fc' '' 'file name of annotation file for feature counts' 'f'
DEFINE_string 'indir' '' 'path to input directory - the windows files' 'i'
DEFINE_boolean 'debug' 'true' '[optional] prints out debug messages.' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'
DEFINE_boolean 'version' 'false' '[optional] print version' 'v'



# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"

# check if mandatory arguments are there
if [[ -z $FLAGS_annot ]]; then
        echo "Parameter -o must be set. (see --help for details)";
        exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_annot_fc ]] ; then
        echo "Parameter -i must be set. (see --help for details)";
        exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_indir ]] ; then
        echo "Parameter -i must be set. (see --help for details)";
        exit $EXIT_MISSING_ARGUMENTS
fi


printParamValues "parameters"

###### collect all amss files into one
all_amss=$FLAGS_indir"/amss.regions"
rm $all_amss
touch $all_amss
echo -e "start_seq_idx\tend_seq_idx\tsum_height\tAUC\tlen\tdir\twindow\n" > $all_amss
for single in $(find $FLAGS_indir -name "*final_amss*"); do
        dir=$(basename $single _final_amss.tsv)
        window=$(echo $single | rev | cut -d/ -f2 | rev)
        cat $single | grep -v "start" | grep -v "NA" | awk -v d=$dir -v w=$window '{print $0 "\t" d "\t" w}' >> $all_amss
done


mkdir -p "$(dirname "$FLAGS_annot")"


rm $FLAGS_annot_fc
touch $FLAGS_annot_fc

for single in $(find $FLAGS_indir -name "annotation_fc_filled_up*"); do
	cat $single >> $FLAGS_annot_fc
done

rm $FLAGS_annot
touch $FLAGS_annot

for single in $(find $FLAGS_indir -name "annotation_filled_up*"); do
        cat $single >> $FLAGS_annot
done




writeParam2File "$FLAGS_returnFilePath" "out" "$FLAGS_out"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK


