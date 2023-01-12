#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'outputDir' '' 'path to output folder' 'o'
DEFINE_string 'bam' '' 'path to bam file' 'b'
DEFINE_string 'gtf' '' 'path to gtf file' 'g'
DEFINE_string 'tss' '' 'path to tss file' 't'
DEFINE_string 'genelist' '' 'path to genes file' 'l'
DEFINE_integer 'promStart' '' 'start position of promoter window' 'a'
DEFINE_integer 'promEnd' '' 'end position of promoter window' 'k'
DEFINE_integer 'bodyStart' '' 'start position of body window' 'c'
DEFINE_integer 'bodyLength' '' 'end position of body window' 'l'
DEFINE_boolean 'debug' 'true' '[optional] prints out debug messages.' 'd'
DEFINE_string 'returnFilePath' '' 'path to the return variables file' 'r'
DEFINE_boolean 'version' 'false' '[optional] print version' 'v'


# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"

# check if mandatory arguments are there
if [[ -z $FLAGS_outputDir ]]; then
	echo "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_bam ]]; then
	echoError "Parameter -b must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_gtf ]]; thencalcPI.sh
	echo "Parameter -g must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
elif [[ -z $FLAGS_tss ]] ; then
	echo "Parameter -t must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

# default values
if [[ -z $FLAGS_promStart ]] ; then
	FLAGS_promStart=0
	FLAGS_promEnd=250
	FLAGS_bodyStart=250
	FLAGS_bodyLength=2000
fi

printParamValues "parameters (after default set if needed)"

### start script

bamname=$(basename $FLAGS_bam .bam)


print $FLAGS_promStart
print $FLAGS_promEnd
print $FLAGS_bodyStart
print $FLAGS_bodyEnd


currentdir=$(pwd)

### 1. create saf files: -one for promoter -one for body
createSaf=$currentdir"/../createBEDandSAF/createBEDandSAF.jar"
genelist=$FLAGS_genelist

bodyEnd=$(("$FLAGS_bodyStart"+"$FLAGS_bodyLength"))

java -jar $createSaf -tss $FLAGS_tss -gtf $FLAGS_gtf -genelist $genelist -saf -outdir $FLAGS_outputDir -name pro_uninf_prom_$bamname -ps $FLAGS_promStart -pe $FLAGS_promEnd -info
java -jar $createSaf -tss $FLAGS_tss -gtf $FLAGS_gtf -genelist $genelist -saf -outdir $FLAGS_outputDir -name pro_uninf_body_$bamname -bs $FLAGS_bodyStart -be $bodyEnd -info


cat $FLAGS_outputDir/mapping.info | cut -f 1,2,3,4,8 | awk '{split($6,coord,"-"); {print $4"\t"$1"\t"$3"\t"$2"\t"substr(coord[2], 1, length(coord[2]-1))-coord[1]+1"\t"coord[1]"\t"substr(coord[2], 1, length(coord[2]-1))}}' > $FLAGS_outputDir/$bamname"_tss2gene.length"


Rscript $currentdir/merge.helper $FLAGS_outputDir/pro_uninf_prom_$bamname.saf $FLAGS_outputDir/$bamname"_tss2gene.length" $FLAGS_outputDir/joined_prom_$bamname.file

promLen=$(($FLAGS_promEnd-$FLAGS_promStart))
totLen=$(($promLen+$FLAGS_bodyLength))

java -jar $currentdir/modifySafWindow.jar -info $FLAGS_outputDir/joined_prom_$bamname.file -start $FLAGS_promStart -end $FLAGS_promEnd -out $FLAGS_outputDir -total_length $totLen -promoter_length $promLen

mv $FLAGS_outputDir/modified_window.saf $FLAGS_outputDir/modified_window_prom_$bamname.saf

Rscript $currentdir/merge.helper $FLAGS_outputDir/pro_uninf_body_$bamname.saf $FLAGS_outputDir/$bamname"_tss2gene.length" $FLAGS_outputDir/joined_body_$bamname.file

bodyEnd=$(($FLAGS_bodyStart+$FLAGS_bodyLength))

java -jar $currentdir/modifySafWindow.jar -info $FLAGS_outputDir/joined_body_$bamname.file -start $FLAGS_bodyStart -end $bodyEnd -out $FLAGS_outputDir -total_length $totLen -promoter_length $promLen

mv $FLAGS_outputDir/modified_window.saf $FLAGS_outputDir/modified_window_body_$bamname.saf







### 2. feature counts for those windows 2 times for promoter and body into 2 diff count files
featurecounts=$currentdir"/../../modules/featureCounts/featureCounts.sh"
#parameter minReadOverlap geht nur mit subread version 1.4

$featurecounts --annotation $FLAGS_outputDir/modified_window_prom_$bamname.saf --input $FLAGS_bam --output $FLAGS_outputDir/prom_$bamname.counts --annotationType SAF --multiMapping --multiCountMetaFeature --minReadOverlap 25 --stranded 2 --moduleVersion 1

$featurecounts --annotation $FLAGS_outputDir/modified_window_body_$bamname.saf --input $FLAGS_bam --output $FLAGS_outputDir/body_$bamname.counts --annotationType SAF --multiMapping --multiCountMetaFeature --minReadOverlap 25 --stranded 2 --moduleVersion 1

rm -f $FLAGS_outputDir/*png
rm -f $FLAGS_outputDir/*stats

rm $FLAGS_outputDir/joined_prom_$bamname.file
rm $FLAGS_outputDir/joined_body_$bamname.file
rm $FLAGS_outputDir/$bamname"_tss2gene.length"






### 3. rpkm aus gene body und promoter counts berechnen
cat $FLAGS_outputDir/prom_$bamname.counts | awk -v pl="$promLen" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/pl}' > $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_$bamname.counts | awk -v bl="$FLAGS_bodyLength" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/bl}' > $FLAGS_outputDir/body_normed_$bamname.counts

cat $FLAGS_outputDir/prom_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/body_normed_$bamname.counts

genes=($(awk '{print $1}' $FLAGS_outputDir/modified_window_prom_$bamname.saf))

totalMappedReads=$(cat $FLAGS_outputDir/prom_$bamname.counts.summary | awk '{s+=$2} END {print s}')

echo -e "gene\tprom-rpkm" > $FLAGS_outputDir/gene_prom_$bamname.rpkm

for gene in "${genes[@]}" ; do
	readcounts=$(cat $FLAGS_outputDir/prom_$bamname.counts | grep "$gene" | cut -f 7)
        readlength=$(cat $FLAGS_outputDir/prom_$bamname.counts | grep "$gene" | cut -f 6)
	rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" 'BEGIN {print ""1000000000*rc/rl/tmr}')
	echo -e "$gene\t$rpkm" >> $FLAGS_outputDir/gene_prom_$bamname.rpkm
done


totalMappedReads=$(cat $FLAGS_outputDir/body_$bamname.counts.summary | awk '{s+=$2} END {print s}')

echo -e "gene\tbody-rpkm" > $FLAGS_outputDir/gene_body_$bamname.rpkm

for gene in ${genes[@]} ; do
        readcounts=$(cat $FLAGS_outputDir/body_$bamname.counts | grep "$gene" | cut -f 7)
        readlength=$(cat $FLAGS_outputDir/body_$bamname.counts | grep "$gene" | cut -f 6)
	rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" 'BEGIN {print ""1000000000*rc/rl/tmr}')
        echo -e "$gene\t$rpkm" >> $FLAGS_outputDir/gene_body_$bamname.rpkm
done

join <(sort $FLAGS_outputDir/gene_prom_$bamname.rpkm) <(sort $FLAGS_outputDir/gene_body_$bamname.rpkm) -t $'\t' > $FLAGS_outputDir/gene_sense_$bamname.rpkm







### 4. pausing index als ratio berechnen
promLen=$(("$FLAGS_promEnd"-"$FLAGS_promStart"))

cat $FLAGS_outputDir/prom_$bamname.counts | cut -f 1,7 | awk -v pl="$promLen" '{print $1"\t"$2/pl}' > $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_$bamname.counts | cut -f 1,7 | awk -v bl="$FLAGS_bodyLength" '{print $1"\t"$2/bl}' > $FLAGS_outputDir/body_normed_$bamname.counts

cat $FLAGS_outputDir/prom_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/body_normed_$bamname.counts

join <(sort $FLAGS_outputDir/prom_normed_$bamname.counts) <(sort $FLAGS_outputDir/body_normed_$bamname.counts) -t $'\t' > $FLAGS_outputDir/gene_$bamname.counts

echo -e "gene\tpromoter\tbody\tPI" > $FLAGS_outputDir/gene_$bamname.pi

while read line ; do
	if [[ $line == *"#"* ]] || [[ $line == *"Gene"* ]] ; then
		continue
	fi
	array=($line)
	body="${array[2]}"
	if [[ $body == "0" ]] ; then
		echo -e "$line\tNA" >> $FLAGS_outputDir/gene_$bamname.pi
	else
		prom="${array[1]}"
		pi=$(awk "BEGIN {printf \"%.2f\",${prom}/${body}}")
		echo -e "$line\t$pi" >> $FLAGS_outputDir/gene_$bamname.pi
	fi
done < $FLAGS_outputDir/gene_$bamname.counts

echo -e "gene\tpromoter\tbody\tPI" > $FLAGS_outputDir/genes_$bamname.pi
cat $FLAGS_outputDir/gene_$bamname.pi | tail -n +2 | sort -k 4,4 -n -r >> $FLAGS_outputDir/genes_$bamname.pi
rm $FLAGS_outputDir/gene_$bamname.pi

cat $FLAGS_outputDir/genes_$bamname.pi | tail -n +2 > $FLAGS_outputDir/tmp.txt && mv $FLAGS_outputDir/tmp.txt $FLAGS_outputDir/genes_$bamname.pi








### 6. result of this modul
# - per gene: sense PI, antisense PI, promoter RPKM, body RPKM, antisense promoter RPKM, antisense body RPKM
# join all middle files in one simply

echo -e "gene\tprom_count\tbody_count\tPI\tprom_rpkm\tbody_rpkm" > $FLAGS_outputDir/result_$bamname.sense
join <(sort $FLAGS_outputDir/genes_$bamname.pi) <(sort $FLAGS_outputDir/gene_sense_$bamname.rpkm) -t $'\t' >> $FLAGS_outputDir/result_$bamname.sense
sort -k4 -r -o $FLAGS_outputDir/result_$bamname.sense $FLAGS_outputDir/result_$bamname.sense
sed -i s/" "/"\t"/g $FLAGS_outputDir/result_$bamname.sense

Rscript $currentdir/classify.R $FLAGS_outputDir/result_$bamname.sense $FLAGS_outputDir

mv $FLAGS_outputDir/result_sense.table $FLAGS_outputDir/$bamname.result

rm $FLAGS_outputDir/mapping*






writeParam2File "$FLAGS_returnFilePath" "pausingindices" "$FLAGS_outputDir"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK
