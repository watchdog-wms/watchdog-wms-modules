#!/bin/bash

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh $@

# define parameters
DEFINE_string 'outputDir' '' 'path to output folder' 'o'
DEFINE_string 'bam' '' 'path to bam file' 'b'
DEFINE_string 'gtf' '' 'path to gtf file' 'g'
DEFINE_string 'tss' '' 'path to tss file' 't'
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
elif [[ -z $FLAGS_gtf ]]; then
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





### 1. create saf files: -one for promoter -one for body

# !!!!!!!!! adapt window search to gene length -> in sense and antisense maybe wit hinterval tree?

#createSaf=/home/proj/software/watchdog/elena_modules/createBEDandSAF/createBEDandSAF.jar
#genelist=/home/proj/projekte/sequencing/Illumina/PolIIPausing/TSS-profiling/workflow/pro_8cluster_onlyuninf.txt

bodyEnd=$(("$FLAGS_bodyStart"+"$FLAGS_bodyLength"))
promLen=$(($FLAGS_promEnd-$FLAGS_promStart))
totLen=$(($promLen+$FLAGS_bodyLength))


#mv $FLAGS_outputDir/modified_window.saf $FLAGS_outputDir/modified_window_prom_$bamname.saf
#mv $FLAGS_outputDir/modified_window.saf $FLAGS_outputDir/modified_window_body_$bamname.saf



### 2. feature counts for those windows 2 times for promoter and body into 2 diff count files

featurecounts=/home/proj/software/watchdog/modules/featureCounts/featureCounts.sh

#parameter minReadOverlap geht nur mit subread version 1.4

$featurecounts --annotation $FLAGS_outputDir/cdk11_prom.saf --input $FLAGS_bam --output $FLAGS_outputDir/prom_$bamname.counts --annotationType SAF --multiMapping --multiCountMetaFeature --minReadOverlap 25 --stranded 2 --moduleVersion 1

$featurecounts --annotation $FLAGS_outputDir/cdk11_body.saf --input $FLAGS_bam --output $FLAGS_outputDir/body_$bamname.counts --annotationType SAF --multiMapping --multiCountMetaFeature --minReadOverlap 25 --stranded 2 --moduleVersion 1

#rm -f $FLAGS_outputDir/*png
#rm -f $FLAGS_outputDir/*stats

#rm $FLAGS_outputDir/joined_prom_$bamname.file
#rm $FLAGS_outputDir/joined_body_$bamname.file
#rm $FLAGS_outputDir/pro_uninf_prom_$bamname.saf
#rm $FLAGS_outputDir/pro_uninf_body_$bamname.saf
#rm $FLAGS_outputDir/$bamname"_tss2gene.length"






### 3. rpkm aus gene body und promoter counts berechnen

# - group genes in expressed + paused, expressed + not paused, not expressed + paused, not expressed + not paused

cat $FLAGS_outputDir/prom_$bamname.counts | awk -v pl="$promLen" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/pl}' > $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_$bamname.counts | awk -v bl="$FLAGS_bodyLength" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/bl}' > $FLAGS_outputDir/body_normed_$bamname.counts

cat $FLAGS_outputDir/prom_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/body_normed_$bamname.counts

genes=($(awk '{print $1}' $FLAGS_outputDir/cdk11_prom.saf))
#echo ${#genes[@]}

totalMappedReads=$(cat $FLAGS_outputDir/prom_$bamname.counts.summary | awk '{s+=$2} END {print s}')

echo -e "gene\tprom-rpkm" > $FLAGS_outputDir/gene_prom_$bamname.rpkm

for gene in "${genes[@]}" ; do
	#echo $gene
	#readcounts=$(cat $FLAGS_outputDir/prom_normed_$bamname.counts | grep "$gene" | cut -f 7)
	#readlength=$(cat $FLAGS_outputDir/prom_normed_$bamname.counts | grep "$gene" | cut -f 6)
	readcounts=$(cat $FLAGS_outputDir/prom_$bamname.counts | grep "$gene" | cut -f 7)
        readlength=$(cat $FLAGS_outputDir/prom_$bamname.counts | grep "$gene" | cut -f 6)
	#rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" '{x=1000000000*rc/rl/tmr  {{print x}}')
	rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" 'BEGIN {print ""1000000000*rc/rl/tmr}')
	#rpkm=$(( 1000000000*$readcounts/$readlength/$totalMappedReads ))
	##echo -e "$readcounts\t$readlength\t$totalMappedReads\t$rpkm"
	#echo $rpkm
	echo -e "$gene\t$rpkm" >> $FLAGS_outputDir/gene_prom_$bamname.rpkm
done


totalMappedReads=$(cat $FLAGS_outputDir/body_$bamname.counts.summary | awk '{s+=$2} END {print s}')

echo -e "gene\tbody-rpkm" > $FLAGS_outputDir/gene_body_$bamname.rpkm

for gene in ${genes[@]} ; do
        readcounts=$(cat $FLAGS_outputDir/body_$bamname.counts | grep "$gene" | cut -f 7)
        readlength=$(cat $FLAGS_outputDir/body_$bamname.counts | grep "$gene" | cut -f 6)
        #rpkm=$(($readcounts*1000000000/$readlength/$totalMappedReads))
	#rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" '{print 1000000000*rc/rl/tmr}')
	rpkm=$(awk -v rc="$readcounts" -v rl="$readlength" -v tmr="$totalMappedReads" 'BEGIN {print ""1000000000*rc/rl/tmr}')
        echo -e "$gene\t$rpkm" >> $FLAGS_outputDir/gene_body_$bamname.rpkm
done

join <(sort $FLAGS_outputDir/gene_prom_$bamname.rpkm) <(sort $FLAGS_outputDir/gene_body_$bamname.rpkm) -t $'\t' > $FLAGS_outputDir/gene_sense_$bamname.rpkm
#rm $FLAGS_outputDir/gene_prom_$bamname.rpkm
#rm $FLAGS_outputDir/gene_body_$bamname.rpkm







### 4. pausing index als ratio berechnen

promLen=$(("$FLAGS_promEnd"-"$FLAGS_promStart"))

cat $FLAGS_outputDir/prom_$bamname.counts | cut -f 1,7 | awk -v pl="$promLen" '{print $1"\t"$2/pl}' > $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_$bamname.counts | cut -f 1,7 | awk -v bl="$FLAGS_bodyLength" '{print $1"\t"$2/bl}' > $FLAGS_outputDir/body_normed_$bamname.counts

cat $FLAGS_outputDir/prom_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/prom_normed_$bamname.counts
cat $FLAGS_outputDir/body_normed_$bamname.counts | tail -n +3 > $FLAGS_outputDir/tmp.file && mv $FLAGS_outputDir/tmp.file $FLAGS_outputDir/body_normed_$bamname.counts

#join $FLAGS_outputDir/prom_normed.counts $FLAGS_outputDir/body_normed.counts | awk '{print $1"\t"$2"\t"$3"\t"$2/$3}' > $FLAGS_outputDir/gene.counts
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
		#echo "count is zero"
	else
		prom="${array[1]}"
		pi=$(awk "BEGIN {printf \"%.2f\",${prom}/${body}}")
		echo -e "$line\t$pi" >> $FLAGS_outputDir/gene_$bamname.pi
	fi
done < $FLAGS_outputDir/gene_$bamname.counts

echo -e "gene\tpromoter\tbody\tPI" > $FLAGS_outputDir/genes_$bamname.pi
cat $FLAGS_outputDir/gene_$bamname.pi | tail -n +2 | sort -k 4,4 -n -r >> $FLAGS_outputDir/genes_$bamname.pi
rm $FLAGS_outputDir/gene_$bamname.pi
#rm $FLAGS_outputDir/gene_$bamname.counts

cat $FLAGS_outputDir/genes_$bamname.pi | tail -n +2 > $FLAGS_outputDir/tmp.txt && mv $FLAGS_outputDir/tmp.txt $FLAGS_outputDir/genes_$bamname.pi







### 5. antisense pi berechnen

# antisense PI gespiegelt, also promoter Bereich -1 bis -250, gen body -250 bis -2250, mit mindestabstand x zu nächstem gen
# - wenn gen früher anfängt is body kleiner
# - wenn nächstes gen im promoter ist ohne mindestabstand, also nicht mind. 251+x weit weg, dann PI = undefiniert
# - interval tree

#java -jar $createSaf -tss $FLAGS_tss -gtf $FLAGS_gtf -genelist $genelist -saf -outdir $FLAGS_outputDir -name pro_uninf_prom -ps $FLAGS_promStart -pe $FLAGS_promEnd -coordinates
#java -jar $createSaf -tss $FLAGS_tss -gtf $FLAGS_gtf -genelist $genelist -saf -outdir $FLAGS_outputDir -name pro_uninf_body -bs $FLAGS_bodyStart -be $bodyEnd -coordinates

#$featurecounts --annotation $FLAGS_outputDir/pro_uninf_prom.saf --input $FLAGS_bam --output $FLAGS_outputDir/prom.counts --annotationType SAF --multiMapping true -- multiCountMetaFeature true -- minReadOverlap 25 --stranded 2
#$featurecounts --annotation $FLAGS_outputDir/pro_uninf_body.saf --input $FLAGS_bam --output $FLAGS_outputDir/body.counts --annotationType SAF --multiMapping true -- multiCountMetaFeature true -- minReadOverlap 25 --stranded 2

#rm $FLAGS_outputDir/*png
#rm $FLAGS_outputDir/*stats







### 6. result of this modul

# - per gene: sense PI, antisense PI, promoter RPKM, body RPKM, antisense promoter RPKM, antisense body RPKM
# join all middle files in one simply

echo -e "gene\tprom_count\tbody_count\tPI\tprom_rpkm\tbody_rpkm" > $FLAGS_outputDir/result_$bamname.sense
join <(sort $FLAGS_outputDir/genes_$bamname.pi) <(sort $FLAGS_outputDir/gene_sense_$bamname.rpkm) -t $'\t' >> $FLAGS_outputDir/result_$bamname.sense
sort -k4 -r -o $FLAGS_outputDir/result_$bamname.sense $FLAGS_outputDir/result_$bamname.sense
sed -i s/" "/"\t"/g $FLAGS_outputDir/result_$bamname.sense

#tail -n +2 "$FLAGS_outputDir/"result_"$bamname.sense" > "$FLAGS_outputDir/res.tmp" && mv "$FLAGS_outputDir/res.tmp" "$FLAGS_outputDir/"result_"$bamname.sense"

#rm $FLAGS_outputDir/gene_sense_$bamname.rpkm
#rm $FLAGS_outputDir/genes_$bamname.pi


Rscript /home/proj/software/watchdog/elena_modules/pausingIndex/classify.R $FLAGS_outputDir/result_$bamname.sense $FLAGS_outputDir
#rm $FLAGS_outputDir/result_$bamname.sense


#bamname=$(basename "$FLAGS_bam" .bam)
mv $FLAGS_outputDir/result_sense.table $FLAGS_outputDir/$bamname.result

#rm $FLAGS_outputDir/body*
#rm $FLAGS_outputDir/prom*
#rm $FLAGS_outputDir/modified*
#rm -r $FLAGS_outputDir/logfiles
#rm $FLAGS_outputDir/mapping*


writeParam2File "$FLAGS_returnFilePath" "pausingindices" "$FLAGS_outputDir"
blockUntilFileIsWritten "$FLAGS_returnFilePath"

exit $EXIT_OK
