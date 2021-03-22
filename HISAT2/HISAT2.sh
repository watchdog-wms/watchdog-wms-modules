#!/bin/bash
## created using moduleMaker

SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='echo:printf:rm:grep:head:wc:hisat2'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_string 'index' '' "Index filename prefix (minus trailing .X.ht2)" 'x'
DEFINE_string 'paired1' '' "Files with #1 mates, paired with files in <m2>. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." '1'
DEFINE_string 'paired2' '' "Files with #2 mates, paired with files in <m1>. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." '2'
DEFINE_string 'output' '' "File for SAM output" 'S'
DEFINE_string 'unpaired' '' "Files with unpaired reads. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." 'U'
DEFINE_boolean 'fastq' '1' "query input files are FASTQ .fq/.fastq (default)" 'q'
DEFINE_boolean 'qseq' '1' "query input files are in Illumina's qseq format" ''
DEFINE_boolean 'fasta' '1' "query input files are (multi-)FASTA .fa/.mfa" 'f'
DEFINE_boolean 'raw' '1' "query input files are raw one-sequence-per-line" 'r'
DEFINE_boolean 'c' '1' "paired1, paired2, unpaired are sequences themselves, not files" 'c'
DEFINE_integer 's' '' "skip the first <int> reads/pairs in the input (none)" 's'
DEFINE_integer 'u' '' "stop after first <int> reads/pairs (no limit)" 'u'
DEFINE_string 'trim5' '' "trim <int> bases from 5'/left end of reads (0)" '5'
DEFINE_string 'trim3' '' "trim <int> bases from 3'/right end of reads (0)" '3'
DEFINE_boolean 'phred33' '1' "qualities are Phred+33 (default)" ''
DEFINE_boolean 'phred64' '1' "qualities are Phred+64" ''
DEFINE_boolean 'intquals' '1' "qualities encoded as space-delimited integers" ''
DEFINE_string 'nceil' '' "func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)" ''
DEFINE_boolean 'ignorequals' '1' "treat all quality values as 30 on Phred scale (off)" ''
DEFINE_boolean 'nofw' '1' "do not align forward (original) version of read (off)" ''
DEFINE_boolean 'norc' '1' "do not align reverse-complement version of read (off)" ''
DEFINE_integer 'pencansplice' '' "penalty for a canonical splice site (0)" ''
DEFINE_integer 'pennoncansplice' '' "penalty for a non-canonical splice site (12)" ''
DEFINE_string 'pencanintronlen' '' "penalty for long introns (G,-8,1) with canonical splice sites" ''
DEFINE_string 'pennoncanintronlen' '' "penalty for long introns (G,-8,1) with noncanonical splice sites" ''
DEFINE_integer 'minintronlen' '' "minimum intron length (20)" ''
DEFINE_integer 'maxintronlen' '' "maximum intron length (500000)" ''
DEFINE_string 'knownsplicesiteinfile' '' "provide a list of known splice sites" ''
DEFINE_string 'novelsplicesiteoutfile' '' "report a list of splice sites" ''
DEFINE_string 'novelsplicesiteinfile' '' "provide a list of novel splice sites" ''
DEFINE_boolean 'notempsplicesite' '1' "disable the use of splice sites found" ''
DEFINE_boolean 'nosplicedalignment' '1' "disable spliced alignment" ''
DEFINE_string 'rnastrandness' '' "Specify strand-specific information (unstranded)" ''
DEFINE_boolean 'tmo' '1' "Reports only those alignments within known transcriptome" ''
DEFINE_boolean 'dta' '1' "Reports alignments tailored for transcript assemblers" ''
DEFINE_boolean 'dtacufflinks' '1' "Reports alignments tailored specifically for cufflinks" ''
DEFINE_integer 'ma' '' "match bonus (0 for --end-to-end, 2 for --local)" ''
DEFINE_string 'mp' '' "max and min penalties for mismatch; lower qual = lower penalty <6,2>" ''
DEFINE_string 'sp' '' "max and min penalties for soft-clipping; lower qual = lower penalty <2,1>" ''
DEFINE_integer 'np' '' "penalty for non-A/C/G/Ts in read/ref (1)" ''
DEFINE_string 'rdg' '' "read gap open, extend penalties (5,3)" ''
DEFINE_string 'rfg' '' "reference gap open, extend penalties (5,3)" ''
DEFINE_string 'scoremin' '' "min acceptable alignment score w/r/t read length (L,0.0,-0.2)" ''
DEFINE_integer 'k' '' "report up to <int> alns per read; MAPQ not meaningful" 'k'
DEFINE_integer 'a' '' "report all alignments; very slow, MAPQ not meaningful" 'a'
DEFINE_boolean 'fr' '1' "-1, -2 mates align fw/rev" ''
DEFINE_boolean 'rf' '1' "-1, -2 mates align rev/fw" ''
DEFINE_boolean 'ff' '1' "-1, -2 mates align fw/fw" ''
DEFINE_boolean 'nomixed' '1' "suppress unpaired alignments for paired reads" ''
DEFINE_boolean 'nodiscordant' '1' "suppress discordant alignments for paired reads" ''
DEFINE_boolean 't' '1' "print wall-clock time taken by search phases" 't'
DEFINE_string 'un' '' "write unpaired reads that didn't align to <path>" ''
DEFINE_string 'al' '' "write unpaired reads that aligned at least once to <path>" ''
DEFINE_string 'unconc' '' "write pairs that didn't align concordantly to <path>" ''
DEFINE_string 'alconc' '' "write pairs that aligned concordantly at least once to <path>" ''
DEFINE_boolean 'quiet' '1' "print nothing to stderr except serious errors" ''
DEFINE_string 'metfile' '' "send metrics to file at <path> (off)" ''
DEFINE_boolean 'metstderr' '1' "send metrics to stderr (off)" ''
DEFINE_integer 'met' '' "report internal counters & metrics every <int> secs (1)" ''
DEFINE_boolean 'nohead' '1' "supppress header lines, i.e. lines starting with @" ''
DEFINE_boolean 'nosq' '1' "supppress @SQ header lines" ''
DEFINE_string 'rgid' '' "set read group id, reflected in @RG line and RG:Z: opt field" ''
DEFINE_string 'rg' '' "add <text> (\"lab:value\") to @RG line of SAM header." ''
DEFINE_boolean 'omitsecseq' '1' "put '*' in SEQ and QUAL fields for secondary alignments." ''
DEFINE_integer 'offrate' '' "override offrate of index; must be >= index's offrate" 'o'
DEFINE_integer 'threads' '' "number of alignment threads to launch (1)" 'p'
DEFINE_boolean 'reorder' '1' "force SAM output order to match order of input reads" ''
DEFINE_boolean 'mm' '1' "use memory-mapped I/O for index; many 'bowtie's can share" ''
DEFINE_boolean 'qcfilter' '1' "filter out reads that are bad according to QSEQ filter" ''
DEFINE_integer 'seed' '' "seed for random number generator (0)" ''
DEFINE_boolean 'nondeterministic' '1' "seed rand. gen. arbitrarily instead of using read attributes" ''
DEFINE_boolean 'removechrname' '1' "remove 'chr' from reference names in alignment" ''
DEFINE_boolean 'addchrname' '1' "add 'chr' to reference names in alignment" ''
DEFINE_boolean 'version' '1' "print version information and quit" ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
        
        MESSAGE=$(hisat2 --version | head -n 1 | awk '{ print $3}');

        echo $MESSAGE
        exit $EXIT_OK
fi


# check if mandatory arguments are there
if [ -z "$FLAGS_index" ]; then 
	echoError "Parameter index must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_output" ]; then 
	echoError "Parameter output must be set. (see --help for details)"; 
	exit $EXIT_MISSING_ARGUMENTS
fi



printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# functional part of the script
flagsAsString=""
if [ "$FLAGS_fastq" == 0 ]; then
 flagsAsString="$flagsAsString -q  "
fi
if [ "$FLAGS_qseq" == 0 ]; then
 flagsAsString="$flagsAsString --qseq "
fi
if [ "$FLAGS_fasta" == 0 ]; then
 flagsAsString="$flagsAsString -f  "
fi
if [ "$FLAGS_raw" == 0 ]; then
 flagsAsString="$flagsAsString -r  "
fi
if [ "$FLAGS_c" == 0 ]; then
 flagsAsString="$flagsAsString -c  "
fi
if [ "$FLAGS_phred33" == 0 ]; then
 flagsAsString="$flagsAsString --phred33 "
fi
if [ "$FLAGS_phred64" == 0 ]; then
 flagsAsString="$flagsAsString --phred64 "
fi
if [ "$FLAGS_intquals" == 0 ]; then
 flagsAsString="$flagsAsString --int-quals "
fi
if [ "$FLAGS_ignorequals" == 0 ]; then
 flagsAsString="$flagsAsString --ignore-quals "
fi
if [ "$FLAGS_nofw" == 0 ]; then
 flagsAsString="$flagsAsString --nofw "
fi
if [ "$FLAGS_norc" == 0 ]; then
 flagsAsString="$flagsAsString --norc "
fi
if [ "$FLAGS_notempsplicesite" == 0 ]; then
 flagsAsString="$flagsAsString --no-temp-splicesite "
fi
if [ "$FLAGS_nosplicedalignment" == 0 ]; then
 flagsAsString="$flagsAsString --no-spliced-alignment "
fi
if [ "$FLAGS_tmo" == 0 ]; then
 flagsAsString="$flagsAsString --tmo "
fi
if [ "$FLAGS_dta" == 0 ]; then
 flagsAsString="$flagsAsString --dta "
fi
if [ "$FLAGS_dtacufflinks" == 0 ]; then
 flagsAsString="$flagsAsString --dta-cufflinks "
fi
if [ "$FLAGS_fr" == 0 ]; then
 flagsAsString="$flagsAsString --fr "
fi
if [ "$FLAGS_nomixed" == 0 ]; then
 flagsAsString="$flagsAsString --no-mixed "
fi
if [ "$FLAGS_nodiscordant" == 0 ]; then
 flagsAsString="$flagsAsString --no-discordant "
fi
if [ "$FLAGS_t" == 0 ]; then
 flagsAsString="$flagsAsString -t  "
fi
if [ "$FLAGS_quiet" == 0 ]; then
 flagsAsString="$flagsAsString --quiet "
fi
if [ "$FLAGS_metstderr" == 0 ]; then
 flagsAsString="$flagsAsString --met-stderr "
fi
if [ "$FLAGS_nohead" == 0 ]; then
 flagsAsString="$flagsAsString --no-head "
fi
if [ "$FLAGS_nosq" == 0 ]; then
 flagsAsString="$flagsAsString --no-sq "
fi
if [ "$FLAGS_omitsecseq" == 0 ]; then
 flagsAsString="$flagsAsString --omit-sec-seq "
fi
if [ "$FLAGS_reorder" == 0 ]; then
 flagsAsString="$flagsAsString --reorder "
fi
if [ "$FLAGS_mm" == 0 ]; then
 flagsAsString="$flagsAsString --mm "
fi
if [ "$FLAGS_qcfilter" == 0 ]; then
 flagsAsString="$flagsAsString --qc-filter "
fi
if [ "$FLAGS_nondeterministic" == 0 ]; then
 flagsAsString="$flagsAsString --non-deterministic "
fi
if [ "$FLAGS_removechrname" == 0 ]; then
 flagsAsString="$flagsAsString --remove-chrname "
fi
if [ "$FLAGS_addchrname" == 0 ]; then
 flagsAsString="$flagsAsString --add-chrname "
fi
if [ "$FLAGS_rf" == 0 ]; then
 flagsAsString="$flagsAsString - "
fi
if [ "$FLAGS_ff" == 0 ]; then
 flagsAsString="$flagsAsString  "
fi
if  [ ! -z "$FLAGS_unpaired" ]; then
	flagsAsString="$flagsAsString -U  $FLAGS_unpaired"
fi
if  [ ! -z "$FLAGS_s" ]; then
	flagsAsString="$flagsAsString -s  $FLAGS_s"
fi
if  [ ! -z "$FLAGS_u" ]; then
	flagsAsString="$flagsAsString -u  $FLAGS_u"
fi
if  [ ! -z "$FLAGS_5" ]; then
	flagsAsString="$flagsAsString -5  $FLAGS_5"
fi
if  [ ! -z "$FLAGS_3" ]; then
	flagsAsString="$flagsAsString -3  $FLAGS_3"
fi
if  [ ! -z "$FLAGS_nceil" ]; then
	flagsAsString="$flagsAsString --n-ceil $FLAGS_nceil"
fi
if  [ ! -z "$FLAGS_pencansplice" ]; then
	flagsAsString="$flagsAsString --pen-cansplice $FLAGS_pencansplice"
fi
if  [ ! -z "$FLAGS_pennoncansplice" ]; then
	flagsAsString="$flagsAsString --pen-noncansplice $FLAGS_pennoncansplice"
fi
if  [ ! -z "$FLAGS_pencanintronlen" ]; then
	flagsAsString="$flagsAsString --pen-canintronlen $FLAGS_pencanintronlen"
fi
if  [ ! -z "$FLAGS_pennoncanintronlen" ]; then
	flagsAsString="$flagsAsString --pen-noncanintronlen $FLAGS_pennoncanintronlen"
fi
if  [ ! -z "$FLAGS_minintronlen" ]; then
	flagsAsString="$flagsAsString --min-intronlen $FLAGS_minintronlen"
fi
if  [ ! -z "$FLAGS_maxintronlen" ]; then
	flagsAsString="$flagsAsString --max-intronlen $FLAGS_maxintronlen"
fi
if  [ ! -z "$FLAGS_knownsplicesiteinfile" ]; then
	flagsAsString="$flagsAsString --known-splicesite-infile $FLAGS_knownsplicesiteinfile"
fi
if  [ ! -z "$FLAGS_novelsplicesiteinfile" ]; then
	flagsAsString="$flagsAsString --novel-splicesite-infile $FLAGS_novelsplicesiteinfile"
fi
if  [ ! -z "$FLAGS_novelsplicesiteoutfile" ]; then
	flagsAsString="$flagsAsString --novel-splicesite-outfile $FLAGS_novelsplicesiteoutfile"
fi
if  [ ! -z "$FLAGS_rnastrandness" ]; then
	flagsAsString="$flagsAsString --rna-strandness $FLAGS_rnastrandness"
fi
if  [ ! -z "$FLAGS_ma" ]; then
	flagsAsString="$flagsAsString --ma $FLAGS_ma"
fi
if  [ ! -z "$FLAGS_mp" ]; then
	flagsAsString="$flagsAsString --mp $FLAGS_mp"
fi
if  [ ! -z "$FLAGS_sp" ]; then
	flagsAsString="$flagsAsString --sp $FLAGS_sp"
fi
if  [ ! -z "$FLAGS_np" ]; then
	flagsAsString="$flagsAsString --np $FLAGS_np"
fi
if  [ ! -z "$FLAGS_rdg" ]; then
	flagsAsString="$flagsAsString --rdg $FLAGS_rdg"
fi
if  [ ! -z "$FLAGS_rfg" ]; then
	flagsAsString="$flagsAsString --rfg $FLAGS_rfg"
fi
if  [ ! -z "$FLAGS_scoremin" ]; then
	flagsAsString="$flagsAsString --score-min $FLAGS_scoremin"
fi
if  [ ! -z "$FLAGS_k" ]; then
	flagsAsString="$flagsAsString -k  $FLAGS_k"
fi
if  [ ! -z "$FLAGS_a" ]; then
	flagsAsString="$flagsAsString -a  $FLAGS_a"
fi
if  [ ! -z "$FLAGS_un" ]; then
	flagsAsString="$flagsAsString --un $FLAGS_un"
fi
if  [ ! -z "$FLAGS_al" ]; then
	flagsAsString="$flagsAsString --al $FLAGS_al"
fi
if  [ ! -z "$FLAGS_unconc" ]; then
	flagsAsString="$flagsAsString --un-conc $FLAGS_unconc"
fi
if  [ ! -z "$FLAGS_alconc" ]; then
	flagsAsString="$flagsAsString --al-conc $FLAGS_alconc"
fi
if  [ ! -z "$FLAGS_ungz" ]; then
	flagsAsString="$flagsAsString --un-gz $FLAGS_ungz"
fi
if  [ ! -z "$FLAGS_metfile" ]; then
	flagsAsString="$flagsAsString --met-file $FLAGS_metfile"
fi
if  [ ! -z "$FLAGS_met" ]; then
	flagsAsString="$flagsAsString --met $FLAGS_met"
fi
if  [ ! -z "$FLAGS_rgid" ]; then
	flagsAsString="$flagsAsString --rg-id $FLAGS_rgid"
fi
if  [ ! -z "$FLAGS_rg" ]; then
	flagsAsString="$flagsAsString --rg $FLAGS_rg"
fi
if  [ ! -z "$FLAGS_offrate" ]; then
	flagsAsString="$flagsAsString -o  $FLAGS_offrate"
fi
if  [ ! -z "$FLAGS_threads" ]; then
	flagsAsString="$flagsAsString -p  $FLAGS_threads"
fi
if  [ ! -z "$FLAGS_seed" ]; then
	flagsAsString="$flagsAsString --seed $FLAGS_seed"
fi
if  [ ! -z "$FLAGS_index" ]; then
	flagsAsString="$flagsAsString -x  $FLAGS_index"
fi
if  [ ! -z "$FLAGS_paired1" ]; then
	flagsAsString="$flagsAsString -1  $FLAGS_paired1"
fi
if  [ ! -z "$FLAGS_paired2" ]; then
	flagsAsString="$flagsAsString -2  $FLAGS_paired2"
fi
if  [ ! -z "$FLAGS_output" ]; then
	flagsAsString="$flagsAsString -S  $FLAGS_output"
fi

echo  $flagsAsString;

# ensure that the parent folder is there
OUT_FOLDER=$(createOutputFolder "$FLAGS_output")

# run it
MESSAGE=$(hisat2 $flagsAsString)


RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile HISAT2)
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
		writeParam2File "$FLAGS_returnFilePath" "SAMFile" "$FLAGS_output"
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
