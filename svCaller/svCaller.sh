#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
# TODO: add tools which are used in this script here in order to check, if they are installed on the system
USED_TOOLS='python3'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'help' '' "show this help message and exit" 'h1'
DEFINE_string 'bed' '' "Path to bedgraph file" 'bed_graph'
DEFINE_integer 'min_cld' '100' "[optional] The mininum distance of two clusters at which they still get combined" 'minimum_cluster_difference'
DEFINE_integer 'min_size' '2' "[optional] Minimum size of a deletion." 'minimum_size'
DEFINE_float 'max_z' '0.0' "[optional] Maximum z score threshold for coverage analysis" 'maximum_z_score'
DEFINE_float 'max_direct' '-2.5' "[optional] Maximum direct z score threshold for coverage analysis" 'maximum_direct_z_score'
DEFINE_float 'max_local' '-6.0' "[optional] Maximum local z score threshold for coverage analysis" 'maximum_local_z_score'
DEFINE_integer 'range' '500' "[optional] Size of range/region before a certain position, used for the determination of local z Score parameters" 'local_range'
DEFINE_integer 'pc' '1' "[optional] Pseudo count for coverages over positions" 'pseudo_count'
DEFINE_float 'tol' '0.8' "[optional] Tolerance of insertion positions mapped to deletions" 'tolerance'
DEFINE_string 'bam' '' "Path to bam file used for clipping patter analysis" 'bamfile'
DEFINE_string 'out_del' '' "Path to output txt file containing deletions" 'outpath_deletions'
DEFINE_string 'out_ins' '' "Path to output txt file containing insertions" 'outpath_insertions'
DEFINE_integer 'max_patt_diff' '10' "[optional] Maximum distance of peaks of clipped reads to count them as insertion" 'maximum_pattern_diff'
DEFINE_float 'min_sur_z' '50.0' "[optional] Minimum local z score for clipping pattern analysis" 'minimum_surounding_z_score'
DEFINE_integer 'ws' '20' "[optional] Size of the window, whose position are controlled to be significantly low" 'window_size'
DEFINE_float 'min_z' '10.0' "[optional] Minimum z score for clipping pattern analysis" 'min_z_score'
DEFINE_string 'get_clp_file' '' "[optional] Set this paramter as a path to get a file containing for each position the number of clipped reads" 'get_clipping_file'
DEFINE_integer 'min_reads' '10' "[optional] Minimum number of reads at which a position is permitted to be a peak" 'minimum_number_of_reads'
DEFINE_integer 'gen_prop' '1000' "[optional] Number of propagations to determine genome start/end" 'genome_propagate'
DEFINE_integer 'gap' '5' "[optional] Maximum number of permitted consecutive gaps/0-coverage positions during the determination of the genome start/end" 'gap_counter'
DEFINE_string 'ref' '' "Path to reference genome" 'reference'
DEFINE_float 'fir_ws' '0.0' "[optional] Primary threshold for the score, which is used for the verification of deletions with clipped sequences" 'first_weight_score'
DEFINE_float 'sec_ws' '1.0' "[optional] Secondary, more stringent threshold for the score, which is used for the verification of deletions with  clipped sequences" 'second_weight_score'
DEFINE_string 'con_path' '' "Path to the file containing the consensus sequences" 'consensus_path'
DEFINE_integer 'mpc' '1' "[optional] Small pseudo count for the log used for the computation of the PWMs" 'matrix_pseudo_count'
DEFINE_integer 'min_length' '10' "[optional] The minimum length of a consensus sequence" 'minimum_length'
DEFINE_integer 'clp_ver_range' '100' "[optional] The range of clipped positons of deletions, where consesus sequences are tried to match on" 'clipping_verfication_range'
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode

if [ "$FLAGS_version" -eq 0 ]; then
	MESSAGE="svCaller 1.0"
	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_bed" ]; then
	echoError "Parameter -bed_graph must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_bam" ]; then
	echoError "Parameter -bamfile must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_ref" ]; then
	echoError "Parameter -reference must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_con_path" ]; then
	echoError "Parameter -consensus_path must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_out_del" ]; then
	echoError "Parameter -outpath_deletions must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ -z "$FLAGS_out_ins" ]; then
	echoError "Parameter -outpath_insertions must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi
if [ "$FLAGS_min_cld" -lt 1 ] ; then
	echoError "Parameter -minimum_cluster_difference must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_min_size" -lt 2 ] ; then
	echoError "Parameter -minimum_size must be greater or equal to 2. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_range" -lt 1 ] ; then
	echoError "Parameter -local_range must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_pc" -lt 1 ] ; then
	echoError "Parameter -pseudo_count must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if (( $(echo "$FLAGS_tol < 0.0" | bc -l) )); then
    echoError "Parameter -tolerance must be greater or equal to 0. (see --help for details)"
    exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_max_patt_diff" -lt 1 ] ; then
	echoError "Parameter -maximum_pattern_diff must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_ws" -lt 1 ] ; then
	echoError "Parameter -window_size must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_min_reads" -lt 1 ] ; then
	echoError "Parameter -minimum_number_of_reads must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_gen_prop" -lt 1 ] ; then
	echoError "Parameter -genome_propagate must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_gap" -lt 0 ] ; then
	echoError "Parameter -gap_counter must be greater or equal to 0. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_mpc" -lt 1 ] ; then
	echoError "Parameter -matrix_pseudo_count must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_min_length" -lt 1 ] ; then
	echoError "Parameter -minimum_length must be greater or equal to 1. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi
if [ "$FLAGS_clp_ver_range" -lt 0 ] ; then
	echoError "Parameter -clipping_verfication_range must be greater or equal to 0. (see --help for details)";
	exit $EXIT_INVALID_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 

flagsAsString=""
if  [ ! -z "$FLAGS_help" ]; then
	flagsAsString="$flagsAsString -h1  $FLAGS_help"
fi
if  [ ! -z "$FLAGS_bed" ]; then
	flagsAsString="$flagsAsString -bed_graph  $FLAGS_bed"
fi
if  [ ! -z "$FLAGS_min_cld" ]; then
	flagsAsString="$flagsAsString -minimum_cluster_difference  $FLAGS_min_cld"
fi
if  [ ! -z "$FLAGS_min_size" ]; then
	flagsAsString="$flagsAsString -minimum_size  $FLAGS_min_size"
fi
if  [ ! -z "$FLAGS_max_z" ]; then
	flagsAsString="$flagsAsString -maximum_z_score  $FLAGS_max_z"
fi
if  [ ! -z "$FLAGS_max_direct" ]; then
	flagsAsString="$flagsAsString -maximum_direct_z_score  $FLAGS_max_direct"
fi
if  [ ! -z "$FLAGS_max_local" ]; then
	flagsAsString="$flagsAsString -maximum_local_z_score  $FLAGS_max_local"
fi
if  [ ! -z "$FLAGS_range" ]; then
	flagsAsString="$flagsAsString -local_range  $FLAGS_range"
fi
if  [ ! -z "$FLAGS_pc" ]; then
	flagsAsString="$flagsAsString -pseudo_count  $FLAGS_pc"
fi
if  [ ! -z "$FLAGS_tol" ]; then
	flagsAsString="$flagsAsString -tolerance  $FLAGS_tol"
fi
if  [ ! -z "$FLAGS_bam" ]; then
	flagsAsString="$flagsAsString -bamfile  $FLAGS_bam"
fi
if  [ ! -z "$FLAGS_out_del" ]; then
	flagsAsString="$flagsAsString -outpath_deletions  $FLAGS_out_del"
fi
if  [ ! -z "$FLAGS_out_ins" ]; then
	flagsAsString="$flagsAsString -outpath_insertions  $FLAGS_out_ins"
fi
if  [ ! -z "$FLAGS_max_patt_diff" ]; then
	flagsAsString="$flagsAsString -maximum_pattern_diff  $FLAGS_max_patt_diff"
fi
if  [ ! -z "$FLAGS_min_sur_z" ]; then
	flagsAsString="$flagsAsString -minimum_surounding_z_score  $FLAGS_min_sur_z"
fi
if  [ ! -z "$FLAGS_ws" ]; then
	flagsAsString="$flagsAsString -window_size  $FLAGS_ws"
fi
if  [ ! -z "$FLAGS_min_z" ]; then
	flagsAsString="$flagsAsString -min_z_score  $FLAGS_min_z"
fi
if  [ ! -z "$FLAGS_get_clp_file" ]; then
	flagsAsString="$flagsAsString -get_clipping_file  $FLAGS_get_clp_file"
	dirname_out_clp=$(dirname "$FLAGS_get_clp_file")
	if [ ! -d "$dirname_out_clp" ]; then
		mkdir -p "$dirname_out_clp"
	fi
fi
if  [ ! -z "$FLAGS_min_reads" ]; then
	flagsAsString="$flagsAsString -minimum_number_of_reads  $FLAGS_min_reads"
fi
if  [ ! -z "$FLAGS_gen_prop" ]; then
	flagsAsString="$flagsAsString -genome_propagate  $FLAGS_gen_prop"
fi
if  [ ! -z "$FLAGS_gap" ]; then
	flagsAsString="$flagsAsString -gap_counter  $FLAGS_gap"
fi
if  [ ! -z "$FLAGS_ref" ]; then
	flagsAsString="$flagsAsString -reference  $FLAGS_ref"
fi
if  [ ! -z "$FLAGS_fir_ws" ]; then
	flagsAsString="$flagsAsString -first_weight_score  $FLAGS_fir_ws"
fi
if  [ ! -z "$FLAGS_sec_ws" ]; then
	flagsAsString="$flagsAsString -second_weight_score  $FLAGS_sec_ws"
fi
if  [ ! -z "$FLAGS_con_path" ]; then
	flagsAsString="$flagsAsString -consensus_path  $FLAGS_con_path"
fi
if  [ ! -z "$FLAGS_mpc" ]; then
	flagsAsString="$flagsAsString -matrix_pseudo_count  $FLAGS_mpc"
fi
if  [ ! -z "$FLAGS_min_length" ]; then
	flagsAsString="$flagsAsString -minimum_length  $FLAGS_min_length"
fi
if  [ ! -z "$FLAGS_clp_ver_range" ]; then
	flagsAsString="$flagsAsString -clipping_verfication_range  $FLAGS_clp_ver_range"
fi
# run it

dirname_out_del=$(dirname "$FLAGS_out_del")
if [ ! -d "$dirname_out_del" ]; then
    mkdir -p "$dirname_out_del"
fi

dirname_out_ins=$(dirname "$FLAGS_out_ins")
if [ ! -d "$dirname_out_ins" ]; then
    mkdir -p "$dirname_out_ins"
fi

dirname_out_con=$(dirname "$FLAGS_con_path")
if [ ! -d "$dirname_out_con" ]; then
    mkdir -p "$dirname_out_con"
fi

MESSAGE=$(python3 $SCRIPT_FOLDER/svCaller.py $flagsAsString)
RET=$?

# check for error
FAIL=0
TMP_FILE=$(getTmpFile svCaller)
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


# exit with exit status 0 or exit code should be there earlier
echoError "Reached end of script! Exit should be performed earlier..."
exit $EXIT_REACHED_END_OF_SCRIPT
