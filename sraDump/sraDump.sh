#!/bin/bash
# todo: add support for new illumina header if needed
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='fastq-dump:sed:tr:basename:awk'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'sraFile' '' "path to the .sra file(s); multiple SRA files muste be separated by ','" 's'
DEFINE_string 'sraID' '' "SRA ids; multiple ids muste be separated by ','" 'i'
DEFINE_string 'rename' '' "new basename for the resulting fastq files; multiple names muste be separated by ','" 'n'
DEFINE_string 'outputFolder' '' 'folder in which the files should be extracted' 'o'
DEFINE_string 'tmpFolder' '/usr/local/storage' 'tmp folder; default: /usr/local/storage' 't'
DEFINE_boolean 'deleteOnSuccess' '1' '[optional] deletes the SRA file when extraction was successfull' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode
BASE=""

# check if mandatory arguments are there
if [ -z "$FLAGS_sraFile" ] && [ -z "$FLAGS_sraID" ]; then
	echoError "Parameter -s or -i must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
else
	# ensure that not both parameters are used
	if [ ! -z "$FLAGS_sraFile" ] && [ ! -z "$FLAGS_sraID" ]; then
		echoError "Parameter -s and -i are exclusive. (see --help for details)";
		exit $EXIT_MISSING_ARGUMENTS
	else
		if [ ! -z "$FLAGS_sraFile" ]; then
			IFS=',' read -ra FILES <<< "$FLAGS_sraFile"
			for I in "${!FILES[@]}"; do
				FILE="${FILES[$I]}"
				FILE=$(echo "$FILE" | sed -E 's!^~!'$HOME'!')
				FILE=$(readlink_own -m $FILE)
				FILES[$I]=$FILE

				ENDING=$(echo "$FILE" | sed 's/^.*\.//')
				if [ "$ENDING" != "sra" ]; then
					echoError "No valid SRA file was found in '$FILE' for parameter --sraFile (see --help for details)";
					exit $EXIT_INVALID_ARGUMENTS
				fi
				# check, if the input file(s) exist
				verifyFileExistence "$FILE"

				# create the filename
				BASE=${BASE}_$(basename ${FILES[0]} ".sra")	
			done
		# work with IDs instead
		else
			IFS=',' read -ra FILES <<< "$FLAGS_sraID"
			for I in "${!FILES[@]}"; do
				FILE="${FILES[$I]}"
				FILES[$I]=$FILE
				BASE=${BASE}_${FILE}
			done
		fi

		# ensure that number of names is valid
		if [ ! -z "$FLAGS_rename" ]; then
			IFS=',' read -ra RENAME <<< "$FLAGS_rename"
			if [ "${#FILES[@]}" -ne "${#RENAME[@]}" ]; then
				echoError "Number of names differ from the number of samples. (see --help for details)";
				exit $INVALID_ARGUMENTS
			fi
		fi
	fi
fi
BASE=$(echo "$BASE" | sed -E 's/^_//') # remove trailing _

if [ -z "$FLAGS_outputFolder" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
# ensure that the parent folder is there
OUT_FOLDER=$(createOutputFolder "$FLAGS_outputFolder/xxx")
rm -f "$FLAGS_outputFolder/${BASE}.fastq" "$FLAGS_outputFolder/${BASE}_R1.fastq" "$FLAGS_outputFolder/${BASE}_R2.fastq" # ensure that these files are emtpy at the beginning
PAIRED_OLD="NO_DATA"

if [ ! -z "$FLAGS_tmpFolder" ]; then
	TMP_OUT_FOLDER=$(createOutputFolder "$FLAGS_tmpFolder/xxx")
else
	TMP_OUT_FOLDER=$OUT_FOLDER
fi

for I in "${!FILES[@]}"; do 
	FILE="${FILES[$I]}"
	echo "extracting $FILE..."

	executeCommand "fastq-dump -F --split-spot --skip-technical --split-3 -O \"$TMP_OUT_FOLDER\" \"$FILE\"" "/dev/null" "fastq-dump"
	B=$(basename $FILE ".sra")

	# check, the files that were written
	if [ -e "$TMP_OUT_FOLDER/${B}.fastq" ]; then
		PAIRED="false"
	else
		if [ -e "$TMP_OUT_FOLDER/${B}_1.fastq" ] && [ -e "$TMP_OUT_FOLDER/${B}_2.fastq" ] && [ ! -e "$TMP_OUT_FOLDER/${B}.fastq" ]; then
			PAIRED="true"
		else
			echoError "Was not able to detect if paired-end or single-end data was extracted."
			exit $EXIT_FAILED
		fi
	fi

	# ensure that all data is paired-end or all is single-end
	if [ "$PAIRED_OLD" != "NO_DATA" ] && [ "$PAIRED_OLD" != "$PAIRED" ]; then
		echoError "All archives must contain either single-end or paired-end data ($FILE vs. ${FILES[$I-1]})!";
		exit $EXIT_MISFORMATED_INPUT
	fi	

	echo "transforming $FILE..."
	NAME=${BASE}
	if [ ! -z "$FLAGS_rename" ]; then
		NAME="${RENAME[$I]}"
	fi
	# clean the data
	if [ "$PAIRED" == "false" ]; then
		mv "$TMP_OUT_FOLDER/${B}.fastq" "$TMP_OUT_FOLDER/${B}.fastq.tmp"
		awk '{if(NR % 4 == 3) {print "+"} else {print $0}}' "$TMP_OUT_FOLDER/${B}.fastq.tmp" >> "$FLAGS_outputFolder/${NAME}.fastq"
		rm "$TMP_OUT_FOLDER/${B}.fastq.tmp"
		CREATED_FILES="$FLAGS_outputFolder/${NAME}.fastq"
	else
		mv "$TMP_OUT_FOLDER/${B}_1.fastq" "$TMP_OUT_FOLDER/${B}_1.fastq.tmp"
		mv "$TMP_OUT_FOLDER/${B}_2.fastq" "$TMP_OUT_FOLDER/${B}_2.fastq.tmp"
		awk '{if(NR % 4 == 3) {print "+"} else {if(NR % 4 == 1) {print $0"/1"} else {print $0}}}' "$TMP_OUT_FOLDER/${B}_1.fastq.tmp" >> "$FLAGS_outputFolder/${NAME}_R1.fastq"
		awk '{if(NR % 4 == 3) {print "+"} else {if(NR % 4 == 1) {print $0"/2"} else {print $0}}}' "$TMP_OUT_FOLDER/${B}_2.fastq.tmp" >> "$FLAGS_outputFolder/${NAME}_R2.fastq"
		rm "$TMP_OUT_FOLDER/${B}_1.fastq.tmp" "$TMP_OUT_FOLDER/${B}_2.fastq.tmp"
		CREATED_FILES="$FLAGS_outputFolder/${NAME}_R1.fastq,$FLAGS_outputFolder/${NAME}_R2.fastq"
	fi

	# deletes the input file if parameter is set
	if [ $FLAGS_deleteOnSuccess -eq 0 ]; then
		rm "$FILE"
	fi
	PAIRED_OLD=$PAIRED
done

# check, if return parameters must be set
if [ ! -z "$FLAGS_returnFilePath" ]; then
	writeParam2File "$FLAGS_returnFilePath" "baseName" "$FLAGS_outputFolder/$NAME"
	writeParam2File "$FLAGS_returnFilePath" "isPairedEnd" "$PAIRED"
	writeParam2File "$FLAGS_returnFilePath" "createdFiles" "$CREATED_FILES"
	blockUntilFileIsWritten "$FLAGS_returnFilePath"
fi

# all seems to be ok
exit $EXIT_OK
