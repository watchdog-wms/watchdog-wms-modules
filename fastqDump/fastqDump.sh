#!/bin/bash
# todo: add support for new illumina header if needed
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh

# check, if used tools are installed
USED_TOOLS='prefetch:fastq-dump:sed:tr:basename:awk'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS" "check_shflag_tools")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoNew "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# define parameters
DEFINE_string 'sraId' '' 'SRA id' 's'
DEFINE_string 'outputFolder' '' 'folder in which the files should be extracted' 'o'
DEFINE_string 'pathToAspera' '' '[optional] use aspera to speedup download' 'a'
DEFINE_boolean 'checkPresent' '1' '[optional] check if files already present and download successfull' ''
DEFINE_string 'returnFilePath' '' 'path to the return variables file' ''
DEFINE_boolean 'version' 'false' '[optional] prints the version' 'v'
DEFINE_boolean 'debug' 'false' '[optional] prints out debug messages.' ''

# parse parameters
FLAGS "$@" || exit $EXIT_INVALID_ARGUMENTS
eval set -- "${FLAGS_ARGV}"
printParamValues "initial parameters" # print param values, if in debug mode
SRA_ID=""

if [ "$FLAGS_version" -eq 0 ]; then
	
	MESSAGE=$(prefetch --version | head -n 2 | tail -n 1 | awk '{print $3}');

	echo $MESSAGE
	exit $EXIT_OK
fi

# check if mandatory arguments are there
if [ -z "$FLAGS_sraId" ]; then
	echoError "Parameter -s must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
else
	SRA_ID="$FLAGS_sraId"
fi

if [ -z "$FLAGS_outputFolder" ]; then
	echoError "Parameter -o must be set. (see --help for details)";
	exit $EXIT_MISSING_ARGUMENTS
fi

printParamValues "parameters before actual script starts" # print param values, if in debug mode
##################################################### START with actual SCRIPT ##################################################### 
# ensure that the parent folder is there
OUT_FOLDER=$(createOutputFolder "$FLAGS_outputFolder/xxx")

RUN=0;


if [ $FLAGS_checkPresent -eq 0 ] && [ -e "$FLAGS_outputFolder/${SRA_ID}.log" ]; then
	file_time_log=$(date -r "$FLAGS_outputFolder/${SRA_ID}.log" +%s);
	file_time_out=0;
	
	if [ -e "$FLAGS_outputFolder/${SRA_ID}_1.fastq" ] && [ -e "$FLAGS_outputFolder/${SRA_ID}_2.fastq" ]; then
		PAIRED="true";
		file_time_out=$(date -r "$FLAGS_outputFolder/${SRA_ID}_2.fastq" +%s);
	elif [ ! -e "$FLAGS_outputFolder/${SRA_ID}.fastq" ]; then
		PAIRED="false";
		file_time_out=$(date -r "$FLAGS_outputFolder/${SRA_ID}.fastq" +%s);
	fi
	
	if [ $file_time_out -ne 0 ] && [ $file_time_out -le $file_time_log ]; then
		COUNT=$(grep -c "Written .* spots for $SRA_ID" "$FLAGS_outputFolder/$SRA_ID.log");
		
		if [ $COUNT -eq 1 ]; then
			RUN=1;	
			echo "Download previously successfull. Will not download again."
		fi
		
	fi
	

fi

if [ $RUN -ne 1 ]; then
	rm -f "$FLAGS_outputFolder/${SRA_ID}.fastq" "$FLAGS_outputFolder/${SRA_ID}_1.fastq" "$FLAGS_outputFolder/${SRA_ID}_2.fastq" # ensure that these files are emtpy at the beginning

	echo "extracting $SRA_ID..."
	
	echo "use prefetch..."
	
	COMMAND_PRE="prefetch";
	
	
	if [ -n "$FLAGS_pathToAspera" ]; then
		echo "use aspera..."
		COMMAND_PRE="$COMMAND_PRE -t ascp -a \"$FLAGS_pathToAspera/bin/ascp|$FLAGS_pathToAspera/etc/asperaweb_id_dsa.openssh\"";
	fi
	
	COMMAND_PRE="$COMMAND_PRE $SRA_ID";
	
	echo $COMMAND_PRE;
	
	MESSAGE=$(eval "$COMMAND_PRE" 2>&1 | tee "$FLAGS_outputFolder/$SRA_ID.log")
	CODE=$?
	
	if [ $CODE -ne 0 ]; then
		echoError "prefetch of $SRA_ID not successfull."
		exit $EXIT_FAILED;									
	fi

	COMMAND="fastq-dump -F --split-3 -O \"$FLAGS_outputFolder\" \"$SRA_ID\"";

	FAIL=0
	MESSAGE=$(eval "$COMMAND" 2>&1 | tee -a "$FLAGS_outputFolder/$SRA_ID.log")
	CODE=$?



	COUNT=$(grep -c "Written .* spots for $SRA_ID" "$FLAGS_outputFolder/$SRA_ID.log")
	
	echo "$COUNT reads downloaded";
		
	if [ $COUNT -ne 1 ] && [ $CODE -ne 0 ]; then
		echoError "fastq-dump of $SRA_ID not successfull."
		exit $EXIT_FAILED;									
	fi

	# check, the files that were written
	if [ -e "$FLAGS_outputFolder/${SRA_ID}.fastq" ]; then
		PAIRED="false"
		READ1FILE="$FLAGS_outputFolder/${SRA_ID}.fastq";
		READ2FILE="$READ1FILE";
	else
		if [ -e "$FLAGS_outputFolder/${SRA_ID}_1.fastq" ] && [ -e "$FLAGS_outputFolder/${SRA_ID}_2.fastq" ] && [ ! -e "$FLAGS_outputFolder/${SRA_ID}.fastq" ]; then
			PAIRED="true"
			READ1FILE="$FLAGS_outputFolder/${SRA_ID}_1.fastq";
			READ2FILE="$FLAGS_outputFolder/${SRA_ID}_2.fastq";
		else
			echoError "Was not able to detect if paired-end or single-end data was extracted."
			exit $EXIT_FAILED
		fi
	fi
fi



# check, if return parameters must be set
if [ ! -z "$FLAGS_returnFilePath" ]; then
	writeParam2File "$FLAGS_returnFilePath" "isPairedEnd" "$PAIRED"
	writeParam2File "$FLAGS_returnFilePath" "readFile1" "$READ1FILE"
	writeParam2File "$FLAGS_returnFilePath" "readFile2" "$READ2FILE"
	blockUntilFileIsWritten "$FLAGS_returnFilePath"
fi

# all seems to be ok
exit $EXIT_OK
