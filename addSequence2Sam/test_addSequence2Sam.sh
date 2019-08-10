#!/bin/bash
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh # include some basic functions
SCRIPT=$(readlink_own -m $SCRIPT_FOLDER/addSequence2Sam.pl)
GLOBAL_TEST_DATA=$SCRIPT_FOLDER/../../test_data
TEST_DATA_FOLDER=$SCRIPT_FOLDER/test_data
FAILED_TESTS=0

TMP_OUT="/tmp"

TEST_S_NEW="$GLOBAL_TEST_DATA/fastq/joined_sam/test_single_new.sam"
TEST_S_OLD="$GLOBAL_TEST_DATA/fastq/joined_sam/test_single_old.sam"
TEST_P_NEW="$GLOBAL_TEST_DATA/fastq/joined_sam/test_paired_new.sam"
TEST_P_OLD="$GLOBAL_TEST_DATA/fastq/joined_sam/test_paired_old.sam"
INVALID_SAM="$GLOBAL_TEST_DATA/fastq/joined_sam/invalid.sam"

TEST_S_F_NEW="$GLOBAL_TEST_DATA/fastq/joined/test_single_new.fastq"
TEST_S_F_OLD="$GLOBAL_TEST_DATA/fastq/joined/test_single_old.fastq"
TEST_P_F_NEW="$GLOBAL_TEST_DATA/fastq/joined/test_paired_new.fastq"
TEST_P_F_OLD="$GLOBAL_TEST_DATA/fastq/joined/test_paired_old.fastq"

# test if all needed modules are there
USED_TOOLS='^perl|Exporter|File::Basename|File::Path|Getopt::Long::Descriptive|Watchdog::ExitCode|strict|warnings'
MESSAGE=$($LIB_SCRIPT_FOLDER/checkUsedTools.sh "$USED_TOOLS")
CODE=$?

if [ $CODE -ne 0 ]; then
	echoError "$MESSAGE"
	exit $EXIT_TOOLS_MISSING
fi

# call with invalid parameter
testExitCode "perl $SCRIPT" "$EXIT_ARGUMENT_PROBLEM_PERL" "Missing parameter test"
testExitCode "perl $SCRIPT -s $TEST_DATA_FOLDER/notExistingFile.super -f $TEST_DATA_FOLDER/notExistingFile.super -o /tmp/addSequence2Sam1.test --returnfilepath /dev/null" "$EXIT_MISSING_INPUT_FILES" "Missing input file test"
testExitCode "perl $SCRIPT -s $GLOBAL_TEST_DATA/test.file.notReadable -f $GLOBAL_TEST_DATA/test.file.notReadable -o /tmp/addSequence2Sam2.test --returnfilepath /dev/null" "$EXIT_MISSING_INPUT_FILES" "Not readable input file test"
testExitCode "perl $SCRIPT -s $TEST_S_OLD -f $GLOBAL_TEST_DATA/fastq/invalid.fastq -o /tmp/addSequence2Sam3.test --returnfilepath /dev/null" "$EXIT_MISFORMATED_INPUT" "Not valid fastq test"  " " "/tmp/addSequence2Sam3.test"
testExitCode "perl $SCRIPT -s $INVALID_SAM -f $TEST_S_F_OLD -o /tmp/addSequence2Sam4.test --returnfilepath /dev/null" "$EXIT_MISFORMATED_INPUT" "Not valid SAM test"  " " "/tmp/addSequence2Sam4.test"

# real calls
testExitCode "perl $SCRIPT -s $TEST_S_OLD -f $TEST_S_F_OLD -o /tmp/addSequence2Sam5.test --returnfilepath /dev/null" "$EXIT_OK" "Simply add test" "/tmp/addSequence2Sam5.test"
testExitCode "perl $SCRIPT -s $TEST_S_OLD -f $TEST_DATA_FOLDER/missing.fastq -o /tmp/addSequence2Sam6.test --returnfilepath /dev/null" "$EXIT_MISFORMATED_INPUT" "Mising fastq entries test" " " "/tmp/addSequence2Sam6.test"
testExitCode "perl $SCRIPT -s $TEST_S_OLD -f $TEST_S_F_OLD -o /tmp/addSequence2Sam7.test --returnfilepath /dev/null" "$EXIT_OK" "Single end old format test" "/tmp/addSequence2Sam7.test"
testExitCode "perl $SCRIPT -s $TEST_S_NEW -f $TEST_S_F_NEW -o /tmp/addSequence2Sam8.test --returnfilepath /dev/null" "$EXIT_OK" "Single end new format test" "/tmp/addSequence2Sam8.test"
testExitCode "perl $SCRIPT -s $TEST_P_OLD -f $TEST_P_F_OLD -o /tmp/addSequence2Sam9.test --returnfilepath /dev/null" "$EXIT_OK" "Paired end old format test" "/tmp/addSequence2Sam9.test"
testExitCode "perl $SCRIPT -s $TEST_P_NEW -f $TEST_P_F_NEW -o /tmp/addSequence2Sam10.test --returnfilepath /dev/null" "$EXIT_OK" "Paired end new format test" "/tmp/addSequence2Sam10.test"

# delete all the temporary file
rm -f $TMP_OUT/addSequence2Sam*.test 2>&1 > /dev/null

# return the number of failed tests
exit $FAILED_TESTS
