#!/bin/bash
#SKIP_TEST
SCRIPT_FOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SCRIPT_FOLDER/../../core_lib/includeBasics.sh # include some basic functions
SCRIPT=$(readlink_own -m $SCRIPT_FOLDER/contextMap.sh)
GLOBAL_TEST_DATA=$SCRIPT_FOLDER/../../test_data
TEST_DATA_FOLDER=$SCRIPT_FOLDER/test_data
FAILED_TESTS=0

# todo: configure theses values to enable the test and remove the #SKIP_TEST line
BWA_PATH="/path/to/bwa"
CM_PATH="/path/to/contextMap/CM.jar"

TMP_OUT="/tmp"
DEFAULT_PARAMS="--alignerName bwa --alignerBin '$BWA_PATH' -i $TEST_DATA_FOLDER/index/rDNA -g $TEST_DATA_FOLDER/ -j '$CM_PATH'"

TEST_S_NEW="$GLOBAL_TEST_DATA/fastq/joined/test_single_new.fastq"
TEST_S_OLD="$GLOBAL_TEST_DATA/fastq/joined/test_single_old.fastq"
TEST_P_NEW="$GLOBAL_TEST_DATA/fastq/joined/test_paired_new.fastq"
TEST_P_OLD="$GLOBAL_TEST_DATA/fastq/joined/test_paired_old.fastq"

# call with invalid parameter
testExitCode "/bin/bash $SCRIPT" "$EXIT_MISSING_ARGUMENTS" "Missing parameter test"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_DATA_FOLDER/notExistingFile.super -o $TMP_OUT/contextMap1.test" "$EXIT_MISSING_INPUT_FILES" "Missing input file test"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $GLOBAL_TEST_DATA/test.file.notReadable -o $TMP_OUT/contextMap2.test" "$EXIT_MISSING_INPUT_FILES" "Not readable input file test"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap3.test --skipsplit true,false" "$EXIT_INVALID_ARGUMENTS" "Wrong skipsplit I test"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap4.test --skipsplit notValid" "$EXIT_INVALID_ARGUMENTS" "Wrong skipsplit II test"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap5.test -j $TEST_DATA_FOLDER/notThere.jar" "$EXIT_MISSING_INPUT_FILES" "Missing jar file test"

# real calls
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap7.test" "$EXIT_OK" "Simple mapping test" "$TMP_OUT/contextMap7.test/test_single_new.sam"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap8.test --threads 4" "$EXIT_OK" "Threading test" "$TMP_OUT/contextMap8.test/test_single_new.sam"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap9.test --polyA" "$EXIT_OK" "PolyA detection test" "$TMP_OUT/contextMap9.test/test_single_new.polyA.bed"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS -r $TEST_S_NEW -o $TMP_OUT/contextMap10.test" "$EXIT_OK" "Single new format mapping test" "$TMP_OUT/contextMap10.test/test_single_new.sam"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS --pairedend -r $TEST_P_NEW -o $TMP_OUT/contextMap12.test" "$EXIT_OK" "Paired new format mapping test" "$TMP_OUT/contextMap12.test/test_paired_new.sam"
testExitCode "/bin/bash $SCRIPT $DEFAULT_PARAMS --pairedend -r $TEST_P_OLD -o $TMP_OUT/contextMap13.test" "$EXIT_OK" "Paired old format mapping test" "$TMP_OUT/contextMap13.test/test_paired_old.sam"

# delete all the temporary file
rm -rf $TMP_OUT/contextMap*.test 2>&1

# return the number of failed tests
exit $FAILED_TESTS
