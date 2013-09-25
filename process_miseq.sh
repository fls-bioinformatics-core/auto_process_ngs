#!/bin/sh -e
#
# process_miseq.sh
#
# Wrapper to manage data and run each stage of the
# auto_process_illumina.sh script for MiSEQ data
#
echo $(basename $0) version 1.0.0
#
DATA_DIR=$1
if [ -z "$DATA_DIR" ] ; then
   echo ERROR no data dir supplied
   exit 1
fi
ANALYSIS_DIR=$(basename $DATA_DIR)_analysis
# Make local copy of primary data
rsync -av -m --include=*/ --include=Data/** --include=RunInfo.xml --include=SampleSheet.csv --exclude=* $DATA_DIR .
# Make analysis dir
DATA_DIR=$(pwd)/$(basename $DATA_DIR)
auto_process_illumina.sh setup miseq $DATA_DIR
# Move into analysis dir and do processing
cd $ANALYSIS_DIR
auto_process_illumina.sh make_fastqs
auto_process_illumina.sh run_qc
# Move out of analysis dir and rsync
cd ..
rsync_seq_data.py --mirror $ANALYSIS_DIR /mnt/kadmon/bcf/ngsdata/Analysis
##
#
