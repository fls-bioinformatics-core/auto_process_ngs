#!/bin/sh -e
#
# process_miseq.sh
#
# Wrapper to manage data and run each stage of the
# auto_process_illumina.sh script for MiSEQ data
#
echo $(basename $0) version 2.0.0
# Get settings
settings_file=$(dirname $0)/process_miseq_setup.sh
if [ ! -f "$settings_file" ] ; then
   echo ERROR no settings file $settings_file
   exit 1
fi
echo Acquiring settings from $settings_file
. $settings_file
echo "ARCHIVE_DIR : $ARCHIVE_DIR"
echo "SEQ_DATA_LOG: $SEQ_DATA_LOG"
# Deal with command line
if [ $# -eq 0 ] ; then
   echo "Usage: $0 DATA_DIR DESCRIPTION"
   exit 1
fi
DATA_DIR=$1
if [ -z "$DATA_DIR" ] ; then
   echo ERROR no data dir supplied
   exit 1
fi
DESCRIPTION=$2
if [ -z "$DESCRIPTION" ] ; then
   echo ERROR no description supplied
   exit 1
fi
ANALYSIS_DIR=$(basename $DATA_DIR)_analysis
# Make local copy of primary data
rsync -av -m --include=*/ --include=Data/** --include=RunInfo.xml --include=SampleSheet.csv --exclude=* $DATA_DIR .
# Make analysis dir
LOCAL_DATA_DIR=$(pwd)/$(basename $DATA_DIR)
auto_process_illumina.sh setup miseq $LOCAL_DATA_DIR
# Move into analysis dir and do processing
cd $ANALYSIS_DIR
auto_process_illumina.sh make_fastqs
auto_process_illumina.sh run_qc
# Remove local copy of primary data
##/bin/rm -rf $LOCAL_DATA_DIR
# Copy results to archive area
cd ..
rsync_seq_data.py --mirror $ANALYSIS_DIR $ARCHIVE_DIR
# Log the run
DESCRIPTION="$(analyse_illumina_run.py --summary $ANALYSIS_DIR) $DESCRIPTION"
log_seq_data.sh $SEQ_DATA_LOG $ARCHIVE_DIR/$(date +%Y)/miseq/$ANALYSIS_DIR "$DESCRIPTION"
##
#
