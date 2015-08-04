#!/bin/bash -e
#
# process_miseq.sh
#
# Wrapper to manage data and run each stage of the
# auto_process_illumina.sh script for MiSEQ data
#
echo $(basename $0) version 2.1.0
# Get settings
settings_file=$(dirname $0)/process_miseq_setup.sh
if [ ! -f "$settings_file" ] ; then
   echo ERROR no settings file $settings_file
   exit 1
fi
echo Acquiring settings from $settings_file
. $settings_file
echo "ARCHIVE_DIR       : $ARCHIVE_DIR"
echo "SEQ_DATA_LOG      : $SEQ_DATA_LOG"
echo "QC_WEB_SERVER     : $QC_WEB_SERVER"
echo "QC_WEB_SERVER_USER: $QC_WEB_SERVER_USER"
echo "QC_WEB_SERVER_DIR : $QC_WEB_SERVER_DIR"
echo "QC_WEB_SERVER_URL : $QC_WEB_SERVER_URL"
if [ -z "$ARCHIVE_DIR" ] || [ -z "$SEQ_DATA_LOG" ] || \
    [ -z "$QC_WEB_SERVER" ] || [ -z "$QC_WEB_SERVER_USER" ] || \
    [ -z "$QC_WEB_SERVER_DIR" ] || [ -z "$QC_WEB_SERVER_URL" ] ; then
   echo ERROR some settings are blank, stopping
   exit 1
fi
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
# Copy results to archive area
cd ..
rsync_seq_data.py --mirror $ANALYSIS_DIR $ARCHIVE_DIR
# Log the run
DESCRIPTION="$DESCRIPTION $(analyse_illumina_run.py --summary $ANALYSIS_DIR)"
log_seq_data.sh $SEQ_DATA_LOG $ARCHIVE_DIR/$(date +%Y)/miseq/$ANALYSIS_DIR "$DESCRIPTION"
# Transfer qc reports to web server
qc_web_server=${QC_WEB_SERVER_USER}@${QC_WEB_SERVER}
echo Copying QC reports to web server: $qc_web_server
ssh $qc_web_server mkdir ${QC_WEB_SERVER_DIR}/$(basename $DATA_DIR)
find $ANALYSIS_DIR -name "qc_report.*.zip" -exec scp {} $qc_web_server:${QC_WEB_SERVER_DIR}/$(basename $DATA_DIR) \;
# Summarise
echo QC reports available at ${QC_WEB_SERVER_URL}/$(basename $DATA_DIR)
grep $ANALYSIS_DIR $SEQ_DATA_LOG
##
#
