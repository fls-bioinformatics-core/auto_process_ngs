# Settings for auto_process module
#
import os
import sys
import JobRunner
from bcf_utils import AttributeDictionary
#
# General parameters
general = AttributeDictionary(
    max_concurrent_jobs=12,
)
#
# bcl2fastq
bcl2fastq = AttributeDictionary(
    nprocessors=1,
)
#
# Define runners for specific jobs
runners = AttributeDictionary(
    bcl2fastq=JobRunner.SimpleJobRunner(join_logs=True),
    qc=JobRunner.SimpleJobRunner(join_logs=True),
    stats=JobRunner.SimpleJobRunner(join_logs=True),
)
#
# Information for archiving analyses
# dirn should be a directory in the form [[user@]host:]path]
archive = AttributeDictionary(
    dirn=None,
    log=None,
    group=None,
    chmod=None
)
#
# Information for uploading QC reports
# dirn should be a directory in the form [[user@]host:]path]
qc_web_server = AttributeDictionary(
    dirn=None,
    url=None,
)
#
# Import site-specific settings from local version
__install_path = os.path.abspath(os.path.normpath(
    os.path.join(os.path.dirname(sys.argv[0]),'..')))
if os.path.exists(os.path.join(__install_path,'auto_process_settings_local.py')):
    try:
        sys.path.append(__install_path)
        from auto_process_settings_local import *
    except ImportError,ex:
        print "Unable to import local settings from %s: %s" % (__install_path,ex)
else:
    print "No local settings file in %s" % __install_path
