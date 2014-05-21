# Settings for auto_process module
#
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
try:
    from auto_process_settings_local import *
except ImportError,ex:
    print "Unable to import local settings: %s" % ex
