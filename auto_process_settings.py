# Settings for auto_process module
#
import JobRunner
from auto_process_utils import AttributeDict
#
# Define runners for specific jobs
runners = AttributeDict(
    bcl2fastq=JobRunner.SimpleJobRunner(),
    qc=JobRunner.SimpleJobRunner(),
)
#
# Information for archiving analyses
archive = AttributeDict(
    dirn=None,
    log=None,
)
#
# Information for uploading QC reports
qc_web_server = AttributeDict(
    user=None,
    server=None,
    webdir=None,
    url=None,
)
#
# Import site-specific settings from local version
try:
    from auto_process_settings_local import *
except ImportError,ex:
    print "Unable to import local settings: %s" % ex
