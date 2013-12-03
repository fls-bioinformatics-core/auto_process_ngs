# Settings for auto_process module
#
import JobRunner
from bcf_utils import AttributeDictionary
#
# Define runners for specific jobs
runners = AttributeDictionary(
    bcl2fastq=JobRunner.SimpleJobRunner(),
    qc=JobRunner.SimpleJobRunner(),
    stats=JobRunner.SimpleJobRunner(),
)
#
# Information for archiving analyses
archive = AttributeDictionary(
    dirn=None,
    log=None,
)
#
# Information for uploading QC reports
qc_web_server = AttributeDictionary(
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
