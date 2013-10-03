# Settings for auto_process module
#
# Customise for local site
import JobRunner
from auto_process_utils import AttributeDict
#
# Define runners for specific jobs
runners = AttributeDict(
    default=JobRunner.SimpleJobRunner(),
    bcl2fastq=JobRunner.SimpleJobRunner(),
    qc=JobRunner.SimpleJobRunner(),
)
