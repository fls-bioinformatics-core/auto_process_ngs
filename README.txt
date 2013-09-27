auto_process
============

Scripts and utilities for automatic processing & management of NGS sequencing
data.

process_miseq.sh: simple wrapper script for auto_process_illumina.sh, to
   automatically process MiSEQ data.
   Requires settings to be added to local copy of process_miseq_setup.sh
   (use process_miseq_setup.sh.sample as starting point).

applications.py: Python module with utilities for generating and running
   NGS-related command line programs.

bclToFastq.py: run the CASAVA bcl to fastq conversion pipeline.
