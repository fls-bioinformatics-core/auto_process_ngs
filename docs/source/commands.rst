Autoprocessing commands
=======================

setup
*****

Creates a new top-level analysis directory for processing data from
a sequencing run.

make_fastqs
***********

Performs Fastq generation from bcl files.

setup_analysis_dirs
*******************

Creates subdirectories populated with Fastq files for each project.

run_qc
******

Run the QC scripts on the Fastqs files for each project.

publish_qc
**********

Copy reports from the QC runs to a webserver or other location.

archive
*******

Copy the final data to an 'archive' location.

