auto_process recipes
====================

Standard protocol
-----------------

Set up the analysis project from the sequencer output directory `<SOURCE>`
(containing the bcl files):

    auto_process.py setup <SOURCE>

Move into the resulting analysis directory and generate FASTQ files from
the bcls: 

    cd <SOURCE>_analysis
    auto_process.py make_fastqs

Edit the `projects.info` file then set up the project directories and run
the QC pipeline:

    auto_process.py setup_analysis_dirs
    auto_process.py run_qc

Copy the QC reports to a location where they can be served via the web:

    auto_process.py publish_qc

Copy the results to the 'archive' location:

    auto_process.py archive


Handling problem situations
---------------------------

### Missing SampleSheet.csv ###

By default the original sample sheet file is assumed to be in the directory

   <RUN_DIR>/Data/Intensities/BaseCalls/SampleSheet.csv

If this isn't found then the `setup` command will fail.

To address this, use the `--sample-sheet` option for `setup` to explicitly
specify the location and name of a non-default sample sheet. (If there is
no sample sheet at all then you will need to fabricate one first, then use
this option to import it on setup.)

### Incomplete run with missing cycles ###

If the sequencing run didn't complete then later cycles in the run won't be
present, and running the `make_fastqs` step will fail.

To address this:

1. Fix the sample sheet: if the run was truncated before the end of the index
   sequences then you will need to create a new sample sheet file with the
   index barcodes truncated to the appropriate length.

   This can be done using the `prep_sample_sheet.py` utility; for example if
   there are only 8bp of a 16bp index sequence then use

       prep_sample_sheet.py --truncate-barcodes=8 \
            -o Samplesheet.8bp.csv \
            SampleSheet.csv

2. Determine the corrected bases mask: the `bases_mask` parameter in
   `auto_process.info` gives the default bases mask, which must be corrected to
    mask out the missing cycles.

   For example if the original bases mask was `y101,I8,I8,y101` but the run
   ended after the first index, then the updated bases mask would be
   `y101,I8,n8,n101`.

3. Run `make_fastqs` specifying the updated sample sheet and bases mask, e.g.

       auto_process.py make_fastqs \
            --sample-sheet=Samplesheet.8bp.csv \
            --use-bases-mask=y101,I8,n8,n101

### Incorrect barcode sequences in sample sheet ###

If one or more barcode sequences given in the original sample sheet were not
correct then demultiplexing will not be successful for the samples associated
with the 'bad' indices. Most commonly this manifests as an unusually small
number of reads for those samples, and a correspondingly larger than usual
number of undetermined reads.

To address this:

1. Determine the true barcode sequences using the `analyse_barcodes` command
   for the lanes with the problem index sequences, e.g.

        auto_process.py analyse_barcodes --lanes=6

   This will list the most common barcode sequences found, and should be
   sufficient to identify the true barcodes by eye, by comparing with the
   barcodes in the original sample sheet file.

2. Create a new sample sheet file containing information on just the lane(s)
   with the bad indices, e.g.

       prep_sample_sheet.py --include-lanes=6 -o SampleSheet.lane6.csv \
            custom_SampleSheet.csv

   Edit the barcodes in the new sample sheet file to replace the bad indices.
   NB don't remove any of the samples.

3. Rerun the Fastq generation using the new sample sheet file, and sending the
   output to a new bcl2fastq directory e.g.

       auto_process.py make_fastqs \
            --sample-sheet=Samplesheet.lane6.csv \
            --output-dir=bcl2fastq.lane6 \
            --stats-file=statistics.lane6.info \
            --skip-rsync

   This will put the Fastq files for the reprocessed lanes into a new
   subdirectory (in this example `bcl2fastq.lane6`) and write the statistics
   to a new file (`statistics.lane6.info`).
  
   `--skip-rsync` uses the existing primary data files to speed up the process.

   Once this has completed the new statistics can be inspected to check that the
   demultiplexing has improved.

4. Merge the new bcl2fastq directory into the old one using the `merge_fastq_dirs`
   command, e.g.

       auto_process.py merge_fastq_dirs

   then update the statistics file, e.g.

       auto_process.py config --set stats_file=statistics.info
       auto_process.py update_fastq_stats

After this has been completed the analysis directory setup and QC steps can be
run as before.
