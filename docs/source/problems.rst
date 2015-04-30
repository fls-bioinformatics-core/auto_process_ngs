Handling problem situations
===========================

There are a variety of relatively common problem situations that can be
encountered; these are outline below along with suggested protocols to
use to solve them.

Missing Sample Sheet
********************

By default the original sample sheet file is assumed to be in the directory::

    <RUN_DIR>/Data/Intensities/BaseCalls/SampleSheet.csv

If this isn't found then the `setup` command will fail.

To address this, use the `--sample-sheet` option for `setup` to explicitly
specify the location and name of a non-default sample sheet. (If there is
no sample sheet at all then you will need to fabricate one first, then use
this option to import it on setup.)

Incomplete run/missing cycles
*****************************

If the sequencing run didn't complete then later cycles in the run won't be
present, and running the `make_fastqs` step will fail.

To address this:

1. Fix the sample sheet: if the run was truncated before the end of the index
   sequences then you will need to create a new sample sheet file with the
   index barcodes truncated to the appropriate length.

   This can be done using the `prep_sample_sheet.py` utility; for example if
   there are only 8bp of a 16bp index sequence then use::

       prep_sample_sheet.py --truncate-barcodes=8 \
            -o Samplesheet.8bp.csv \
            SampleSheet.csv

2. Determine the corrected bases mask: the `bases_mask` parameter in
   `auto_process.info` gives the default bases mask, which must be corrected to
   mask out the missing cycles.

   For example if the original bases mask was `y101,I8,I8,y101` but the run
   ended after the first index, then the updated bases mask would be
   `y101,I8,n8,n101`.

3. Run `make_fastqs` specifying the updated sample sheet and bases mask, e.g.::

       auto_process.py make_fastqs \
            --sample-sheet=Samplesheet.8bp.csv \
            --use-bases-mask=y101,I8,n8,n101

Incorrect barcode sequences in sample sheet
*******************************************

If one or more barcode sequences given in the original sample sheet were not
correct then demultiplexing will not be successful for the samples associated
with the 'bad' indices. Most commonly this manifests as an unusually small
number of reads for those samples, and a correspondingly larger than usual
number of undetermined reads.

To address this:

1. Determine the true barcode sequences using the `analyse_barcodes` command
   for the lanes with the problem index sequences, e.g.::

        auto_process.py analyse_barcodes --lanes=6

   This will list the most common barcode sequences found, and should be
   sufficient to identify the true barcodes by eye, by comparing with the
   barcodes in the original sample sheet file.

2. Create a new sample sheet file containing information on just the lane(s)
   with the bad indices, e.g.::

       prep_sample_sheet.py --include-lanes=6 -o SampleSheet.lane6.csv \
            custom_SampleSheet.csv

   Edit the barcodes in the new sample sheet file to replace the bad indices.
   NB don't remove any of the samples.

3. Rerun the Fastq generation using the new sample sheet file, and sending the
   output to a new bcl2fastq directory e.g.::

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

4. Merge the new bcl2fastq directory into the old one using the
   `merge_fastq_dirs` command, e.g.::

       auto_process.py merge_fastq_dirs

   then update the statistics file, e.g.::

       auto_process.py config --set stats_file=statistics.info
       auto_process.py update_fastq_stats

After this has been completed the analysis directory setup and QC steps can be
run as before; note that you will need to reset where `auto_process.py` looks
for the Fastq files before running `setup_project_dirs`, i.e.::

    auto_process.py config --set unaligned_dir=bcl2fastq


Skip demultiplexing in ``make_fastqs`` stage
********************************************

.. warning::

    This section is still under development and might be inaccurate or even
    completely wrong!

The demultiplexing can be skipped in one of two ways.

To process each lane without any demultiplexing, edit the sample sheet so
that there is only one "sample" defined for each lane, and remove any barcode
index sequence.

For example::

    FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
    FC1,1,Lane1,,,,,,,AllReads

Then update the bases mask so that the index sequences are either ignored or
are collected as part of the reads.

For example, if the initial bases mask was ``y300,I8,I8,y300`` then set this to
``y300,n8,n8,y300`` to ignore them (in which case index sequences will be lost)
or to e.g. ``y316,y300`` (in which case the last 16 bases of each R1 read will
be the index sequence).

Note that in either case, the index sequence will not appear in the header for
each read.

Alternatively a pseudo-demultiplexing approach can be used, by specifying a single
"sample" in the sample sheet but this time including an appropriate length index
sequence which cannot be matched::

    FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
    FC1,1,Lane1,,AAAAAAAA-AAAAAAAA,,,,,AllReads

Using this approach should put all the reads into the "undetermined" project;
however this way the index sequences should still have been captured in the read
headers.

