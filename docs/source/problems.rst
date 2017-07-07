Handling problem situations
===========================

There are a variety of relatively common problem situations that can be
encountered; these are outline below along with suggested protocols to
use to solve them.

 * :ref:`problem-missing-sample-sheet`
 * :ref:`problem-split-processing`
 * :ref:`problem-incomplete-run`
 * :ref:`problem-incorrect-barcodes`
 * :ref:`problem-skip-demultiplexing`
 * :ref:`problem-handling-inline-barcodes`
 * :ref:`problem-no-adapter-trimming`
 * :ref:`problem-missing-fastqs-after-demultiplexing`

.. _problem-missing-sample-sheet:

Missing Sample Sheet
********************

By default the original sample sheet file is assumed to be in the directory::

    <RUN_DIR>/Data/Intensities/BaseCalls/SampleSheet.csv

If this isn't found then the ``setup`` command will fail.

To address this, use the ``--sample-sheet`` option for ``setup`` to explicitly
specify the location and name of a non-default sample sheet. (If there is
no sample sheet at all then you will need to fabricate one first, then use
this option to import it on setup.)

.. _problem-split-processing:

Splitting processing for specific lanes with different parameters
*****************************************************************

This is the most common 'problem' situation, where a subset of lanes need
to be processed using different parameters compared to the others (for
example because one set of samples has different barcode lengths to the
others).

One procedure for handling this is to split the sample sheet into subsets
of lanes which share common processing parameters, run the Fastq
generation on each separately, and then "merge" the outputs before
proceeding with the QC etc.

For example: say we have a HISeq run where the majority of samples use
dual indexes, but samples in lanes 5 and 6 need to be demultiplexed using
only the first part of the dual index (i.e. first 8 bases only, for 16 base
dual indexes).

1. **Edit the sample sheet**

   First edit the sample sheet to make any necessary changes to
   to the lanes needing special treatment (e.g. updating barcodes).
   It might be preferrable to do this in a copy of the original sample
   sheet.

2. **Process the lanes which use default parameters first**

   Use the ``--lanes`` option of the ``make_fastqs`` command to process
   the first set of lanes, for example::

       auto_process.py make_fastqs --lanes=1-4,7-8

3. **Process the remaining 'problem' lanes**

   In this case we would need to use an updated bases mask to tell the
   demultiplexer to ignore the trailing 8 bases of the barcodes. We
   also want to output the Fastq files to a new directory to avoid clashing
   with the output from the first round of processing.

   For example::

       auto_process.py make_fastqs \
            --lanes=5-6 \
            --output-dir=bcl2fastq.lanes56 \
            --use-bases-mask=y101,I8,n8,y101 \
            --stats-file=statistics.lanes56.info \
            --skip-rsync

   Using ``--skip-rsync`` means that the processing doesn't try to fetch
   the raw data again.

4. **Merge fastq output directories and regenerate statistics**

   Assuming that the processing has completed okay there will be an
   initial ``bcl2fastq`` directory with the first set of Fastqs and
   another called ``bcl2fastq.lanes56``. The contents of this second
   directory can be merged into the first using the ``merge_fastq_dirs``
   command::

       auto_process.py merge_fastq_dirs

   To regenerate the statistics use the ``update_fastq_stats``::

       auto_process.py params --set stats_file=statistics.info
       auto_process.py update_fastq_stats

The remaining processing should then be performed as normal.

.. note::

   This process can be adapted to work with multiple subsets of
   lanes, not just two, by repeating step 3 above for each subset.

.. note::

   In older versions of the autoprocess software, it would be necessary
   to first split the sample sheet into two copies each containing with
   just one set of lanes (e.g. using
   ``prep_sample_sheet.py --include-lanes=...``), and making appropriate
   changes (e.g. to correct the barcodes).

.. _problem-incomplete-run:

Incomplete run/missing cycles
*****************************

If the sequencing run didn't complete then later cycles in the run won't be
present, and running the ``make_fastqs`` step will fail.

To address this:

1. **Fix the sample sheet:** if the run was truncated before the end of the
   index sequences then you will need to create a new sample sheet file with
   the index barcodes truncated to the appropriate length.

   This can be done using the ``prep_sample_sheet.py`` utility; for example if
   there are only 8bp of a 16bp index sequence then use::

       prep_sample_sheet.py --truncate-barcodes=8 \
            -o Samplesheet.8bp.csv \
            SampleSheet.csv

2. **Determine the corrected bases mask:** the ``bases_mask`` parameter in
   ``auto_process.info`` gives the default bases mask, which must be corrected
   to mask out the missing cycles.

   For example if the original bases mask was ``y101,I8,I8,y101`` but the run
   ended after the first index, then the updated bases mask would be
   ``y101,I8,n8,n101``.

3. **Generate the fastqs:** run ``make_fastqs`` specifying the updated sample
   sheet and bases mask, e.g.::

       auto_process.py make_fastqs \
            --sample-sheet=Samplesheet.8bp.csv \
            --use-bases-mask=y101,I8,n8,n101

.. _problem-incorrect-barcodes:

Incorrect barcode sequences in sample sheet
*******************************************

If one or more barcode sequences given in the original sample sheet were not
correct then demultiplexing will not be successful for the samples associated
with the 'bad' indices. Most commonly this manifests as an unusually small
number of reads for those samples, and a correspondingly larger than usual
number of undetermined reads.

To address this:

1. **Determine the actual barcode sequences:** use the ``analyse_barcodes``
   command for the lanes with the problem index sequences, e.g.::

        auto_process.py analyse_barcodes --lanes=6

   This will list the most common barcode sequences found, and should be
   sufficient to identify the true barcodes by eye, by comparing with the
   barcodes in the original sample sheet file.

2. **Reprocess the subset of lane(s):** use the procedure outlined in
   :ref:`problem-split-processing` to create a new sample sheet file for
   just the lane(s) with the bad indices, e.g.::

       prep_sample_sheet.py --include-lanes=6 -o SampleSheet.lane6.csv \
            custom_SampleSheet.csv

   Edit the barcodes in the new sample sheet file to replace the bad indices.
   NB don't remove any of the samples.

   Then rerun the Fastq generation using the new sample sheet file, and
   merge the outputs as described elsewhere.

.. _problem-skip-demultiplexing:

Skip demultiplexing in ``make_fastqs`` stage
********************************************

.. warning::

    This section is still under development

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

.. _problem-handling-inline-barcodes:

Handling inline barcodes
************************

.. warning::

    This section is still under development.

.. warning::

    Currently this is only implemented for **single-ended FASTQs**

In this situation the barcode index sequences are part of each read (e.g.
the first five bases of the first read), so ``bcl2fastq``'s standard
demultiplexing process can't be used.

In this case the following procedure can be used:

 * **Perform ``bcl`` to ``fastq`` conversion without demultiplexing**:
   put all the reads into a single fastq file by following the approach
   outlined in :ref:`problem-skip-demultiplexing` to avoid assigning
   index sequences to each read.

 * **Extract and assign inline barcodes**: use the ``assign_barcodes.py``
   utility to extract the barcode sequences from each read from the fastq
   file produced by the previous step and assign these to the read header,
   for example::

       assign_barcodes.py -n 5 all_S1_R1_001.fastq.gz all_barcoded_S1_R1_001.fastq.gz

 * **Split into separate fastq files by barcode sequence**: use the
   ``barcode_splitter.py`` utility to assign reads to individual fastqs,
   for example::

       barcode_splitter.py -b ATACC -b TCTAG -b GCAGC all_barcoded_S1_R1_001.fastq.gz

.. _problem-no-adapter-trimming:

Tuning or turning off adapter trimming and masking
**************************************************

.. note::

   This only applies when using ``bcl2fastq`` version 2.

By default ``bcl2fastq`` version 2 performs adapter trimming and masking
on the reads in the output FASTQ files, using the adapter sequences that
are provided in the input sample sheet file.

The default procedure it uses is:

 * Reads that contain sequence matching the adapters are trimmed to remove
   the matching sequence and all subsequent bases;

 * If a trimmed read is less than 35 bases long, it is padded with ``N``'s
   to make the length back up to 35 bases (this length can be modified
   using the ``--minimum-trimmed-read-length`` option of ``make_fastqs``);

 * If there are fewer than 22 non-``N`` bases in the read then the entire
   read is masked with ``N``'s (this length can be modified using the
   ``--mask-short-adapter-reads`` option of ``make_fastqs``).

There is no explicit switch to turn off the trimming and adapter masking,
however this can effectively be done by setting the adapter sequences in the
sample sheet to empty strings, for example::

    prep_sample_sheet.py -o SampleSheet.csv --set-adapter='' --set-adapter2='' SampleSheet.csv

.. _problem-missing-fastqs-after-demultiplexing:

Missing FASTQ files after demultiplexing by bcl2fastq
*****************************************************

If no reads match an index sequence in the sample sheet file, ``bcl2fastq``
will not produce a FASTQ for that sample, leading to a verification
failure when the auto processor sees that some expected output FASTQs
are missing.

To workaround this use the ``--create-empty-fastqs`` option when
(re)running the ``make_fastqs`` command. This will create an empty
'placeholder' FASTQ for each missing file, which enables verification to
complete successfully.

.. note::

   Before using this option it is recommended to check that the missing
   FASTQs are not due to some other problem or error in the data or
   pipeline.

.. warning::

   Be aware that the empty FASTQs may not be treated as valid input to
   some external downstream software packages.
