Troubleshooting unusual situations
==================================

There are a variety of problematic situations that might be
encountered when generating Fastqs from BCL data using the
``make_fastqs`` command; some of these are outlined below along
with suggestions on how to address them:

 * :ref:`troubleshooting-incomplete-run`
 * :ref:`troubleshooting-skip-demultiplexing`
 * :ref:`troubleshooting-handling-inline-barcodes`
 * :ref:`troubleshooting-no-adapter-trimming`
 * :ref:`troubleshooting-missing-fastqs-after-demultiplexing`

.. _troubleshooting-incomplete-run:

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

.. _troubleshooting-skip-demultiplexing:

Skip demultiplexing in ``make_fastqs`` stage
********************************************

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

.. _troubleshooting-handling-inline-barcodes:

Handling inline barcodes
************************

.. warning::

    Currently this is only implemented for **single-ended Fastqs**

In this situation the barcode index sequences are part of each read (e.g.
the first five bases of the first read), so ``bcl2fastq``'s standard
demultiplexing process can't be used.

In this case the following procedure can be used:

 * **Perform ``bcl`` to ``fastq`` conversion without demultiplexing**:
   put all the reads into a single fastq file by following the approach
   outlined in :ref:`troubleshooting-skip-demultiplexing` to avoid assigning
   index sequences to each read.

 * **Extract and assign inline barcodes**: use the ``assign_barcodes.py``
   utility to extract the barcode sequences from each read from the Fastq
   file produced by the previous step and assign these to the read header,
   for example::

       assign_barcodes.py -n 5 all_S1_R1_001.fastq.gz all_barcoded_S1_R1_001.fastq.gz

 * **Split into separate Fastq files by barcode sequence**: use the
   ``barcode_splitter.py`` utility to assign reads to individual Fastqs,
   for example::

       barcode_splitter.py -b ATACC -b TCTAG -b GCAGC all_barcoded_S1_R1_001.fastq.gz

.. _troubleshooting-no-adapter-trimming:

Tuning or turning off adapter trimming and masking
**************************************************

.. note::

   This only applies when using ``bcl2fastq`` version 2.

By default ``bcl2fastq`` version 2 performs adapter trimming and masking
on the reads in the output Fastq files, using the adapter sequences that
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

.. _troubleshooting-missing-fastqs-after-demultiplexing:

Missing Fastq files after demultiplexing by bcl2fastq
*****************************************************

If no reads match an index sequence in the sample sheet file, ``bcl2fastq``
will not produce a Fastq for that sample, leading to a verification
failure when the auto processor sees that some expected output Fastqs
are missing.

To workaround this use the ``--create-empty-fastqs`` option when
(re)running the ``make_fastqs`` command. This will create an empty
'placeholder' Fastq for each missing file, which enables verification to
complete successfully.

.. note::

   Before using this option it is recommended to check that the missing
   Fastqs are not due to some other problem or error in the data or
   pipeline.

.. warning::

   Be aware that the empty Fastqs may not be treated as valid input to
   some external downstream software packages.
