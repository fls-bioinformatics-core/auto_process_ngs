Barcode Analysis
================

The ``analyse_barcodes.py`` utility can perform various analyses of the index
sequences (a.k.a. barcode sequences) in the read headers in one or more FASTQ
files.

Analysis can be done in two steps:

 1. Count the raw index sequences in the target FASTQs

 2. Analyse the raw counts to determine most numerous index sequences,
    group into similar sequences, and match groups to sequences in a
    sample sheet.

The first stage can be very time consuming so the counts can be output to an
intermediate ``counts`` file (using the ``-o FILE`` option), which can be
reloaded into the program (using the ``-c FILE`` option) to do the analyses
from the second stage more efficiently.

.. note::

   ``analyse_barcodes.py`` replaces the ``count_barcodes.py`` program from
   earlier versions.

Making the analyses quicker
---------------------------

For large numbers of reads the grouping analyses can be extremely time
consuming.

One way to make the grouping stage quicker is to used one or both of the
``--cutoff`` and ``--coverage`` options, each of which attempts to reduce
the number of reads that are used in the grouping analysis by throwing
away part of the "tail" of the barcode distribution:

 * ``--cutoff`` excludes barcodes which appear in fewer than a set
   percentage of reads, for example ``--cutoff 0.001`` excludes those
   which appear in less than 0.1% of reads.

 * ``--coverage`` limits the analysis to barcodes in the top percentage
   of reads, for example ``--coverage 0.9`` takes the barcodes from the
   top 90% of all reads.

