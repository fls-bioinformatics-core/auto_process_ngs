===================
Updating CellRanger
===================

This is an overview of updates that need to be made to the code to
add support for new versions of 10x Genomics ``cellranger`` family
of pipelines, and a suggested procedure for doing this.

It is also recommended to scan the release notes for the new
version to see if any significant changes to existing options or
outputs are flagged up.

**Step 1: Check and update version detection**

Run ``cellranger --version`` to see what the version string output
looks like, for example for CellRanger 7.1.0:

::

   $ cellranger --version
   cellranger cellranger-7.1.0

If the format differs from the previous version then the
``cellranger_info`` function in ``auto_process_ngs.tenx.utils``)
will need to be updated; add a new to ``TestCellrangerInfo`` and
update ``cellranger_info`` so that the test passes.

The version returned by the ``Mock10xPackageExe`` class in the
``mock`` module should also be updated in this case; the default
version should also be updated.

**Step 2: Run Fastq generation**

Next, update the configuration to use the new version for Fastq
generation, then run the pipeline on a test dataset and look
for errors.

Any failures will flag up issues with changes to the command line
for the ``mkfastq`` command; once the differences are identified
then the appropriate parts of the ``Mock10xPackageExe`` class in
the ``mock`` module should be updated, to mimick the behaviour for
new, modified or removed command line options (the CellRanger
documentation can help here).

A new version-specific unit test case for the ``MakeFastqs``
pipeline should also be added (e.g.
``test_makefastqs_10x_chromium_sc_protocol_710``).

**Step 3: Run single library analysis (count)**

Once any issues with the Fastq generation have been addressed,
run the QC pipeline to perform the single library analyses.
Again, any failures will flag up issues with changes to both
the command line and the outputs from ``count``: depending on
what has changed, these may require updates to:

* The ``Mock10xPackageExe`` implementation of ``count``
* CellRanger-specific elements of ``qc.pipeline``
* Code for identification and handling of output files and
  metrics (in ``tenx.metrics`` and ``qc.cellranger``).

Even if the pipeline runs without issues, it is recommended to
check the contents of the metric files in case the format has
changed; if so then example versions of the new contents can
be added to ``mock10xdata`` and new test cases created to
check that the appropriate classes handle it correctly.

(The ``Mock10xPackageExe`` should also be updated to use any
new example output data for its mock outputs.)

**Step 4: Run multiplexing analysis (multi)**

This essentially repeats step 3, with the ``multi`` command.
