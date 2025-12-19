===================
Updating CellRanger
===================

This is an overview of updates that need to be made to the code to
add support for new versions of 10x Genomics ``cellranger`` family
of pipelines, and a suggested procedure for doing this.

It is also recommended to scan the release notes for the new
version to see if any significant changes to existing options or
outputs are flagged up.

Background: what needs to be updated
------------------------------------

* ``tenx`` module: the default CellRanger version should be set to the
  new version;
* ``mock10xdata`` module: add constants with example ``metrics_summary.csv``
  outputs from CellRanger's ``count`` command, and from the ``multi``
  command (for both CellPlex and Flex);
* ``mock`` module: the ``Mock10xPackage`` should be extended to mimick
  the new version:

  - returning the correct version string
  - add support for new options for the ``count`` and ``multi`` (and
    removing any that are no longer available)
  - generate the appropriate outputs (using the examples added to the
    ``mock10xdata`` module)

* ``mockqc`` module: the ``cellranger_multi`` method of the ``MockQCOutputs``
  class should be updated to recognise the new version and create
  appropriate outputs (also using the examples added to ``mock10xdata``)
* ``tenx.utils`` module: add new tests for the ``cellranger_info``
  function to handle output from ``cellranger --version`` for the new
  version (and update the function if necessary)
* ``tenx.metrics`` module: implement tests to check that the
  ``GexSummary`` and ``MultiplexSummary`` classes can handle the metrics
  output from the new version (and update the classes if necessary)
* ``qc.modules.cellranger_count`` module: implement tests for the
  ``cellranger_count`` QC module to handle the new CellRanger version,
  and update if necessary
* ``qc.modules.cellranger_multi`` module: perform similar updates for the
  ``cellranger_multi`` QC module.

Steps to perform the updates
----------------------------

In practice it isn't always convenient to implement the changes outlined
above in the specific order they're listed. The following steps give a
more pragmatic guide to gathering the information and making the changes.

**Step 1: Check and update version detection**

Run ``cellranger --version`` to see what the version string output
looks like, for example for CellRanger 10.0.0:

::

   $ cellranger --version
   cellranger 10.0.0

If the format differs from the previous version then the
``cellranger_info`` function in ``tenx.utils``) will need to be updated; add
a new test to ``TestCellrangerInfo`` and update ``cellranger_info`` so that
the test passes.

Although not actually required until later, it's a good idea to also
update default version of CellRanger which is used for mocking and unit tests:
this is done by changing the ``DEFAULT_CELLRANGER_VERSION`` in
``tenx.__init__`` to match the new version, and updating the handling of the
``--version`` argument in the ``Mock10xPackageExe`` class in the ``mock``
module (so that the correctly formatted version string is returned).

**Step 2: Run Fastq generation**

Update the configuration to use the new version for Fastq generation
and run the ``make_fastq`` command on a test dataset to look for errors.

Any failures will flag up issues with changes to the command line
for the ``mkfastq`` command; once the differences are identified
then the appropriate parts of the ``Mock10xPackageExe`` class in
the ``mock`` module should be updated, to mimick the behaviour for
new, modified or removed command line options (the CellRanger
documentation can help here).

A new version-specific unit test case for the ``MakeFastqs``
pipeline should also be added (e.g.
``test_makefastqs_10x_chromium_sc_protocol_10_0_0``).

.. note::

   Since CellRanger 9.0.0 the use of ``mkfastq`` is deprecated in
   favour of running ``bcl2fastq`` or ``bcl-convert`` directly,
   although it appears to be supported still in version 10.0.0.
   However these checks may become obsolete in future.

**Step 3: Run single library analysis (count)**

Once any issues with the Fastq generation have been addressed,
run the QC pipeline to perform the single library analyses.
Again, any failures will flag up issues with changes to both
the command line and the outputs from ``count``: depending on
what has changed, these may require updates to:

* The ``Mock10xPackageExe`` implementation of ``count``
* CellRanger-specific elements of ``qc.pipeline``
* Code for identification and handling of output files and
  metrics in ``tenx.metrics``

.. note::

   If there are significant changes to the structure of the
   CellRanger ``count`` outputs then updates may also be
   needed to code in the ``qc.apps.cellranger`` module.

Even if the pipeline runs without issues, it is recommended to
check the contents of the metric files in case the format has
changed; if so then example versions of the new contents can
be added to ``mock10xdata`` and new test cases created to
check that the appropriate classes handle it correctly.

(The ``Mock10xPackageExe`` should also be updated to use any
new example output data for its mock outputs.)

**Step 4: Run multiplexing analysis (multi)**

This essentially repeats step 3, with the ``multi`` command.

Also the ``cellranger_multi`` method of the ``MockQCOutputs``
class in the ``mockqc`` module should be updated to generate
appropriate outputs for the new CellRanger version.

**Step 5: Run all unit tests**

This is a final sanity check to ensure that the updates haven't
broken existing functionality, and there are no unexpected bugs.

.. note::

   It's likely that there will some failures if the sample
   ``count`` and ``multi`` outputs have different cell counts etc
   from the new version compared with the previous one - these
   will need to be updated.