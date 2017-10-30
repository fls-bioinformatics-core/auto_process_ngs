Standard protocol
=================

Set up the analysis project from the sequencer output directory
``<SOURCE>`` (containing the bcl files):

::

    auto_process.py setup <SOURCE>

Move into the resulting analysis directory and generate FASTQ files from
the bcls:

::

    cd <SOURCE>_analysis
    auto_process.py make_fastqs

Edit the ``projects.info`` file then set up the project directories and
run the QC pipeline:

::

    auto_process.py setup_analysis_dirs
    auto_process.py run_qc

Copy the QC reports to a location where they can be served via the web:

::

    auto_process.py publish_qc

Stage the results to a "pending" directory in the 'archive' location:

::

    auto_process.py archive

Copy the results to the final archive location:

::

    auto_process.py archive --final
