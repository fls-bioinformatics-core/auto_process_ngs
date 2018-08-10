Moving data to the archive location using ``auto_process archive``
==================================================================

Copy the final data to an 'archive' location::

   auto_process.py archive [ANALYSIS_DIR]

Stage the results to a "pending" directory in the 'archive' location:

::

    auto_process.py archive

Copy the results to the final archive location:

::

    auto_process.py archive --final

