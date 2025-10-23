#!/usr/bin/env python3
#
#     bcl2fastq.apps.py: Bcl to Fastq conversion command line apps
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Provides methods to create ``Command`` instances for command line
applications used in Bcl to Fastq conversion.

The following functions are available:

* configureBclToFastq
* bcl2fastq2
* bclconvert
"""

from ..command import Command
from .utils import get_nmismatches

#######################################################################
# Functions
#######################################################################

def configureBclToFastq(basecalls_dir, sample_sheet, output_dir="Unaligned",
                        mismatches=None,
                        bases_mask=None,
                        force=False,
                        ignore_missing_bcl=False,
                        ignore_missing_stats=False,
                        ignore_missing_control=False,
                        configure_bcl_to_fastq_exe=None):
    """Generate Command instance for 'configureBclToFastq.pl' script

    Creates a Command instance to run the CASAVA 'configureBclToFastq.pl'
    script (which generates a Makefile to perform the bcl to fastq
    conversion).

    Arguments:
      basecalls_dir: path to the top-level directory holding the bcl
        files (typically 'Data/Intensities/Basecalls/' subdirectory)
      sample_sheet: path to the sample sheet file to use
      output_dir: optional, path to the output directory. Defaults to
        'Unaligned'. If this directory already exists then the
        conversion will fail unless the force option is set to True
      mismatches: optional, specify maximum number of mismatched bases
        allowed for matching index sequences during multiplexing.
        Recommended values are zero for indexes shorter than 6 base
        pairs, 1 for indexes of 6 or longer
        (If not specified and bases_mask is supplied then mismatches
        will be derived automatically from the bases mask string)
      bases_mask: optional, specify string indicating how to treat
        each cycle within each read e.g. 'y101,I6,y101'
      force: optional, if True then force overwrite of an existing
        output directory (default is False)
      ignore_missing_bcl: optional, if True then interpret missing bcl
        files as no call (default is False)
      ignore_missing_stats: optional, if True then fill in with zeroes
        when *.stats files are missing (default is False)
      ignore_missing_control: optional, if True then interpret missing
        control files as not-set control bits (default is False)
      configure_bcl_to_fastq_exe: optional, if set then will be taken
        as the name/path for the 'configureBclToFastq.pl' script

    Returns:
      Command object.
    """
    if configure_bcl_to_fastq_exe is None:
        configure_bcl_to_fastq_exe = 'configureBclToFastq.pl'
    configure_cmd = Command(configure_bcl_to_fastq_exe,
                            '--input-dir', basecalls_dir,
                            '--output-dir', output_dir,
                            '--sample-sheet', sample_sheet,
                            '--fastq-cluster-count', '-1')
    if bases_mask is not None:
        configure_cmd.add_args('--use-bases-mask', bases_mask)
    if mismatches is not None:
        configure_cmd.add_args('--mismatches', mismatches)
    else:
        # Nmismatches not supplied, derive from bases mask
        if bases_mask is not None:
            configure_cmd.add_args('--mismatches', get_nmismatches(bases_mask))
    if force:
        configure_cmd.add_args('--force')
    if ignore_missing_bcl:
        configure_cmd.add_args('--ignore-missing-bcl')
    if ignore_missing_stats:
        configure_cmd.add_args('--ignore-missing-stats')
    if ignore_missing_control:
        configure_cmd.add_args('--ignore-missing-control')
    return configure_cmd


def bcl2fastq2(run_dir, sample_sheet, output_dir="Unaligned",
               mismatches=None,
               bases_mask=None,
               ignore_missing_bcls=False,
               no_lane_splitting=False,
               minimum_trimmed_read_length=None,
               mask_short_adapter_reads=None,
               create_fastq_for_index_reads=False,
               find_adapters_with_sliding_window=False,
               loading_threads=None,
               demultiplexing_threads=None,
               processing_threads=None,
               writing_threads=None,
               bcl2fastq_exe=None):
    """
    Generate Command instance for 'bcl2fastq' program (v2.*)

    Creates a Command instance to run the Illumina 'bcl2fastq'
    program (for versions 2.*).

    Arguments:
      run_dir: path to the top-level directory for the run
      sample_sheet: path to the sample sheet file to use
      output_dir: optional, path to the output directory. Defaults to
        'Unaligned'
      mismatches: optional, specify maximum number of mismatched bases
        allowed for matching index sequences during multiplexing.
        Recommended values are zero for indexes shorter than 6 base
        pairs, 1 for indexes of 6 or longer
        (If not specified and bases_mask is supplied then mismatches
        will be derived automatically from the bases mask string)
      bases_mask: optional, specify string indicating how to treat
        each cycle within each read e.g. 'y101,I6,y101'
      ignore_missing_bcls: optional, if True then interpret missing bcl
        files as no call (default is False)
      no_lane_splitting: optional, if True then don't split FASTQ
        files by lane (--no-lane-splitting) (default is False)
      minimum_trimmed_read_length: optional, specify minimum length
        for reads after adapter trimming (shorter reads will be padded
        with Ns to make them long enough)
      mask_short_adapter_reads: optional, specify the minimum length
        of ACGT bases that must be present in a read after adapter
        trimming for it not to be masked completely with Ns.
      create_fastq_for_index_reads: optional, if True then also create
        Fastq files for index reads (default, don't create index read
        Fastqs) (--create-fastq-for-index-reads)
      find_adapters_with_sliding_window: optional, if True then
        use the sliding window algorithm rather than string matching
        when identifying adapter sequences for trimming (default,
        don't use sliding window algorithm)
        (--find-adapters-with-sliding-window)
      loading_threads: optional, specify number of threads to use
        for loading bcl data (--loading-threads)
      demultiplexing_threads: optional, specify number of threads to
        use for demultiplexing (--demultiplexing-threads)
      processing_threads: optional, specify number of threads to
        use for processing (--processing-threads)
      writing_threads: optional, specify number of threads to
        use for writing FASTQ data (--writing-threads)
      bcl2fastq_exe: optional, if set then specifies the name/path
        of the bcl2fastq executable to use

    Returns:
      Command object.
    """
    if bcl2fastq_exe is None:
        bcl2fastq_exe = 'bcl2fastq'
    bcl2fastq_cmd = Command(bcl2fastq_exe,
                            '--runfolder-dir', run_dir,
                            '--output-dir', output_dir,
                            '--sample-sheet', sample_sheet)
    if bases_mask is not None:
        bcl2fastq_cmd.add_args('--use-bases-mask', bases_mask)
    if mismatches is not None:
        bcl2fastq_cmd.add_args('--barcode-mismatches', mismatches)
    else:
        # Nmismatches not supplied, derive from bases mask
        if bases_mask is not None:
            bcl2fastq_cmd.add_args('--barcode-mismatches',
                                   get_nmismatches(bases_mask))
    if ignore_missing_bcls:
        bcl2fastq_cmd.add_args('--ignore-missing-bcls')
    if no_lane_splitting:
        bcl2fastq_cmd.add_args('--no-lane-splitting')
    if minimum_trimmed_read_length is not None:
        bcl2fastq_cmd.add_args('--minimum-trimmed-read-length',
                               minimum_trimmed_read_length)
    if mask_short_adapter_reads is not None:
        bcl2fastq_cmd.add_args('--mask-short-adapter-reads',
                               mask_short_adapter_reads)
    if create_fastq_for_index_reads:
        bcl2fastq_cmd.add_args('--create-fastq-for-index-read')
    if find_adapters_with_sliding_window:
        bcl2fastq_cmd.add_args('--find-adapters-with-sliding-window')
    if loading_threads is not None:
        bcl2fastq_cmd.add_args('-r', loading_threads)
    if demultiplexing_threads is not None:
        bcl2fastq_cmd.add_args('-d', demultiplexing_threads)
    if processing_threads is not None:
        bcl2fastq_cmd.add_args('-p', processing_threads)
    if writing_threads is not None:
        bcl2fastq_cmd.add_args('-w', writing_threads)
    return bcl2fastq_cmd


def bclconvert(run_dir, output_dir, sample_sheet=None,
               lane=None, no_lane_splitting=False,
               sampleproject_subdirectories=False,
               num_parallel_tiles=None,
               num_conversion_threads=None,
               num_compression_threads=None,
               num_decompression_threads=None,
               bclconvert_exe=None):
    """
    Generate Command instance for 'bcl-convert' program (v3.*)

    Creates a Command instance to run the Illumina 'bcl-convert'
    program (for versions 3.*).

    Arguments:
      run_dir: path to the top-level directory for the run
      output_dir: path to the output directory
      sample_sheet: optional, path to the sample sheet file to use
        (must be present in top-level of input directory if not
        specified here)
      lane (integer): restrict processing to single lane (sample
        sheet must only contain this lane) (--bcl-only-lane)
      no_lane_splitting: optional, if True then don't split FASTQ
        files by lane (--no-lane-splitting) (default is False)
      sampleproject_subdirectories: optional, if True then create
        subdirectories with project names in output (default is
        False) (--bcl-sampleproject-subdirectories)
      num_parallel_tiles: optional, specify the number of tiles
        being converted to Fastqs in parallel
        (--bcl-num-parallel-tiles)
      num_conversion_threads: optional, specify the number of threads
        to use for conversion per tile (--bcl-num-conversion-threads)
      num_compression_threads: optional, specify the number of
        threads for compressing output Fastq files
        (--bcl-num-compression-threads)
      num_decompression_threads: optional, specify the number of
        threads for decompression input bcl files
        (--bcl-num-decompression-threads)
      bclconvert_exe: optional, if set then specifies the name/path
        of the bcl-convert executable to use

    Returns:
      Command object.
    """
    if bclconvert_exe is None:
        bclconvert_exe = 'bcl-convert'
    bclconvert_cmd = Command(bclconvert_exe,
                             '--bcl-input-directory', run_dir,
                             '--output-dir', output_dir)
    if sample_sheet is not None:
        bclconvert_cmd.add_args('--sample-sheet', sample_sheet)
    if lane is not None:
        bclconvert_cmd.add_args('--bcl-only-lane', lane)
    if no_lane_splitting:
        bclconvert_cmd.add_args('--no-lane-splitting',
                                'true')
    if sampleproject_subdirectories:
        bclconvert_cmd.add_args('--bcl-sampleproject-subdirectories',
                                'true')
    if num_parallel_tiles is not None:
        bclconvert_cmd.add_args('--bcl-num-parallel-tiles',
                                num_parallel_tiles)
    if num_conversion_threads is not None:
        bclconvert_cmd.add_args('--bcl-num-conversion-threads',
                                num_conversion_threads)
    if num_compression_threads is not None:
        bclconvert_cmd.add_args('--bcl-num-compression-threads',
                                num_compression_threads)
    if num_decompression_threads is not None:
        bclconvert_cmd.add_args('--bcl-num-decompression-threads',
                                num_decompression_threads)
    return bclconvert_cmd