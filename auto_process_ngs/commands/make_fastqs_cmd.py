#!/usr/bin/env python
#
#     make_fastqs_cmd.py: implement auto process make_fastqs command
#     Copyright (C) University of Manchester 2018-2020 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import time
import shutil
import gzip
import logging
from ..bcl2fastq.pipeline import PROTOCOLS
from ..bcl2fastq.pipeline import MakeFastqs

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

BCL2FASTQ_VERSIONS = ('1.8','2.17','2.20',)

BCL2FASTQ_DEFAULTS = {
    "minimum_trimmed_read_length": 35,
    "mask_short_adapter_reads": 22,
}

#######################################################################
# Command functions
#######################################################################

def make_fastqs(ap,protocol='standard',platform=None,
                unaligned_dir=None,sample_sheet=None,
                lanes=None,lane_subsets=None,
                icell8_well_list=None,
                ignore_missing_bcl=False,ignore_missing_stats=False,
                skip_rsync=False,remove_primary_data=False,
                nprocessors=None,require_bcl2fastq_version=None,
                bases_mask=None,no_lane_splitting=None,
                minimum_trimmed_read_length=None,
                mask_short_adapter_reads=None,
                trim_adapters=True,
                adapter_sequence=None,
                adapter_sequence_read2=None,
                create_fastq_for_index_read=False,
                generate_stats=True,stats_file=None,
                per_lane_stats_file=None,
                analyse_barcodes=True,barcode_analysis_dir=None,
                skip_fastq_generation=False,
                only_fetch_primary_data=False,
                force_copy_of_primary_data=False,
                create_empty_fastqs=None,runner=None,
                icell8_swap_i1_and_i2=False,
                icell8_reverse_complement=None,
                cellranger_jobmode=None,
                cellranger_mempercore=None,
                cellranger_maxjobs=None,
                cellranger_jobinterval=None,
                cellranger_localcores=None,
                cellranger_localmem=None,
                cellranger_ignore_dual_index=False):
    """
    Create and summarise FASTQ files

    Wrapper for operations related to FASTQ file generation and analysis.
    The operations are typically:
 
    - get primary data (BCL files)
    - run bcl-to-fastq conversion
    - generate statistics
    - analyse barcodes

    If the number of processors and the job runner are not explicitly
    specified then these are taken from the settings for the bcl2fastq
    and the statistics generation steps, which may differ from each other.
    However if either of these values are set explicitly then the same
    values will be used for both steps.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      protocol (str): if set then specifies the protocol to use
        for fastq generation, otherwise use the 'standard' bcl2fastq
        protocol
      platform (str): if set then specifies the sequencing platform
        (otherwise platform will be determined from the primary data)
      unaligned_dir (str): if set then use this as the output directory
        for bcl-to-fastq conversion. Default is 'bcl2fastq' (unless
        an alternative is already specified in the config file)
      sample_sheet (str): if set then use this as the input samplesheet
      lanes (list): (optional) specify a list of lane numbers to
        use in the processing; lanes not in the list will be excluded
        (default is to include all lanes)
      lane_subsets (list): (optional) specify a list of lane subsets
        to process separately before merging at the end; each subset
        is a dictionary which should be generated using the 'subset'
        function, and can include custom values for processing
        parameters (e.g. protocol, trimming and masking options etc)
        to override the defaults for this lane. Lanes not in a subset
        will still be processed unless excluded via the 'lanes'
        keyword
      icell8_well_list (str): well list file for ICELL8 platforms
        (required for ICELL8 processing protocols)
      nprocessors (int) : number of processors to use
      ignore_missing_bcl (bool): if True then run bcl2fastq with
        --ignore-missing-bcl
      ignore_missing_stats (bool): if True then run bcl2fastq with
        --ignore-missing-stats
      skip_rsync (bool): if True then don't rsync primary data at the
        start of bcl2fastq conversion
      remove_primary_data (bool): if True then remove primary data at
        the end of bcl2fastq conversion (default is to keep it)
      generate_stats (bool): if True then (re)generate statistics file
        for fastqs
      analyse_barcodes (bool): if True then (re)analyse barcodes for
        fastqs
      require_bcl2fastq_version (str): (optional) specify bcl2fastq
        version to use. Should be a string of the form '1.8.4' or
        '>2.0'. Set to None to automatically determine required
        bcl2fastq version.
      bases_mask (str): if set then use this as an alternative bases
        mask setting
      no_lane_splitting (bool): if True then run bcl2fastq with
        --no-lane-splitting
      minimum_trimmed_read_length (int): if set then specify minimum
        length for reads after adapter trimming (shorter reads will
        be padded with Ns to make them long enough)
      mask_short_adapter_reads (int): if set then specify the minimum
        length of ACGT bases that must be present in a read after
        adapter trimming for it not to be masked completely
        with Ns.
      trim_adapters (boolean): if True (the default) then pass
        adapter sequence(s) to bcl2fastq to perform adapter trimming;
        otherwise remove adapter sequences
      adapter_sequence (str): if not None then specifies adapter
        sequence to use instead of any sequences already set in the
        samplesheet (nb will be ignored if 'trim_adapters' is False)
      adapter_sequence_read2 (str): if not None then specifies adapter
        sequence to use for read2 instead of any sequences already set
        in the samplesheet (nb will be ignored if 'trim_adapters' is
        False)
      create_fastq_for_index_reads (boolean): if True then also create
        Fastq files for index reads (default, don't create index read
        Fastqs)
      stats_file (str): if set then use this as the name of the output
        per-fastq stats file.
      per_lane_stats_file (str): if set then use this as the name of
        the output per-lane stats file.
      barcode_analysis_dir (str): if set then specifies path to the
        output directory for barcode analysis
      skip_fastq_generation (bool): if True then don't perform fastq
        generation
      only_fetch_primary_data (bool): if True then fetch primary data,
        don't do anything else
      force_copy_of_primary_data (bool): if True then force primary
        data to be copied (rsync'ed) even if it's on the local system
        (default is to link to primary data unless it's on a remote
        filesystem).
      create_empty_fastqs (bool): if True then create empty 'placeholder'
        fastq files for any missing fastqs after bcl2fastq
        (must have completed with zero exit status)
      runner (JobRunner): (optional) specify a non-default job runner
        to use for fastq generation
      icell8_swap_i1_and_i2 (bool): if True then swap I1 and I2 reads
        when matching to barcodes in the ICELL8 well list (ICELL8 ATAC
        data only)
      icell8_reverse_complement (str): one of 'i1', 'i2', 'both', or
        None; if set then the specified index reads will be reverse
        complemented when matching to barcodes in the ICELL8 well list
        (ICELL8 ATAC data only)
      cellranger_jobmode (str): (optional) job mode to run cellranger in
        (10xGenomics Chromium SC data only)
      cellranger_mempercore (int): (optional) memory assumed per core
        (in Gbs) (10xGenomics Chromium SC data only)
      cellranger_maxjobs (int): (optional) maxiumum number of concurrent
         jobs to run (10xGenomics Chromium SC data only)
      cellranger_jobinterval (int): (optional) how often jobs are
         submitted (in ms) (10xGenomics Chromium SC data only)
      cellranger_localcores (int): (optional) maximum number of cores
         cellranger can request in jobmode 'local' (10xGenomics Chromium
         SC data only)
      cellranger_localmem (int): (optional) maximum memory cellranger
         can request in jobmode 'local' (10xGenomics Chromium SC data
         only)
      cellranger_ignore_dual_index (bool): (optional) on a dual-indexed
         flowcell where the second index was not used for the 10x
         sample, ignore it (10xGenomics Chromium SC data only)
    """
    # Report protocol
    print("Protocol              : %s" % protocol)
    if protocol not in PROTOCOLS:
        raise Exception("Unknown protocol: '%s' (must be one of "
                        "%s)" % (protocol,','.join(PROTOCOLS)))

    # Trap for unsupported options
    if only_fetch_primary_data:
        raise Exception("'Only fetch primary data' not supported")
    if skip_fastq_generation:
        raise Exception("'Skip fastq generation' not supported")
    if ignore_missing_bcl:
        raise Exception("'Ignore missing bcl' not suppported")
    if ignore_missing_stats:
        raise Exception("'Ignore missing stats' not supported")
    if remove_primary_data:
        raise Exception("'Remove primary data' not supported")

    # Output (unaligned) dir
    if unaligned_dir is not None:
        ap.params['unaligned_dir'] = unaligned_dir
    elif ap.params['unaligned_dir'] is None:
        ap.params['unaligned_dir'] = 'bcl2fastq'
    print("Output dir            : %s" % ap.params.unaligned_dir)

    # Sample sheet
    if sample_sheet is None:
        sample_sheet = ap.params.sample_sheet
    if not os.path.isabs(sample_sheet):
        sample_sheet = os.path.join(ap.analysis_dir,sample_sheet)
    if not os.path.isfile(sample_sheet):
        raise Exception("Missing sample sheet '%s'" % sample_sheet)
    ap.params['sample_sheet'] = sample_sheet
    print("Source sample sheet   : %s" % ap.params.sample_sheet)

    # Check requested lanes are actually present
    print("Lanes                 : %s" % ('all' if lanes is None
                                          else
                                          ','.join([str(l)
                                                    for l in lanes])))
    if lanes is not None:
        s = IlluminaData.SampleSheet(ap.params.sample_sheet)
        if not s.has_lanes:
            raise Exception("Requested subset of lanes but "
                            "samplesheet doesn't contain any "
                            "lane information")
        samplesheet_lanes = list(set([l['Lane'] for l in s]))
        for l in lanes:
            if l not in samplesheet_lanes:
                raise Exception("Requested lane '%d' not present "
                                "in samplesheet" % l)

    # Barcode analysis
    if barcode_analysis_dir is not None:
        ap.params['barcode_analysis_dir'] = barcode_analysis_dir
    elif ap.params.barcode_analysis_dir is None:
        ap.params['barcode_analysis_dir'] = 'barcode_analysis'
    barcode_analysis_dir = ap.params.barcode_analysis_dir
    if not os.path.isabs(barcode_analysis_dir):
        barcode_analysis_dir = os.path.join(ap.params.analysis_dir,
                                            barcode_analysis_dir)

    # Statistics files
    if stats_file is None:
        if ap.params['stats_file'] is not None:
            stats_file = ap.params['stats_file']
        else:
            stats_file='statistics.info'
    if per_lane_stats_file is None:
        if ap.params['per_lane_stats_file'] is not None:
            per_lane_stats_file = ap.params['per_lane_stats_file']
        else:
            per_lane_stats_file='per_lane_statistics.info'

    # Log dir
    log_dir = 'make_fastqs'
    if protocol != 'standard':
        log_dir += "_%s" % protocol
    if lanes:
        log_dir += "_L%s" % ''.join([str(l) for l in sorted(lanes)])
    ap.set_log_dir(ap.get_log_subdir(log_dir))

    # Pipeline log file
    pipeline_log = os.path.join(ap.log_dir,"make_fastqs.log")

    # Deal with platform information
    if not platform:
        platform = ap.metadata.platform
    # Bases mask
    if bases_mask is not None:
        ap.params['bases_mask'] = bases_mask
    bases_mask = ap.params.bases_mask

    # Default trimming/masking values
    if minimum_trimmed_read_length is None:
        minimum_trimmed_read_length = \
                BCL2FASTQ_DEFAULTS['minimum_trimmed_read_length']
    if mask_short_adapter_reads is None:
        mask_short_adapter_reads = \
                BCL2FASTQ_DEFAULTS['mask_short_adapter_reads']

    # Require specific bcl2fastq version
    if require_bcl2fastq_version is None:
        # Look for platform-specific requirement
        try:
            require_bcl2fastq_version = \
                ap.settings.platform[ap.metadata.platform].bcl2fastq
            print("Bcl2fastq version %s required for platform '%s'" %
                  (ap.metadata.platform,require_bcl2fastq_version))
        except (KeyError,AttributeError):
            pass
    if require_bcl2fastq_version is None:
        # Look for default requirement
        require_bcl2fastq_version = ap.settings.bcl2fastq.default_version
        print("Bcl2fastq version %s required by default" %
              require_bcl2fastq_version)
    if require_bcl2fastq_version is not None:
        # No version requirement
        print("No bcl2fastq version explicitly specified")

    # Set up pipeline runners
    runners = {
        'rsync_runner': ap.settings.runners.rsync,
        'bcl2fastq_runner': ap.settings.runners.bcl2fastq,
        'demultiplex_icell8_atac_runner': ap.settings.runners.bcl2fastq,
        'cellranger_runner': ap.settings.runners.cellranger,
        'cellranger_atac_runner': ap.settings.runners.cellranger,
        'stats_runner': ap.settings.runners.stats,
    }
    if runner is not None:
        # Override configured runners
        for r in runners:
            runner[r] = runner

    # Set up pipeline environment modules
    envmodules = {}
    for name in ('bcl2fastq',
                 'cellranger_mkfastq',
                 'cellranger_atac_mkfastq',):
        try:
            envmodules[name] = ap.settings.modulefiles[name]
        except KeyError:
            try:
                envmodules[name] = ap.settings.modulefiles['make_fastqs']
            except KeyError:
                envmodules[name] = None

    # Other pipeline settings
    poll_interval = ap.settings.general.poll_interval

    # Construct and run pipeline
    make_fastqs = MakeFastqs(ap.params.data_dir,
                             ap.params.sample_sheet,
                             protocol=protocol,
                             bases_mask=bases_mask,
                             platform=platform,
                             icell8_well_list=icell8_well_list,
                             minimum_trimmed_read_length=\
                             minimum_trimmed_read_length,
                             mask_short_adapter_reads=\
                             mask_short_adapter_reads,
                             adapter_sequence=adapter_sequence,
                             adapter_sequence_read2=\
                             adapter_sequence_read2,
                             icell8_atac_swap_i1_and_i2=\
                             icell8_swap_i1_and_i2,
                             icell8_atac_reverse_complement=\
                             icell8_reverse_complement,
                             lane_subsets=lane_subsets,
                             lanes=lanes,
                             trim_adapters=trim_adapters,
                             fastq_statistics=generate_stats,
                             analyse_barcodes=analyse_barcodes)
    status = make_fastqs.run(ap.analysis_dir,
                             out_dir=ap.params.unaligned_dir,
                             barcode_analysis_dir=barcode_analysis_dir,
                             primary_data_dir=ap.params.primary_data_dir,
                             force_copy_of_primary_data=\
                             force_copy_of_primary_data,
                             no_lane_splitting=no_lane_splitting,
                             create_fastq_for_index_read=\
                             create_fastq_for_index_read,
                             create_empty_fastqs=create_empty_fastqs,
                             stats_file=stats_file,
                             per_lane_stats=per_lane_stats_file,
                             nprocessors=nprocessors,
                             require_bcl2fastq=require_bcl2fastq_version,
                             cellranger_jobmode=cellranger_jobmode,
                             cellranger_mempercore=cellranger_mempercore,
                             cellranger_maxjobs=cellranger_maxjobs,
                             cellranger_jobinterval=\
                             cellranger_jobinterval,
                             cellranger_localcores=cellranger_localcores,
                             cellranger_localmem=cellranger_localmem,
                             runners=runners,
                             envmodules=envmodules,
                             log_dir=ap.log_dir,
                             log_file=pipeline_log,
                             poll_interval=poll_interval)

    # Update the parameters
    ap.params['primary_data_dir'] = make_fastqs.output.primary_data_dir
    ap.params['acquired_primary_data'] = \
                                    make_fastqs.output.acquired_primary_data
    ap.params['stats_file'] = make_fastqs.output.stats_file
    ap.params['per_lane_stats_file'] = make_fastqs.output.per_lane_stats
    for param in ('stats_file','per_lane_stats_file'):
        filen = ap.params[param]
        if filen is not None:
            if filen.startswith(ap.analysis_dir):
                ap.params[param] = os.path.relpath(filen,ap.analysis_dir)

    # Update the metadata
    if status == 0:
        ap.metadata['platform'] = make_fastqs.output.platform
        ap.metadata['bcl2fastq_software'] = make_fastqs.output.bcl2fastq_info
        ap.metadata['cellranger_software'] = make_fastqs.output.cellranger_info

    # Make a file listing missing Fastqs
    if make_fastqs.output.missing_fastqs:
        missing_fastqs_log = os.path.join(ap.log_dir,
                                          "missing_fastqs.log")
        with open(missing_fastqs_log,'wt') as fp:
            for fq in make_fastqs.output.missing_fastqs:
                fp.write("%s\n" % fq)
        print("Wrote list of missing Fastq files to '%s'" %
              missing_fastqs_log)
    # Raise exception on failure
    if status != 0:
        raise Exception("Fastq generation failed")

    # Make or update 'projects.info' metadata file
    if lanes:
        ap.update_project_metadata_file()
    else:
        ap.make_project_metadata_file()

    # Finish
    return status
