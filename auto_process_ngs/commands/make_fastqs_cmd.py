#!/usr/bin/env python
#
#     make_fastqs_cmd.py: implement auto process make_fastqs command
#     Copyright (C) University of Manchester 2018-2025 Peter Briggs
#
#########################################################################

#######################################################################
# Imports
#######################################################################

import os
import shutil
import gzip
import ast
import logging
from ..bcl2fastq.pipeline import PROTOCOLS
from ..bcl2fastq.pipeline import MakeFastqs
from bcftbx.IlluminaData import SampleSheet

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

BCL2FASTQ_DEFAULTS = {
    "bcl_converter": 'bcl2fastq',
    "minimum_trimmed_read_length": 35,
    "mask_short_adapter_reads": 22,
}

#######################################################################
# Command functions
#######################################################################

def make_fastqs(ap,protocol='standard',platform=None,
                unaligned_dir=None,sample_sheet=None,
                name=None,lanes=None,lane_subsets=None,
                icell8_well_list=None,
                nprocessors=None,bcl_converter=None,
                bases_mask=None,r1_length=None,r2_length=None,
                r3_length=None,no_lane_splitting=None,
                minimum_trimmed_read_length=None,
                mask_short_adapter_reads=None,
                trim_adapters=True,
                adapter_sequence=None,
                adapter_sequence_read2=None,
                create_fastq_for_index_read=None,
                find_adapters_with_sliding_window=None,
                generate_stats=True,stats_file=None,
                per_lane_stats_file=None,
                analyse_barcodes=True,barcode_analysis_dir=None,
                force_copy_of_primary_data=False,
                create_empty_fastqs=False,
                ignore_missing_bcls=False,runner=None,
                icell8_swap_i1_and_i2=False,
                icell8_reverse_complement=None,
                cellranger_jobmode=None,
                cellranger_mempercore=None,
                cellranger_maxjobs=None,
                cellranger_jobinterval=None,
                cellranger_localcores=None,
                cellranger_localmem=None,
                cellranger_ignore_dual_index=False,
                spaceranger_rc_i2_override=None,
                max_jobs=None,max_cores=None,batch_limit=None,
                enable_conda=None,conda_env_dir=None,
                use_conda_for_bcl2fastq=None,
                verbose=False,working_dir=None):
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
      name (str): (optional) identifier for outputs that are not
        set explicitly
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
      generate_stats (bool): if True then (re)generate statistics file
        for fastqs
      analyse_barcodes (bool): if True then (re)analyse barcodes for
        fastqs
      bcl_converter (str): default BCL-to-Fastq conversion software to
        use; optionally can include a version specification (e.g.
        "bcl2fastq>2.0" or "bcl-convert=3.7.5"). Defaults to "bcl2fastq"
      bases_mask (str): if set then use this as an alternative bases
        mask setting
      r1_length (int): explicitly specify length to truncate R1 reads
         to (ignored if bases mask is set)
      r2_length (int): explicitly specify length to truncate R2 reads
         to (ignored if bases mask is set, or if there is no R2 read)
      r3_length (int): explicitly specify length to truncate R3 reads
         to (ignored if bases mask is set, or if there is no R2 read)
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
      create_fastq_for_index_reads (bool): if True then also create
        Fastq files for index reads (default, don't create index read
        Fastqs)
      ignore_missing_bcls (bool): if True then tell BCL conversion
        software to ignore missing or corrupted BCLs (default: False,
        don't ignore missing or corrupted BCL files)
      find_adapters_with_sliding_window (boolean): if True then use
        sliding window algorithm to identify adapter sequences for
        trimming
      stats_file (str): if set then use this as the name of the output
        per-fastq stats file.
      per_lane_stats_file (str): if set then use this as the name of
        the output per-lane stats file.
      barcode_analysis_dir (str): if set then specifies path to the
        output directory for barcode analysis
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
      spaceranger_rc_i2_override (bool): (optional) if set then value
         is passed to Spaceranger's '--rc-i2-override' option (True for
         reverse complement workflow B, False for forward complement
         workflow A). If not set then Spaceranger will be left to
         determine the workflow automatically
      max_jobs (int): maximum number of concurrent jobs allowed
      max_cores (int): maximum number of cores available
      batch_limit (int): if set then run commands in each task in
         batches, with the batch size set dyanmically so as not to
         exceed this limit
      working_dir (str): path to a working directory (defaults to
         temporary directory in the current directory)
      enable_conda (bool): if True then use conda to resolve
        dependencies declared on tasks in the pipeline
      conda_env_dir (str): path to non-default directory for conda
        environments
      use_conda_for_bcl2fastq (bool): if True then use conda
        packages for 'bcl2fastq' dependency resolution (NB ignored
        unless 'enable_conda' is also True)
      verbose (bool): if True then report additional information for
         pipeline diagnostics
    """
    # Report protocol
    print("Protocol              : %s" % protocol)
    if protocol not in PROTOCOLS:
        raise Exception("Unknown protocol: '%s' (must be one of "
                        "%s)" % (protocol,','.join(PROTOCOLS)))

    # Output (unaligned) dir
    if not unaligned_dir and name:
        unaligned_dir = 'bcl2fastq_%s' % name
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
        s = SampleSheet(ap.params.sample_sheet)
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
    if not barcode_analysis_dir and name:
        barcode_analysis_dir = 'barcode_analysis_%s' % name
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
        if name:
            stats_file='statistics.%s.info' % name
        elif ap.params['stats_file'] is not None:
            stats_file = ap.params['stats_file']
        else:
            stats_file='statistics.info'
    if per_lane_stats_file is None:
        if name:
            per_lane_stats_file='per_lane_statistics.%s.info' % name
        elif ap.params['per_lane_stats_file'] is not None:
            per_lane_stats_file = ap.params['per_lane_stats_file']
        else:
            per_lane_stats_file='per_lane_statistics.info'

    # Log dir
    log_dir = 'make_fastqs'
    if name:
        log_dir += "_%s" % name
    if protocol != 'standard':
        log_dir += "_%s" % protocol
    if lanes:
        log_dir += "_L%s" % ''.join([str(l) for l in sorted(lanes)])
    ap.set_log_dir(ap.get_log_subdir(log_dir))

    # Pipeline log file
    pipeline_log = os.path.join(ap.log_dir,"make_fastqs.log")

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

    # Get platform
    if not platform:
        platform = ap.metadata.platform

    # Set options from supplied arguments, platform-specific settings
    # and configured defaults
    defaults = {
        'bcl_converter': bcl_converter,
        'nprocessors': nprocessors,
        'no_lane_splitting': no_lane_splitting,
        'create_empty_fastqs': create_empty_fastqs,
        'ignore_missing_bcls': ignore_missing_bcls
    }
    for item in ('bcl_converter',
                 'nprocessors',
                 'no_lane_splitting',
                 'create_empty_fastqs',):
        if defaults[item] is None:
            value = None
            if platform in ap.settings.platform:
                value = ap.settings.platform[platform][item]
            if value is None:
                value = ap.settings.bcl_conversion[item]
            defaults[item] = value

    # BCL converter
    if defaults['bcl_converter'] is None:
        bcl_converter = BCL2FASTQ_DEFAULTS['bcl_converter']
    else:
        bcl_converter = defaults['bcl_converter']

    # Number of processors
    nprocessors = defaults['nprocessors']

    # Split Fastqs across lanes
    no_lane_splitting = defaults['no_lane_splitting']

    # Create empty Fastqs
    create_empty_fastqs = defaults['create_empty_fastqs']

    # Ignore missing/corrupted BCL files
    ignore_missing_bcls = defaults['ignore_missing_bcls']

    # Set up pipeline runners
    default_runner = ap.settings.general.default_runner
    runners = {
        'rsync_runner': ap.settings.runners.rsync,
        'barcode_analysis_runner': ap.settings.runners.barcode_analysis,
        'bcl2fastq_runner': ap.settings.runners.bcl2fastq,
        'bclconvert_runner': ap.settings.runners.bcl_convert,
        'demultiplex_icell8_atac_runner': ap.settings.runners.bcl2fastq,
        'cellranger_runner': ap.settings.runners.cellranger_mkfastq,
        'cellranger_atac_runner': ap.settings.runners.cellranger_mkfastq,
        'cellranger_arc_runner': ap.settings.runners.cellranger_mkfastq,
        'merge_fastqs_runner': ap.settings.runners.merge_fastqs,
        'spaceranger_runner': ap.settings.runners.cellranger_mkfastq,
        'stats_runner': ap.settings.runners.stats,
    }
    if runner is not None:
        # Override configured runners
        default_runner = runner
        for r in runners:
            runner[r] = runner

    # Set up pipeline environment modules
    envmodules = {}
    for envmod in ('bcl2fastq',
                   'bcl_convert',
                   'cellranger_mkfastq',
                   'cellranger_atac_mkfastq',
                   'cellranger_arc_mkfastq',
                   'spaceranger_mkfastq',):
        try:
            envmodules[envmod] = ap.settings.modulefiles[envmod]
        except KeyError:
            try:
                envmodules[envmod] = ap.settings.modulefiles['make_fastqs']
            except KeyError:
                envmodules[envmod] = None

    # Conda dependency resolution
    if enable_conda is None:
        enable_conda = ap.settings.conda.enable_conda
    if conda_env_dir is None:
        conda_env_dir = ap.settings.conda.env_dir

    # Other pipeline settings
    poll_interval = ap.settings.general.poll_interval

    # Construct and run pipeline
    make_fastqs = MakeFastqs(ap.params.data_dir,
                             ap.params.sample_sheet,
                             protocol=protocol,
                             bases_mask=bases_mask,
                             r1_length=r1_length,
                             r2_length=r2_length,
                             r3_length=r3_length,
                             bcl_converter=bcl_converter,
                             platform=platform,
                             icell8_well_list=icell8_well_list,
                             minimum_trimmed_read_length=\
                             minimum_trimmed_read_length,
                             mask_short_adapter_reads=\
                             mask_short_adapter_reads,
                             adapter_sequence=adapter_sequence,
                             adapter_sequence_read2=\
                             adapter_sequence_read2,
                             spaceranger_rc_i2_override=\
                             spaceranger_rc_i2_override,
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
                             name=name,
                             out_dir=ap.params.unaligned_dir,
                             barcode_analysis_dir=barcode_analysis_dir,
                             primary_data_dir=ap.params.primary_data_dir,
                             force_copy_of_primary_data=\
                             force_copy_of_primary_data,
                             no_lane_splitting=no_lane_splitting,
                             create_fastq_for_index_read=\
                             create_fastq_for_index_read,
                             find_adapters_with_sliding_window=\
                             find_adapters_with_sliding_window,
                             create_empty_fastqs=create_empty_fastqs,
                             ignore_missing_bcls=ignore_missing_bcls,
                             stats_file=stats_file,
                             per_lane_stats=per_lane_stats_file,
                             nprocessors=nprocessors,
                             default_runner=default_runner,
                             cellranger_jobmode=cellranger_jobmode,
                             cellranger_mempercore=cellranger_mempercore,
                             cellranger_maxjobs=cellranger_maxjobs,
                             cellranger_jobinterval=\
                             cellranger_jobinterval,
                             cellranger_localcores=cellranger_localcores,
                             cellranger_localmem=cellranger_localmem,
                             runners=runners,
                             enable_conda=enable_conda,
                             conda_env_dir=conda_env_dir,
                             use_conda_for_bcl2fastq=use_conda_for_bcl2fastq,
                             envmodules=envmodules,
                             log_dir=ap.log_dir,
                             log_file=pipeline_log,
                             max_jobs=max_jobs,
                             max_slots=max_cores,
                             batch_limit=batch_limit,
                             poll_interval=poll_interval,
                             working_dir=working_dir,
                             verbose=verbose)

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
        # Platform
        ap.metadata['platform'] = make_fastqs.output.platform
        # Flow cell mode
        ap.metadata['flow_cell_mode'] = make_fastqs.output.flow_cell_mode
        # Software used for processing
        try:
            processing_software = ast.literal_eval(
                ap.metadata.processing_software)
        except ValueError:
            processing_software = dict()
        outputs = make_fastqs.output
        if outputs.bcl2fastq_info:
            processing_software['bcl2fastq'] = outputs.bcl2fastq_info
        if outputs.bclconvert_info:
            processing_software['bcl-convert'] = outputs.bclconvert_info
        if outputs.cellranger_info:
            processing_software['cellranger'] = outputs.cellranger_info
        if outputs.cellranger_atac_info:
            processing_software['cellranger-atac'] = \
                                                outputs.cellranger_atac_info
        if outputs.cellranger_arc_info:
            processing_software['cellranger-arc'] = \
                                                outputs.cellranger_arc_info
        if outputs.spaceranger_info:
            processing_software['spaceranger'] = outputs.spaceranger_info
        ap.metadata['processing_software'] = processing_software
        # Legacy metadata items
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
    if not ap.params.project_metadata:
        ap.make_project_metadata_file()
    else:
        ap.update_project_metadata_file()
    ap.save_data()

    # Finish
    return status
