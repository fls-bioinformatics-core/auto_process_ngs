#!/usr/bin/env python
#
#     make_fastqs_cmd.py: implement auto process make_fastqs command
#     Copyright (C) University of Manchester 2018-2019 Peter Briggs
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
from ..bcl2fastq_utils import make_custom_sample_sheet
from ..bcl2fastq_utils import get_bases_mask
from ..bcl2fastq_utils import bases_mask_is_valid
from ..bcl2fastq_utils import available_bcl2fastq_versions
from ..bcl2fastq_utils import bcl_to_fastq_info
from ..bcl2fastq_utils import get_required_samplesheet_format
from ..bcl2fastq_utils import get_nmismatches
from ..bcl2fastq_utils import check_barcode_collisions
from ..barcodes.pipeline import AnalyseBarcodes
from ..icell8.utils import get_bases_mask_icell8
from ..icell8.utils import get_bases_mask_icell8_atac
from ..icell8.utils import ICell8WellList
from ..tenx_genomics_utils import has_chromium_sc_indices
from ..tenx_genomics_utils import get_bases_mask_10x_atac
from ..tenx_genomics_utils import cellranger_info
from ..tenx_genomics_utils import run_cellranger_mkfastq
from ..samplesheet_utils import SampleSheetLinter
from ..applications import Command
from ..applications import general as general_apps
from ..simple_scheduler import SchedulerJob
from ..qc.processing import report_processing_qc
from bcftbx import IlluminaData
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

MAKE_FASTQS_PROTOCOLS = ('standard',
                         'icell8',
                         'icell8_atac',
                         '10x_chromium_sc',
                         '10x_chromium_sc_atac')

#######################################################################
# Command functions
#######################################################################

def make_fastqs(ap,protocol='standard',platform=None,
                unaligned_dir=None,sample_sheet=None,lanes=None,
                icell8_well_list=None,
                ignore_missing_bcl=False,ignore_missing_stats=False,
                skip_rsync=False,remove_primary_data=False,
                nprocessors=None,require_bcl2fastq_version=None,
                bases_mask=None,no_lane_splitting=None,
                minimum_trimmed_read_length=None,
                mask_short_adapter_reads=None,
                create_fastq_for_index_reads=False,
                generate_stats=True,stats_file=None,
                per_lane_stats_file=None,
                analyse_barcodes=True,barcode_analysis_dir=None,
                skip_fastq_generation=False,
                only_fetch_primary_data=False,
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
    """Create and summarise FASTQ files

    Wrapper for operations related to FASTQ file generation and analysis.
    The operations are typically:
 
    - get primary data (BCL files)
    - run bcl-to-fastq conversion
    - generate statistics

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
      icell8_well_list (str): well list file for ICELL8 platforms
        (required for ICELL8 processing protocols)
      nprocessors (int) : number of processors to run bclToFastq.py with
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
    print "Protocol              : %s" % protocol
    if protocol not in MAKE_FASTQS_PROTOCOLS:
        raise Exception("Unknown protocol: '%s' (must be one of "
                        "%s)" % (protocol,
                                 ','.join([MAKE_FASTQS_PROTOCOLS])))
    # Unaligned dir
    if unaligned_dir is not None:
        ap.params['unaligned_dir'] = unaligned_dir
    elif ap.params['unaligned_dir'] is None:
        ap.params['unaligned_dir'] = 'bcl2fastq'
    print "Output dir            : %s" % ap.params.unaligned_dir
    # Sample sheet
    if sample_sheet is None:
        sample_sheet = ap.params.sample_sheet
    if not os.path.isabs(sample_sheet):
        sample_sheet = os.path.join(ap.analysis_dir,sample_sheet)
    if not os.path.isfile(sample_sheet):
        raise Exception("Missing sample sheet '%s'" % sample_sheet)
    ap.params['sample_sheet'] = sample_sheet
    print "Source sample sheet   : %s" % ap.params.sample_sheet
    # Check requested lanes are actually present
    print "Lanes                 : %s" % ('all' if lanes is None
                                          else
                                          ','.join([str(l)
                                                    for l in lanes]))
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
    # Make a temporary sample sheet
    if lanes:
        lanes_id = ".L%s" % ''.join([str(l) for l in lanes])
    else:
        lanes_id = ""
    sample_sheet = os.path.join(ap.tmp_dir,
                                "SampleSheet%s.%s.csv" %
                                (lanes_id,
                                 time.strftime("%Y%m%d%H%M%S")))
    make_custom_sample_sheet(ap.params.sample_sheet,
                             sample_sheet,
                             lanes=lanes)
    # Check the temporary sample sheet
    print "Checking temporary sample sheet"
    invalid_barcodes = SampleSheetLinter(
        sample_sheet_file=sample_sheet).has_invalid_barcodes()
    if invalid_barcodes:
        logger.error("Invalid barcodes detected")
        for line in invalid_barcodes:
            logger.critical("%s" % line)
    invalid_characters = SampleSheetLinter(
        sample_sheet_file=sample_sheet).has_invalid_characters()
    if invalid_characters:
        logger.critical("Invalid non-printing/non-ASCII characters "
                        "detected")
    if invalid_barcodes or invalid_characters:
        raise Exception("Errors detected in generated sample sheet")
    # Adjust verification settings for 10xGenomics Chromium SC
    # data if necessary
    verify_include_sample_dir = False
    if has_chromium_sc_indices(sample_sheet):
        if protocol in ('10x_chromium_sc','10x_chromium_sc_atac',):
            # Force inclusion of sample-name subdirectories
            # when verifying Chromium SC data
            print "Sample sheet includes Chromium SC indices"
            verify_include_sample_dir = True
        else:
            # Chromium SC indices detected but not using
            # 10x_chromium_sc protocol
            raise Exception("Detected 10xGenomics Chromium SC indices "
                            "in generated sample sheet but protocol "
                            "'%s' has been specified; use an "
                            "appropriate '10x_...' protocol for these "
                            "indices" % protocol)
    # Check for pre-existing Fastq outputs
    if verify_fastq_generation(
            ap,
            unaligned_dir=ap.params.unaligned_dir,
            icell8_well_list=icell8_well_list,
            lanes=lanes,
            include_sample_dir=verify_include_sample_dir):
        print "Expected Fastq outputs already present"
        skip_rsync = True
        skip_fastq_generation = True
    # Check if there's anything to do
    if (skip_rsync and skip_fastq_generation) and \
       not (generate_stats or analyse_barcodes):
        print "Nothing to do"
        return
    # Log dir
    log_dir = 'make_fastqs'
    if protocol != 'standard':
        log_dir += "_%s" % protocol
    if lanes:
        log_dir += "_L%s" % ''.join([str(l) for l in sorted(lanes)])
    ap.set_log_dir(ap.get_log_subdir(log_dir))
    # Fetch primary data
    if not skip_rsync and not ap.params.acquired_primary_data:
        if get_primary_data(ap) != 0:
            logger.error("Failed to acquire primary data")
            raise Exception("Failed to acquire primary data")
        else:
            ap.params['acquired_primary_data'] = True
    if only_fetch_primary_data:
        return
    # Deal with platform information
    if not platform:
        platform = ap.metadata.platform
    # Do fastq generation using the specified protocol
    if not skip_fastq_generation:
        # Set primary data location and report info
        primary_data_dir = os.path.join(
            ap.params.primary_data_dir,
            os.path.basename(ap.params.data_dir))
        print "Primary data dir      : %s" % primary_data_dir
        try:
            illumina_run = IlluminaData.IlluminaRun(primary_data_dir,
                                                    platform=platform)
        except IlluminaData.IlluminaDataPlatformError as ex:
            logger.critical("Error loading primary data: %s" % ex)
            if platform is None:
                logger.critical("Try specifying platform using --platform?")
            else:
                logger.critical("Check specified platform is valid (or "
                                "omit --platform")
            raise Exception("Error determining sequencer platform")
        print "Platform              : %s" % illumina_run.platform
        print "Bcl format            : %s" % illumina_run.bcl_extension
        # Set platform in metadata
        ap.metadata['platform'] = illumina_run.platform
        # Bases mask
        if bases_mask is not None:
            ap.params['bases_mask'] = bases_mask
        bases_mask = ap.params.bases_mask
        print "Bases mask setting    : %s" % bases_mask
        if protocol not in ('icell8_atac',
                            '10x_chromium_sc',
                            '10x_chromium_sc_atac',):
            if bases_mask == "auto":
                print "Determining bases mask from RunInfo.xml"
                bases_mask = get_bases_mask(illumina_run.runinfo_xml,
                                            sample_sheet)
                if not bases_mask_is_valid(bases_mask):
                    raise Exception("Invalid bases mask: '%s'" %
                                    bases_mask)
        # Update for variants of standard protocol
        if protocol == 'icell8':
            # ICELL8 single-cell RNA-seq
            # Update bcl2fastq settings appropriately
            print "Updating read trimming and masking for ICELL8"
            minimum_trimmed_read_length = 21
            mask_short_adapter_reads = 0
            # Reset the default bases mask
            bases_mask = IlluminaData.IlluminaRunInfo(
                illumina_run.runinfo_xml).bases_mask
            bases_mask = get_bases_mask_icell8(bases_mask,
                                               sample_sheet=sample_sheet)
            if not bases_mask_is_valid(bases_mask):
                raise Exception("Invalid bases mask: '%s'" %
                                bases_mask)
            # Switch to standard protocol
            protocol = 'standard'
        # Do fastq generation according to protocol
        if protocol == 'standard':
            # Standard protocol
            try:
                exit_code = bcl_to_fastq(
                    ap,
                    unaligned_dir=ap.params.unaligned_dir,
                    sample_sheet=sample_sheet,
                    primary_data_dir=primary_data_dir,
                    require_bcl2fastq=require_bcl2fastq_version,
                    bases_mask=bases_mask,
                    ignore_missing_bcl=ignore_missing_bcl,
                    ignore_missing_stats=ignore_missing_stats,
                    no_lane_splitting=no_lane_splitting,
                    minimum_trimmed_read_length=minimum_trimmed_read_length,
                    mask_short_adapter_reads=mask_short_adapter_reads,
                    create_fastq_for_index_reads=create_fastq_for_index_reads,
                    nprocessors=nprocessors,
                    runner=runner)
            except Exception as ex:
                raise Exception("Bcl2fastq stage failed: '%s'" % ex)
        elif protocol == 'icell8_atac':
            # ICELL8 single-cell ATAC-seq
            # Check for well list
            if icell8_well_list is None:
                raise Exception("Need to provide an ICELL8 well list "
                                "file")
            # Reset the default bases mask
            if bases_mask == "auto":
                bases_mask = get_bases_mask_icell8_atac(
                    illumina_run.runinfo_xml)
            print "Bases mask for ICELL8 ATAC: %s" % bases_mask
            if not bases_mask_is_valid(bases_mask):
                raise Exception("Invalid bases mask: '%s'" %
                                bases_mask)
            # Perform bcl to fastq conversion and demultiplexing
            try:
                exit_code = bcl_to_fastq_icell8_atac(
                    ap,unaligned_dir=ap.params.unaligned_dir,
                    sample_sheet=sample_sheet,
                    well_list=icell8_well_list,
                    primary_data_dir=primary_data_dir,
                    bases_mask=bases_mask,
                    ignore_missing_bcl=ignore_missing_bcl,
                    ignore_missing_stats=ignore_missing_stats,
                    no_lane_splitting=no_lane_splitting,
                    minimum_trimmed_read_length=minimum_trimmed_read_length,
                    mask_short_adapter_reads=mask_short_adapter_reads,
                    swap_i1_and_i2=icell8_swap_i1_and_i2,
                    reverse_complement=icell8_reverse_complement,
                    nprocessors=nprocessors,
                    runner=runner)
            except Exception as ex:
                raise Exception("ICELL8 scATAC-seq Fastq generation failed: "
                                "'%s'" % ex)
            # Turn off barcode analysis
            analyse_barcodes = False
        elif protocol == '10x_chromium_sc':
            # 10xGenomics Chromium SC
            if bases_mask == 'auto':
                bases_mask = None
            try:
                # Check we have cellranger
                cellranger = find_program('cellranger')
                if not cellranger:
                    raise Exception("No cellranger package found")
                cellranger_software_info = cellranger_info(cellranger)
                print "Using cellranger %s: %s" % \
                    (cellranger_software_info[-1],
                     cellranger)
                # Check we have bcl2fastq
                bcl2fastq = find_program('bcl2fastq')
                if not bcl2fastq:
                    raise Exception("No bcl2fastq package found")
                bcl2fastq = available_bcl2fastq_versions(
                    paths=(os.path.dirname(bcl2fastq),),
                    reqs='>=2.17')
                if not bcl2fastq:
                    raise Exception("No appropriate bcl2fastq software "
                                    "located")
                bcl2fastq = bcl2fastq[0]
                bcl2fastq_info = bcl_to_fastq_info(bcl2fastq)
                print "Using bcl2fastq %s: %s" % (bcl2fastq_info[-1],
                                                  bcl2fastq)
                # Store info on bcl2fastq package
                ap.metadata['bcl2fastq_software'] = bcl2fastq_info
                # Store info on cellranger package
                ap.metadata['cellranger_software'] = cellranger_software_info
                # Put a copy of sample sheet in the log directory
                shutil.copy(sample_sheet,ap.log_dir)
                # Determine output directory absolute path
                output_dir = ap.params.unaligned_dir
                if not os.path.isabs(output_dir):
                    output_dir = os.path.join(ap.analysis_dir,
                                              output_dir)
                # Run cellranger mkfastq
                exit_code = run_cellranger_mkfastq(
                    sample_sheet=sample_sheet,
                    primary_data_dir=primary_data_dir,
                    output_dir=output_dir,
                    lanes=(None if lanes is None
                           else ','.join([str(l) for l in lanes])),
                    bases_mask=bases_mask,
                    cellranger_exe=cellranger,
                    cellranger_jobmode=cellranger_jobmode,
                    cellranger_maxjobs=cellranger_maxjobs,
                    cellranger_mempercore=cellranger_mempercore,
                    cellranger_jobinterval=cellranger_jobinterval,
                    cellranger_localcores=cellranger_localcores,
                    cellranger_localmem=cellranger_localmem,
                    working_dir=ap.analysis_dir,
                    log_dir=ap.log_dir)
            except Exception as ex:
                raise Exception("'cellranger mkfastq' stage failed: "
                                "'%s'" % ex)
            # Turn off barcode analysis
            analyse_barcodes = False
        elif protocol == '10x_chromium_sc_atac':
            # 10xGenomics Chromium scATAC-seq
            exit_code = bcl_to_fastq_10x_chromium_sc_atac(
                ap,
                output_dir=ap.params.unaligned_dir,
                sample_sheet=sample_sheet,
                primary_data_dir=primary_data_dir,
                lanes=lanes,
                bases_mask=bases_mask,
                cellranger_jobmode=cellranger_jobmode,
                cellranger_maxjobs=cellranger_maxjobs,
                cellranger_mempercore=cellranger_mempercore,
                cellranger_jobinterval=cellranger_jobinterval,
                cellranger_localcores=cellranger_localcores,
                cellranger_localmem=cellranger_localmem,
                log_dir=ap.log_dir
            )
            # Turn off barcode analysis
            analyse_barcodes = False
        else:
            # Unknown protocol
            raise Exception("Unknown protocol '%s'" % protocol)
        # Check the outputs
        if exit_code != 0:
            raise Exception("Fastq generation finished with error: "
                            "exit code %d" % exit_code)
        if not verify_fastq_generation(
                ap,
                lanes=lanes,
                icell8_well_list=icell8_well_list,
                include_sample_dir=verify_include_sample_dir):
            # Check failed
            logger.error("Failed to verify output Fastqs against "
                         "sample sheet")
            # Try to load the data from unaligned dir
            try:
                illumina_data = IlluminaData.IlluminaData(
                    ap.analysis_dir,
                    unaligned_dir=ap.params.unaligned_dir)
            except IlluminaData.IlluminaDataError as ex:
                raise Exception("Unable to load data from %s: %s"
                                % (ap.params.unaligned_dir,ex))
            # Generate a list of missing Fastqs
            missing_fastqs = IlluminaData.list_missing_fastqs(
                illumina_data,
                sample_sheet,
                include_sample_dir=verify_include_sample_dir)
            assert(len(missing_fastqs) > 0)
            missing_fastqs_file = os.path.join(ap.log_dir,
                                               "missing_fastqs.log")
            print "Writing list of missing Fastq files to %s" % \
                missing_fastqs_file
            with open(missing_fastqs_file,'w') as fp:
                for fq in missing_fastqs:
                    fp.write("%s\n" % fq)
            # Create empty FASTQs
            if create_empty_fastqs is None:
                try:
                    create_empty_fastqs = \
                        ap.settings.platform[ap.metadata.platform].\
                        create_empty_fastqs
                except (KeyError,AttributeError):
                    pass
            if create_empty_fastqs is None:
                create_empty_fastqs = \
                    ap.settings.bcl2fastq.create_empty_fastqs
            if create_empty_fastqs:
                logger.warning("Making 'empty' placeholder Fastqs")
                for fq in missing_fastqs:
                    fastq = os.path.join(ap.analysis_dir,
                                         ap.params.unaligned_dir,fq)
                    print "-- %s" % fastq
                    if not os.path.exists(os.path.dirname(fastq)):
                        mkdirs(os.path.dirname(fastq))
                    with gzip.GzipFile(filename=fastq,mode='wb') as fp:
                        fp.write('')
            else:
                raise Exception("Fastq generation failed to produce "
                                "expected outputs")
    # Generate statistics
    if generate_stats:
        fastq_statistics(ap,
                         stats_file=stats_file,
                         per_lane_stats_file=per_lane_stats_file,
                         unaligned_dir=ap.params.unaligned_dir,
                         nprocessors=nprocessors,
                         runner=runner)
    # Run barcode analysis
    if analyse_barcodes:
        # Determine output directory
        if barcode_analysis_dir is not None:
            ap.params['barcode_analysis_dir'] = barcode_analysis_dir
        elif ap.params.barcode_analysis_dir is None:
            ap.params['barcode_analysis_dir'] = 'barcode_analysis'
        barcode_analysis_dir = ap.params.barcode_analysis_dir
        if not os.path.isabs(barcode_analysis_dir):
            barcode_analysis_dir = os.path.join(ap.params.analysis_dir,
                                                barcode_analysis_dir)
        # Report title
        title = "Barcode analysis for %s" % ap.metadata.run_name
        # Log file
        log_file = os.path.join(ap.log_dir,"analyse_barcodes.log")
        # Set up runner
        if runner is None:
            runner = ap.settings.general.default_runner
        runner.set_log_dir(ap.log_dir)
        # Get scheduler parameters
        max_jobs = ap.settings.general.max_concurrent_jobs
        poll_interval = ap.settings.general.poll_interval
        # Create and run barcode analysis pipeline
        barcode_analysis = AnalyseBarcodes(
            os.path.join(
                ap.params.analysis_dir,
                ap.params.unaligned_dir))
        barcode_analysis.run(
            barcode_analysis_dir,
            title=title,
            lanes=lanes,
            sample_sheet=sample_sheet,
            log_file=log_file,
            runner=runner,
            max_jobs=max_jobs,
            poll_interval=poll_interval,
            verbose=False)
    # Make a 'projects.info' metadata file
    if lanes:
        ap.update_project_metadata_file()
    else:
        ap.make_project_metadata_file()
    # Remove primary data
    if remove_primary_data:
        remove_primary_data(ap)

# TODO: remove dependency on AutoProcessor instance for
# TODO: runner, source and target destinations etc
def get_primary_data(ap,runner=None):
    """
    Acquire the primary sequencing data (i.e. BCL files)

    Copies the primary sequencing data (bcl files etc) to a local area
    using rsync.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      runner (JobRunner): (optional) specify a non-default job runner
        to use for primary data rsync
    """
    # Source and target directories
    data_dir = ap.params.data_dir
    ap.params["primary_data_dir"] = ap.add_directory('primary_data')
    # Set up runner
    if runner is None:
        runner = ap.settings.runners.rsync
    runner.set_log_dir(ap.log_dir)
    # Run rsync command
    rsync_cmd = general_apps.rsync(data_dir,
                                   ap.params.primary_data_dir,
                                   prune_empty_dirs=True,
                                   extra_options=('--copy-links',
                                                  '--include=*/',
                                                  '--include=Data/**',
                                                  '--include=RunInfo.xml',
                                                  '--include=SampleSheet.csv',
                                                  '--include=RTAComplete.txt',
                                                  '--include=runParameters.xml',
                                                  '--include=RunParameters.xml',
                                                  '--exclude=*'))
    print "Running %s" % rsync_cmd
    rsync_job = SchedulerJob(runner,
                             rsync_cmd.command_line,
                             name='rsync.primary_data',
                             working_dir=os.getcwd())
    rsync_job.start()
    try:
        rsync_job.wait(
            poll_interval=ap.settings.general.poll_interval
        )
    except KeyboardInterrupt:
        logger.warning("Keyboard interrupt, terminating primary data "
                       "rsync operation")
        rsync_job.terminate()
        return -1
    exit_code = rsync_job.exit_code
    print "rsync of primary data completed: exit code %s" % exit_code
    if exit_code != 0:
        logger.error("Failed to acquire primary data (non-zero "
                     "exit code returned)")
    return exit_code

#TODO: remove dependency on AutoProcessor instance
#TODO: and use 'primary_data_dir' parameter instead of generating
#TODO: location on the fly
def remove_primary_data(ap):
    """Remove primary data

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to remove data from
    """
    primary_data = os.path.join(ap.params.primary_data_dir,
                                os.path.basename(ap.params.data_dir))
    if os.path.isdir(primary_data):
        print "Removing copy of primary data in %s" % primary_data
        shutil.rmtree(primary_data)

def bcl_to_fastq(ap,unaligned_dir,sample_sheet,primary_data_dir,
                 require_bcl2fastq=None,bases_mask=None,
                 ignore_missing_bcl=False,ignore_missing_stats=False,
                 no_lane_splitting=None,minimum_trimmed_read_length=None,
                 mask_short_adapter_reads=None,
                 create_fastq_for_index_reads=False,
                 nprocessors=None,runner=None):
    """
    Generate FASTQ files from the raw BCL files

    Performs FASTQ generation from raw BCL files produced by an Illumina
    sequencer, by running the external 'bclToFastq.py' program (which
    wraps the 'configureBclToFastq' and 'make' steps).

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      unaligned_dir (str): output directory for bcl-to-fastq conversion
      sample_sheet (str): path to input sample sheet file
      primary_data_dir (str): path to the top-level directory holding
        the sequencing data
      require_bcl2fastq (str): if set then should be a string of the form
        '1.8.4' or '>2.0' explicitly specifying the version of
        bcl2fastq to use. (Default to use specifications from the
        settings)
      bases_mask (str): if set then use this as an alternative bases mask
        setting
      ignore_missing_bcl (bool): if True then run bcl2fastq with
        --ignore-missing-bcl
      ignore_missing_stats (bool): if True then run bcl2fastq with
        --ignore-missing-stats
      no_lane_splitting (bool): if True then run bcl2fastq with
        --no-lane-splitting
      minimum_trimmed_read_length (int): if set then supply to bcl2fastq
        via --minimum-trimmed-read-length
      mask_short_adapter_reads (int): if set then supply to bcl2fastq via
        --mask-short-adapter-reads
      create_fastq_for_index_reads (boolean): if True then also create
        Fastq files for index reads (default, don't create index read
        Fastqs)
      nprocessors (int): number of processors to run bclToFastq.py with
      runner (JobRunner): (optional) specify a non-default job runner to
        use for fastq generation
    """
    # Directories
    analysis_dir = ap.params.analysis_dir
    # Bases mask
    if bases_mask is None:
        bases_mask = ap.params.bases_mask
    # Check for basic information needed to do bcl2fastq conversion
    if ap.params.data_dir is None:
        raise Exception("No source data directory")
    if bases_mask is None:
        raise Exception("No bases mask")
    # Number of cores
    if nprocessors is None:
        nprocessors = ap.settings.bcl2fastq.nprocessors
    # Whether to use lane splitting
    if no_lane_splitting is None:
        try:
            no_lane_splitting = ap.settings.platform[ap.metadata.platform].no_lane_splitting
        except (KeyError,AttributeError):
            pass
        if no_lane_splitting is None:
            no_lane_splitting = ap.settings.bcl2fastq.no_lane_splitting
    # Determine which bcl2fastq software to use
    if require_bcl2fastq is None:
        try:
            require_bcl2fastq = ap.settings.platform[ap.metadata.platform].bcl2fastq
        except (KeyError,AttributeError):
            pass
        if require_bcl2fastq is None:
            require_bcl2fastq = ap.settings.bcl2fastq.default_version
    if require_bcl2fastq is not None:
        print "Platform '%s' requires bcl2fastq version %s" \
            % (ap.metadata.platform,require_bcl2fastq)
    else:
        logger.warning("No bcl2fastq version explicitly specified")
    bcl2fastq = available_bcl2fastq_versions(require_bcl2fastq)
    if bcl2fastq:
        bcl2fastq_exe = bcl2fastq[0]
        bcl2fastq_info = bcl_to_fastq_info(bcl2fastq_exe)
    else:
        raise Exception("No appropriate bcl2fastq software located")
    # Store info on bcl2fastq package
    ap.metadata['bcl2fastq_software'] = bcl2fastq_info
    # Generate temporary sample sheet with required format
    fmt = get_required_samplesheet_format(bcl2fastq_info[2])
    tmp_sample_sheet = os.path.join(ap.tmp_dir,
                                    "SampleSheet.%s.%s.csv" %
                                    (fmt,
                                     time.strftime("%Y%m%d%H%M%S")))
    print "Generating '%s' format sample sheet: %s" % (fmt,tmp_sample_sheet)
    make_custom_sample_sheet(sample_sheet,
                             tmp_sample_sheet,
                             fmt=fmt)
    # Put a copy in the log directory
    shutil.copy(tmp_sample_sheet,ap.log_dir)
    # Create bcl2fastq directory
    bcl2fastq_dir = ap.add_directory(unaligned_dir)
    # Determine initial number of mismatches
    nmismatches = get_nmismatches(bases_mask)
    # Check for barcode collisions
    collisions = check_barcode_collisions(tmp_sample_sheet,
                                          nmismatches)
    if collisions:
        # Report problem barcodes
        logger.warning("Barcode collisions detected using %d mismatches"
                       % nmismatches)
        for collision in collisions:
            logger.warning("Barcode collision for barcodes: %s, %s"
                           % (collision[0],collision[1]))
        # Reduce mismatches to try and address the collisions
        logger.warning("Attempting to address by adjusting #mismatches")
        while nmismatches > 0 and collisions:
            nmismatches -= 1
            collisions = check_barcode_collisions(tmp_sample_sheet,
                                                  nmismatches)
            if not collisions:
                print "No collisions using %d mismatches" % nmismatches
                break
        else:
            # Unable to address collisions, bail out
            raise Exception("Barcode collisions with zero mismatches "
                            "(duplicated indexes?): unable to proceed")
    else:
        print "No barcode collisions detected using %d mismatches" % \
            nmismatches
    # Report values and settings
    print "Bcl-to-fastq exe      : %s" % bcl2fastq_exe
    print "Bcl-to-fastq version  : %s %s" % (bcl2fastq_info[1],
                                             bcl2fastq_info[2])
    print "Sample sheet          : %s" % os.path.basename(tmp_sample_sheet)
    print "Bases mask            : %s" % bases_mask
    print "Nmismatches           : %d" % nmismatches
    print "Nprocessors           : %s" % nprocessors
    print "Ignore missing bcl    : %s" % ignore_missing_bcl
    print "Ignore missing stats  : %s" % ignore_missing_stats
    print "No lane splitting     : %s" % no_lane_splitting
    print "Min trimmed read len  : %s" % minimum_trimmed_read_length
    print "Mask short adptr reads: %s" % mask_short_adapter_reads
    print "Create index Fastqs   : %s" % create_fastq_for_index_reads
    # Set up runner
    if runner is None:
        runner = ap.settings.runners.bcl2fastq
    runner.set_log_dir(ap.log_dir)
    # Run bcl2fastq
    bcl2fastq = Command('bclToFastq.py',
                        '--nprocessors',nprocessors,
                        '--use-bases-mask',bases_mask,
                        '--nmismatches',nmismatches,
                        '--ignore-missing-control')
    if ignore_missing_bcl:
        bcl2fastq.add_args('--ignore-missing-bcl')
    if ignore_missing_stats:
        bcl2fastq.add_args('--ignore-missing-stats')
    if no_lane_splitting:
        bcl2fastq.add_args('--no-lane-splitting')
    if minimum_trimmed_read_length is not None:
        bcl2fastq.add_args('--minimum-trimmed-read-length',
                           minimum_trimmed_read_length)
    if mask_short_adapter_reads is not None:
        bcl2fastq.add_args('--mask-short-adapter-reads',
                           mask_short_adapter_reads)
    if create_fastq_for_index_reads:
        bcl2fastq.add_args('--create-fastq-for-index-reads')
    bcl2fastq.add_args('--platform',
                       ap.metadata.platform,
                       '--bcl2fastq_path',
                       bcl2fastq_exe,
                       primary_data_dir,
                       bcl2fastq_dir,
                       tmp_sample_sheet)

    print "Running %s" % bcl2fastq
    bcl2fastq_job = SchedulerJob(runner,
                                 bcl2fastq.command_line,
                                 name='bclToFastq',
                                 working_dir=os.getcwd())
    bcl2fastq_job.start()
    try:
        bcl2fastq_job.wait(
            poll_interval=ap.settings.general.poll_interval
        )
    except KeyboardInterrupt,ex:
        logger.warning("Keyboard interrupt, terminating bcl2fastq")
        bcl2fastq_job.terminate()
        raise ex
    exit_code = bcl2fastq_job.exit_code
    print "bcl2fastq completed: exit code %s" % exit_code
    if exit_code != 0:
        logger.error("bcl2fastq exited with an error")
    return exit_code

def bcl_to_fastq_10x_chromium_sc_atac(ap,output_dir,sample_sheet,
                                      primary_data_dir,lanes=None,
                                      bases_mask=None,
                                      cellranger_jobmode=None,
                                      cellranger_maxjobs=None,
                                      cellranger_mempercore=None,
                                      cellranger_jobinterval=None,
                                      cellranger_localcores=None,
                                      cellranger_localmem=None,
                                      log_dir=None):
    """
    Generate FASTQ files for 10xGenomics single-cell ATAC-seq run

    Performs FASTQ generation from raw BCL files produced by an
    Illumina sequencer using the 10xGenomics Chromium single-cell
    (sc) ATAC-seq protocol, by running 'cellranger-atac mkfastq'.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      output_dir (str): output directory for bcl-to-fastq conversion
      sample_sheet (str): path to input sample sheet file
      primary_data_dir (str): path to the top-level directory holding
        the sequencing data
      bases_mask (str): if set then use this as an alternative bases
        mask setting (default is to acquire from the autoprocessor
        parameters)
      ...TBD...
    """
    # Load run data
    illumina_run = IlluminaData.IlluminaRun(primary_data_dir,
                                            platform=ap.metadata.platform)
    # Deal with bases mask
    if bases_mask is None:
        bases_mask = ap.params.bases_mask
    if bases_mask == 'auto':
        # Update bases mask to only use first 8 bases from
        # first index e.g. I8nnnnnnnn and convert second index
        # to read e.g. Y16
        print "Determining bases mask from RunInfo.xml"
        bases_mask = get_bases_mask_10x_atac(illumina_run.runinfo_xml)
        print "Bases mask: %s (updated for 10x scATAC-seq)" % bases_mask
        if not bases_mask_is_valid(bases_mask):
            raise Exception("Invalid bases mask: '%s'" %
                            bases_mask)
    # Check we have cellranger-atac
    cellranger_atac = find_program('cellranger-atac')
    if not cellranger_atac:
        raise Exception("No cellranger package found")
    cellranger_package_info = cellranger_info(cellranger_atac)
    print "Using cellranger-atac %s: %s" % \
        (cellranger_package_info[-1],
         cellranger_atac)
    # Check we have bcl2fastq
    bcl2fastq = find_program('bcl2fastq')
    if not bcl2fastq:
        raise Exception("No bcl2fastq package found")
    bcl2fastq = available_bcl2fastq_versions(
        paths=(os.path.dirname(bcl2fastq),),
        reqs='>=2.17')
    if not bcl2fastq:
        raise Exception("No appropriate bcl2fastq software "
                        "located")
    bcl2fastq = bcl2fastq[0]
    bcl2fastq_info = bcl_to_fastq_info(bcl2fastq)
    print "Using bcl2fastq %s: %s" % (bcl2fastq_info[-1],
                                      bcl2fastq)
    # Store info on bcl2fastq package
    ap.metadata['bcl2fastq_software'] = bcl2fastq_info
    # Store info on cellranger package
    ap.metadata['cellranger_software'] = cellranger_package_info
    # Put a copy of sample sheet in the log directory
    shutil.copy(sample_sheet,log_dir)
    # Determine output directory absolute path
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(ap.analysis_dir,
                                  output_dir)
    # Working directory (set to analysis dir)
    working_dir = ap.analysis_dir
    # Report values and settings
    print "Cellranger-atac exe    : %s" % cellranger_atac
    print "Cellranger-atac version: %s %s" % (cellranger_package_info[1],
                                              cellranger_package_info[2])
    print "Bcl-to-fastq exe       : %s" % bcl2fastq
    print "Bcl-to-fastq version   : %s %s" % (bcl2fastq_info[1],
                                              bcl2fastq_info[2])
    print "Sample sheet           : %s" % os.path.basename(sample_sheet)
    print "Bases mask             : %s" % bases_mask
    print "Cellranger jobmode     : %s" %cellranger_jobmode
    print "Cellranger maxjobs     : %s" % cellranger_maxjobs
    print "Cellranger mempercore  : %s" % cellranger_mempercore
    print "Cellranger jobinterval : %s" % cellranger_jobinterval
    print "Cellranger localcores  : %s" % cellranger_localcores
    print "Cellranger localmem    : %s" % cellranger_localmem
    print "Working directory      : %s" % working_dir
    print "Log directory          : %s" % log_dir
    # Run cellranger-atac mkfastq
    try:
        return run_cellranger_mkfastq(
            sample_sheet=sample_sheet,
            primary_data_dir=primary_data_dir,
            output_dir=output_dir,
            lanes=(None if lanes is None
                   else ','.join([str(l) for l in lanes])),
            bases_mask=bases_mask,
            cellranger_exe=cellranger_atac,
            cellranger_jobmode=cellranger_jobmode,
            cellranger_maxjobs=cellranger_maxjobs,
            cellranger_mempercore=cellranger_mempercore,
            cellranger_jobinterval=cellranger_jobinterval,
            cellranger_localcores=cellranger_localcores,
            cellranger_localmem=cellranger_localmem,
            working_dir=working_dir,
            log_dir=log_dir)
    except Exception as ex:
        raise Exception("'cellranger-atac mkfastq' failed: "
                        "'%s'" % ex)

def bcl_to_fastq_icell8_atac(ap,unaligned_dir,sample_sheet,
                             well_list,primary_data_dir,
                             bases_mask=None,
                             ignore_missing_bcl=False,
                             ignore_missing_stats=False,
                             no_lane_splitting=None,
                             minimum_trimmed_read_length=None,
                             mask_short_adapter_reads=None,
                             swap_i1_and_i2=False,
                             reverse_complement=None,
                             nprocessors=None,runner=None):
    """
    Generate FASTQ files for ICELL8 scATAC-seq data

    Performs FASTQ generation from raw BCL files produced by an Illumina
    sequencer from ICELL8 scATAC-seq samples.

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      unaligned_dir (str): output directory for bcl-to-fastq conversion
      sample_sheet (str): path to input sample sheet file
      well_list (str): path to the ICELL8 well list file to use
      primary_data_dir (str): path to the top-level directory holding
        the sequencing data
      bases_mask (str): if set then use this as an alternative bases mask
        setting
      ignore_missing_bcl (bool): if True then run bcl2fastq with
        --ignore-missing-bcl
      ignore_missing_stats (bool): if True then run bcl2fastq with
        --ignore-missing-stats
      no_lane_splitting (bool): if True then run bcl2fastq with
        --no-lane-splitting
      minimum_trimmed_read_length (int): if set then supply to bcl2fastq
        via --minimum-trimmed-read-length
      mask_short_adapter_reads (int): if set then supply to bcl2fastq via
        --mask-short-adapter-reads
      swap_i1_and_i2 (bool): if True then swap I1 and I2 reads when
        matching to barcodes in the ICELL8 well list
      reverse_complement (str): one of 'i1', 'i2', 'both', or None; if
        set then the specified index reads will be reverse complemented
        when matching to barcodes in the ICELL8 well list
      nprocessors (int): number of processors to run bclToFastq.py with
      runner (JobRunner): (optional) specify a non-default job runner to
        use for fastq generation
    """
    # Set up runner
    if runner is None:
        runner = ap.settings.runners.bcl2fastq
    runner.set_log_dir(ap.log_dir)
    # Number of cores
    if nprocessors is None:
        nprocessors = ap.settings.bcl2fastq.nprocessors
    # Do Fastq generation to a temporary directory
    icell8_tmp = os.path.join(ap.analysis_dir,"__icell8")
    if os.path.exists(icell8_tmp):
        # Remove existing working directory
        logger.warning("Found existing working directory '%s '"
                       "(will be removed)" % icell8_tmp)
        shutil.rmtree(icell8_tmp)
    mkdirs(icell8_tmp)
    icell8_unaligned = os.path.join("__icell8","bcl2fastq")
    tmp_bcl2fastq = os.path.join(ap.analysis_dir,icell8_unaligned)
    if not os.path.exists(tmp_bcl2fastq):
        exit_code = bcl_to_fastq(
            ap,
            unaligned_dir=tmp_bcl2fastq,
            sample_sheet=sample_sheet,
            primary_data_dir=primary_data_dir,
            require_bcl2fastq=">=2",
            bases_mask=bases_mask,
            ignore_missing_bcl=ignore_missing_bcl,
            ignore_missing_stats=ignore_missing_stats,
            no_lane_splitting=no_lane_splitting,
            minimum_trimmed_read_length=minimum_trimmed_read_length,
            mask_short_adapter_reads=mask_short_adapter_reads,
            create_fastq_for_index_reads=True,
            nprocessors=nprocessors,
            runner=runner)
    # Load data from the temporary directory and get Fastqs
    illumina_data = IlluminaData.IlluminaData(
            ap.analysis_dir,
            unaligned_dir=icell8_unaligned)
    fastqs = []
    for project in illumina_data.projects:
        for sample in project.samples:
            for fq in sample.fastq:
                fastqs.append(os.path.join(sample.dirn,fq))
    if not fastqs:
        logger.error("No Fastqs were produced, cannot proceed")
        return 1
    # Do demultiplexing into samples based on well list
    tmp_demultiplex_dir = os.path.join(icell8_tmp,"demultiplexed")
    if not os.path.exists(tmp_demultiplex_dir):
        demultiplexer = Command('demultiplex_icell8_atac.py',
                                '--mode=samples',
                                '--output-dir',tmp_demultiplex_dir,
                                '-n',nprocessors,
                                '--update-read-headers')
        if swap_i1_and_i2:
            demultiplexer.add_args('--swap-i1-i2')
        if reverse_complement:
            demultiplexer.add_args('--reverse-complement=%s' %
                                   reverse_complement)
        demultiplexer.add_args(well_list)
        demultiplexer.add_args(*fastqs)
        print "Running %s" % demultiplexer
        demultiplexer_job = SchedulerJob(runner,
                                         demultiplexer.command_line,
                                         name='demultiplex_icell8',
                                         working_dir=os.getcwd())
        demultiplexer_job.start()
        try:
            demultiplexer_job.wait(
                poll_interval=ap.settings.general.poll_interval
            )
        except KeyboardInterrupt,ex:
            logger.warning("Keyboard interrupt, terminating demultiplexer")
            demultiplexer_job.terminate()
            raise ex
        exit_code = demultiplexer_job.exit_code
        print "demultiplexer completed: exit code %s" % exit_code
        if exit_code != 0:
            logger.error("demultiplexer exited with an error")
            return exit_code
    # Build bcl2fastq2-style output dir
    bcl2fastq_dir = os.path.join(ap.tmp_dir,"bcl2fastq")
    print "Building temporary output dir: %s" % bcl2fastq_dir
    for d in ('Stats','Reports',):
        mkdirs(os.path.join(bcl2fastq_dir,d))
    project = illumina_data.projects[0]
    mkdirs(os.path.join(bcl2fastq_dir,project.name))
    # Copy the reports and JSON file
    for f in ('icell8_atac_stats.xlsx',
              'icell8_atac_stats.json'):
        os.link(os.path.join(tmp_demultiplex_dir,f),
                os.path.join(bcl2fastq_dir,'Reports',f))
    # Copy fastqs to final location
    for f in os.listdir(tmp_demultiplex_dir):
        if f.endswith(".fastq.gz"):
            # Undetermined goes to top level
            if f.startswith("Undetermined_S0_"):
                os.link(os.path.join(tmp_demultiplex_dir,f),
                        os.path.join(bcl2fastq_dir,f))
            else:
                os.link(os.path.join(tmp_demultiplex_dir,f),
                        os.path.join(bcl2fastq_dir,project.name,f))
    print "Moving %s to final destination: %s" % (bcl2fastq_dir,
                                                  unaligned_dir)
    os.rename(bcl2fastq_dir,unaligned_dir)
    # Remove the intermediate directories
    shutil.rmtree(icell8_tmp)
    # Finish
    return 0

def verify_fastq_generation(ap,unaligned_dir=None,lanes=None,
                            icell8_well_list=None,
                            include_sample_dir=False):
    """Check that generated Fastqs match sample sheet predictions

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to do Fastqs verification on
      unaligned_dir (str): explicitly specify the bcl2fastq output
        directory to check
      lanes (list): specify a list of lane numbers (integers) to
        check (others will be ignored)
      icell8_well_list (str): path to the ICELL8 well list file to
        verify against, instead of samplesheet (for ICELL8 data)
      include_sample_dir (bool): if True then include a
        'sample_name' directory level when checking for
        bcl2fastq2 outputs, even if one shouldn't be present

     Returns:
       True if outputs match sample sheet, False otherwise.
    """
    if unaligned_dir is None:
        if ap.params.unaligned_dir is not None:
            unaligned_dir = ap.params.unaligned_dir
        else:
            raise Exception("Bcl2fastq output directory not defined")
    print "Checking bcl2fastq output directory '%s'" % unaligned_dir
    bcl_to_fastq_dir = os.path.join(ap.analysis_dir,unaligned_dir)
    if not os.path.isdir(bcl_to_fastq_dir):
        # Directory doesn't exist
        return False
    # Make a temporary sample sheet to verify against
    tmp_sample_sheet = os.path.join(ap.tmp_dir,
                                    "SampleSheet.verify.%s.csv" %
                                    time.strftime("%Y%m%d%H%M%S"))
    make_custom_sample_sheet(ap.params.sample_sheet,
                             tmp_sample_sheet,
                             lanes=lanes)
    # Verify against ICELL8 well list file
    if icell8_well_list is not None:
        # Get the project name from the sample sheet
        sample_sheet = IlluminaData.SampleSheet(tmp_sample_sheet)
        project = sample_sheet.data[0][sample_sheet.sample_project_column]
        # Check Fastqs for each sample
        missing_fastqs = list()
        for ii,sample in enumerate(ICell8WellList(icell8_well_list).samples(),
                                   start=1):
            for read in ('R1','R2','I1','I2'):
                fq = os.path.join(bcl_to_fastq_dir,
                                  project,
                                  "%s_S%d_%s_001.fastq.gz" % (sample,
                                                              ii,
                                                              read))
                if not os.path.exists(fq):
                    logger.warning("Missing: %s" % fq)
                    missing_fastqs.append(fq)
        # Return status of verification
        return (len(missing_fastqs) == 0)
    # Try to create an IlluminaData object
    try:
        illumina_data = IlluminaData.IlluminaData(
            ap.analysis_dir,
            unaligned_dir=unaligned_dir)
    except IlluminaData.IlluminaDataError as ex:
        # Failed to initialise
        logger.warning("Failed to get information from %s: %s" %
                        (bcl_to_fastq_dir,ex))
        return False
    # Do check
    return IlluminaData.verify_run_against_sample_sheet(
        illumina_data,
        tmp_sample_sheet,
        include_sample_dir=include_sample_dir)

def fastq_statistics(ap,stats_file=None,per_lane_stats_file=None,
                     unaligned_dir=None,sample_sheet=None,add_data=False,
                     nprocessors=None,runner=None):
    """Generate statistics for Fastq files

    Generates statistics for all Fastq files found in the
    'unaligned' directory, by running the 'fastq_statistics.py'
    program.

    Arguments
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to create Fastqs for
      stats_file (str): path of a non-default file to write the
        statistics to (defaults to 'statistics.info' unless
        over-ridden by local settings)
      per_lane_stats_file (str): path for per-lane statistics
        output file (defaults to 'per_lane_statistics.info'
        unless over-ridden by local settings)
      unaligned_dir (str): output directory for bcl-to-fastq
        conversion
      sample_sheet (str): path to sample sheet file used in
        bcl-to-fastq conversion
      add_data (bool): if True then add stats to the existing
        stats files (default is to overwrite existing stats
        files)
      nprocessors (int): number of cores to use when running
        'fastq_statistics.py'
      runner (JobRunner): (optional) specify a non-default job
        runner to use for running 'fastq_statistics.py'
    """
    # Get file names for output files
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
    # Sort out unaligned_dir
    if unaligned_dir is None:
        if ap.params.unaligned_dir is None:
            ap.params['unaligned_dir'] = 'bcl2fastq'
        unaligned_dir = ap.params.unaligned_dir
    if not os.path.exists(os.path.join(ap.params.analysis_dir,unaligned_dir)):
        logger.error("Unaligned dir '%s' not found" % unaligned_dir)
    # Check for sample sheet
    if sample_sheet is None:
        sample_sheet = ap.params['sample_sheet']
    # Check if any Fastqs are newer than stats files
    newest_mtime = 0
    for f in (stats_file,per_lane_stats_file,):
        try:
            newest_mtime = max(newest_mtime,
                               os.path.getmtime(f))
        except OSError:
            # Missing file
            newest_mtime = 0
            break
    illumina_data = IlluminaData.IlluminaData(ap.params.analysis_dir,
                                              unaligned_dir)
    if newest_mtime > 0:
        regenerate_stats = False
        for project in illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq:
                    if (os.path.getmtime(os.path.join(sample.dirn,fq)) >
                        newest_mtime):
                        regenerate_stats = True
                        break
        if regenerate_stats:
            logger.warning("Fastqs are newer than stats files")
        else:
            # Don't rerun the stats, just regenerate the report
            logger.warning("Stats files are newer than Fastqs")
            processing_qc_html = os.path.join(ap.analysis_dir,
                                              "processing_qc.html")
            report_processing_qc(ap,processing_qc_html)
            return
    # Set up runner
    if runner is None:
        runner = ap.settings.runners.stats
    runner.set_log_dir(ap.log_dir)
    # Number of cores
    if nprocessors is None:
        nprocessors = ap.settings.fastq_stats.nprocessors
    # Generate statistics
    fastq_statistics_cmd = Command(
        'fastq_statistics.py',
        '--unaligned',unaligned_dir,
        '--sample-sheet',sample_sheet,
        '--output',os.path.join(ap.params.analysis_dir,stats_file),
        '--per-lane-stats',os.path.join(ap.params.analysis_dir,
                                        per_lane_stats_file),
        ap.params.analysis_dir,
        '--nprocessors',nprocessors
    )
    if add_data:
        fastq_statistics_cmd.add_args('--update')
    print "Generating statistics: running %s" % fastq_statistics_cmd
    fastq_statistics_job = SchedulerJob(
        runner,
        fastq_statistics_cmd.command_line,
        name='fastq_statistics',
        working_dir=ap.analysis_dir
    )
    fastq_statistics_job.start()
    try:
        fastq_statistics_job.wait(
            poll_interval=ap.settings.general.poll_interval
        )
    except KeyboardInterrupt as ex:
        logger.warning("Keyboard interrupt, terminating fastq_statistics")
        fastq_statistics_job.terminate()
        raise ex
    exit_code = fastq_statistics_job.exit_code
    print "fastq_statistics completed: exit code %s" % exit_code
    if exit_code != 0:
        raise Exception("fastq_statistics exited with an error")
    ap.params['stats_file'] = stats_file
    ap.params['per_lane_stats_file'] = per_lane_stats_file
    print "Statistics generation completed: %s" % ap.params.stats_file
    print "Generating processing QC report"
    processing_qc_html = os.path.join(ap.analysis_dir,
                                      "processing_qc.html")
    report_processing_qc(ap,processing_qc_html)
