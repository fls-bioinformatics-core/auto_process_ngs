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
from ..icell8.utils import get_icell8_bases_mask
from ..tenx_genomics_utils import has_chromium_sc_indices
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

MAKE_FASTQS_PROTOCOLS = ('standard','icell8','10x_chromium_sc')

#######################################################################
# Command functions
#######################################################################

def make_fastqs(ap,protocol='standard',platform=None,
                unaligned_dir=None,sample_sheet=None,lanes=None,
                ignore_missing_bcl=False,ignore_missing_stats=False,
                skip_rsync=False,remove_primary_data=False,
                nprocessors=None,require_bcl2fastq_version=None,
                bases_mask=None,no_lane_splitting=None,
                minimum_trimmed_read_length=None,
                mask_short_adapter_reads=None,
                generate_stats=True,stats_file=None,
                per_lane_stats_file=None,
                skip_fastq_generation=False,
                only_fetch_primary_data=False,
                create_empty_fastqs=None,runner=None,
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
      stats_file (str): if set then use this as the name of the output
        per-fastq stats file.
      per_lane_stats_file (str): if set then use this as the name of
        the output per-lane stats file.
      skip_fastq_generation (bool): if True then don't perform fastq
        generation
      only_fetch_primary_data (bool): if True then fetch primary data,
        don't do anything else
      create_empty_fastqs (bool): if True then create empty 'placeholder'
        fastq files for any missing fastqs after bcl2fastq
        (must have completed with zero exit status)
      runner (JobRunner): (optional) specify a non-default job runner
        to use for fastq generation
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
        if protocol == '10x_chromium_sc':
            # Force inclusion of sample-name subdirectories
            # when verifying Chromium SC data
            print "Sample sheet includes Chromium SC indices"
            verify_include_sample_dir = True
        else:
            # Chromium SC indices detected but not using
            # 10x_chromium_sc protocol
            raise Exception("Detected 10xGenomics Chromium SC indices "
                            "in generated sample sheet but protocol "
                            "'%s' has been specified; must use "
                            "'10x_chromium_sc' for these indices" %
                            protocol)
    # Check for pre-existing Fastq outputs
    if verify_fastq_generation(
            ap,
            unaligned_dir=ap.params.unaligned_dir,
            lanes=lanes,
            include_sample_dir=verify_include_sample_dir):
        print "Expected Fastq outputs already present"
        skip_rsync = True
        skip_fastq_generation = True
    # Check if there's anything to do
    if (skip_rsync and skip_fastq_generation) and not generate_stats:
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
    if not skip_rsync:
        if get_primary_data(ap) != 0:
            logger.error("Failed to acquire primary data")
            raise Exception, "Failed to acquire primary data"
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
        if bases_mask is None:
            bases_mask = ap.params.bases_mask
        print "Bases mask setting    : %s" % bases_mask
        if protocol != '10x_chromium_sc':
            if bases_mask == "auto":
                print "Determining bases mask from RunInfo.xml"
                bases_mask = get_bases_mask(illumina_run.runinfo_xml,
                                            sample_sheet)
                if not bases_mask_is_valid(bases_mask):
                    raise Exception("Invalid bases mask: '%s'" %
                                    bases_mask)
        ap.params.bases_mask = bases_mask
        # Do fastq generation according to protocol
        if protocol == 'icell8':
            # ICell8 data
            # Update bcl2fastq settings appropriately
            print "Updating read trimming and masking for ICell8"
            minimum_trimmed_read_length = 21
            mask_short_adapter_reads = 0
            # Reset the default bases mask
            bases_mask = IlluminaData.IlluminaRunInfo(
                illumina_run.runinfo_xml).bases_mask
            bases_mask = get_icell8_bases_mask(bases_mask,
                                               sample_sheet=sample_sheet)
            if not bases_mask_is_valid(bases_mask):
                raise Exception("Invalid bases mask: '%s'" %
                                bases_mask)
            ap.params.bases_mask = bases_mask
            # Switch to standard protocol
            protocol = 'standard'
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
                    nprocessors=nprocessors,
                    runner=runner)
            except Exception as ex:
                raise Exception("Bcl2fastq stage failed: '%s'" % ex)
        elif protocol == '10x_chromium_sc':
            # 10xGenomics Chromium SC
            if bases_mask == 'auto':
                bases_mask = None
            try:
                # Check we have cellranger
                cellranger = find_program('cellranger')
                if not cellranger:
                    raise Exception("No cellranger package found")
                print "Using cellranger %s: %s" % (
                    cellranger_info(cellranger)[-1],cellranger)
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
    if runner is not None:
        runner = fetch_runner(runner)
    else:
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
        rsync_job.wait()
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
                 mask_short_adapter_reads=None,nprocessors=None,
                 runner=None):
    """Generate FASTQ files from the raw BCL files

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
      nprocessors (int): number of processors to run bclToFastq.py with
      runner (JobRunner): (optional) specify a non-default job runner to
        use for fastq generation
    """
    # Directories
    analysis_dir = ap.params.analysis_dir
    # Bases mask
    if bases_mask is None:
        bases_mask = ap.params.bases_mask
    else:
        ap.params['bases_mask'] = bases_mask
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
    # Set up runner
    if runner is not None:
        runner = fetch_runner(runner)
    else:
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
        bcl2fastq_job.wait()
    except KeyboardInterrupt,ex:
        logger.warning("Keyboard interrupt, terminating bcl2fastq")
        bcl2fastq_job.terminate()
        raise ex
    exit_code = bcl2fastq_job.exit_code
    print "bcl2fastq completed: exit code %s" % exit_code
    if exit_code != 0:
        logger.error("bcl2fastq exited with an error")
    return exit_code

def verify_fastq_generation(ap,unaligned_dir=None,lanes=None,
                            include_sample_dir=False):
    """Check that generated Fastqs match sample sheet predictions

    Arguments:
      ap (AutoProcessor): autoprocessor pointing to the analysis
        directory to do Fastqs verification on
      unaligned_dir (str): explicitly specify the bcl2fastq output
        directory to check
      lanes (list): specify a list of lane numbers (integers) to
        check (others will be ignored)
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
            logger.debug("Bcl2fastq output directory not defined")
            return False
    else:
        logger.warning("Checking custom bcl2fastq output directory '%s'" %
                       unaligned_dir)
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
                     unaligned_dir=None,add_data=False,nprocessors=None,
                     runner=None):
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
    # Set up runner
    if runner is not None:
        runner = fetch_runner(runner)
    else:
        runner = ap.settings.runners.stats
    runner.set_log_dir(ap.log_dir)
    # Number of cores
    if nprocessors is None:
        nprocessors = ap.settings.fastq_stats.nprocessors
    # Generate statistics
    fastq_statistics_cmd = Command(
        'fastq_statistics.py',
        '--unaligned',unaligned_dir,
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
        fastq_statistics_job.wait()
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
