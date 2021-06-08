#!/usr/bin/env python
#
#     bcl2fastq.pipeline.py: pipelines for Fastq generation
#     Copyright (C) University of Manchester 2020-2021 Peter Briggs
#

"""
Pipeline components for generating Fastqs from Bcl files.

Pipeline classes:

- MakeFastqs

Pipeline task classes:

- FetchPrimaryData
- MakeSampleSheet
- GetBcl2Fastq
- GetBclConvert
- RestoreBackupDirectory
- RunBcl2Fastq
- GetBasesMaskIcell8
- GetBasesMaskIcell8Atac
- Get10xPackage
- DemultiplexIcell8Atac
- MergeFastqs
- MergeFastqDirs
- Run10xMkfastq
- FastqStatistics
- ReportProcessingQC

Utility functions:

- subset
"""

######################################################################
# Imports
######################################################################

import os
import gzip
import shutil
import tempfile
import time
from bcftbx.IlluminaData import IlluminaRun
from bcftbx.IlluminaData import IlluminaRunInfo
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import samplesheet_index_sequence
from bcftbx.IlluminaData import verify_run_against_sample_sheet
from bcftbx.IlluminaData import list_missing_fastqs
from bcftbx.utils import find_program
from bcftbx.utils import mkdirs
from bcftbx.utils import walk
from ..analysis import AnalysisFastq
from ..applications import general as general_apps
from ..applications import bcl2fastq as bcl2fastq_apps
from ..barcodes.pipeline import AnalyseBarcodes
from ..bcl2fastq.utils import available_bcl2fastq_versions
from ..bcl2fastq.utils import bases_mask_is_valid
from ..bcl2fastq.utils import bcl_to_fastq_info
from ..bcl2fastq.utils import bclconvert_info
from ..bcl2fastq.utils import check_barcode_collisions
from ..bcl2fastq.utils import get_bases_mask
from ..bcl2fastq.utils import get_nmismatches
from ..bcl2fastq.utils import get_sequencer_platform
from ..bcl2fastq.utils import make_custom_sample_sheet
from ..command import Command
from ..fastq_utils import group_fastqs_by_name
from ..fileops import Location
from ..icell8.utils import get_bases_mask_icell8
from ..icell8.utils import get_bases_mask_icell8_atac
from ..samplesheet_utils import barcode_is_10xgenomics
from ..samplesheet_utils import SampleSheetLinter
from ..pipeliner import Pipeline
from ..pipeliner import PipelineTask
from ..pipeliner import PipelineFunctionTask
from ..pipeliner import PipelineCommandWrapper
from ..pipeliner import PipelineFailure
from ..pipeliner import PipelineParam as Param
from ..pipeliner import FunctionParam
from ..pipeliner import ListParam
from ..pipeliner import PathExistsParam
from ..pipeliner import PathJoinParam
from ..pipeliner import resolve_parameter
from ..tenx_genomics_utils import add_cellranger_args
from ..tenx_genomics_utils import cellranger_info
from ..tenx_genomics_utils import spaceranger_info
from ..tenx_genomics_utils import get_bases_mask_10x_atac
from ..tenx_genomics_utils import make_qc_summary_html
from .reporting import ProcessingQCReport

# Module specific logger
import logging
logger = logging.getLogger(__name__)

######################################################################
# Constants
######################################################################

# Protocols

PROTOCOLS = ('standard',
             'mirna',
             'icell8',
             'icell8_atac',
             '10x_chromium_sc',
             '10x_atac',
             '10x_visium',
             '10x_multiome',)

# 10xGenomics protocols
PROTOCOLS_10X = ('10x_chromium_sc',
                 '10x_atac',
                 '10x_visium',
                 '10x_multiome',)

# Valid attribute names for lane subsets

LANE_SUBSET_ATTRS = (
    'protocol',
    'bases_mask',
    'trim_adapters',
    'adapter_sequence',
    'adapter_sequence2',
    'minimum_trimmed_read_length',
    'mask_short_adapter_reads',
    'create_fastq_for_index_read',
    'find_adapters_with_sliding_window',
    'no_lane_splitting',
    'icell8_well_list',
    'icell8_atac_swap_i1_and_i2',
    'icell8_atac_reverse_complement',
    'analyse_barcodes',
    'masked_index',
    'bcl_converter',
)

######################################################################
# Pipeline classes
######################################################################

class MakeFastqs(Pipeline):
    """
    Run the Fastq generation pipeline on one or more lane subsets

    Pipeline to run Fastq generation on multiple projects.

    Example usage for processing a standard run:

    >>> make_fastqs = MakeFastqs(run_dir,sample_sheet)
    >>> make_fastqs.run()

    Example for splitting a run to use different protocols for
    different lanes:

    >>> make_fastqs = MakeFastqs(run_dir,sample_sheet,
    ...                          lane_subsets=(
    ...                             subset(lanes=[1,2,3,4,5,6],
    ...                                    protocol="standard"),
    ...                             subset(lanes=[7,8],
    ...                                    protocol="10x_chromium_sc")))
    >>> make_fastqs.run()

    In this case subsets of lanes are defined by calling the
    'subset' function; each subset is processed separately using
    the protocol specified for that subset, before being merged
    into a single output directory.

    Parameters defined in the lane subsets override those defined
    globally in the pipleine.

    On completion the pipeline makes the follow outputs availble:

    - platform: the platform assigned to the primary data
    - primary_data_dir: the directory containing the primary data
    - acquired_primary_data: boolean indicating if the primary
        data exists
    - bcl2fastq_info: tuple with information on the bcl2fastq
        software used
    - cellranger_info: tuple with information on the cellranger
        software used
    - stats_file: path to the statistics file
    - stats_full: path to the full statistics file
    - per_lane_stats: path to the per-lane statistics file
    - per_lane_sample_stats: path to the per-lane per-sample
        statistics file
    - missing_fastqs: list of Fastq files that bcl2fastq failed
        to generate
    """
    def __init__(self,run_dir,sample_sheet,protocol='standard',
                 bases_mask="auto",bcl_converter='bcl2fastq',
                 platform=None,icell8_well_list=None,
                 minimum_trimmed_read_length=None,
                 mask_short_adapter_reads=None,
                 adapter_sequence=None,adapter_sequence_read2=None,
                 icell8_atac_swap_i1_and_i2=None,
                 icell8_atac_reverse_complement=None,
                 lanes=None,trim_adapters=True,fastq_statistics=True,
                 analyse_barcodes=True,lane_subsets=None):
        """
        Create a new MakeFastqs pipeline instance

        Arguments:
          run_dir (str): path to directory with raw Bcl data
            from sequencer run
          sample_sheet (str): path to sample sheet file to use
          protocol (str): default protocol to use (defaults to
            "standard")
          bases_mask (str): default bases mask to use (defaults
            to "auto")
          bcl_converter (str): default BCL-to-Fastq conversion
            software to use (defaults to "bcl2fastq")
          platform (str): optionally specify the platform for
            the sequencer run (e.g. 'miseq', 'nextseq' etc)
          icell8_well_list (str): optionally specify path to a
            well list file for ICELL8 data
          minimum_trimmed_read_length (int): optionally specify
            the minimum length of reads after adapter trimming;
            trimmed reads shorter than this length will be
            padded with Ns (bcl2fastq --minimum-trimmed-read-length
            option)
          mask_short_adapter_reads (int): optionally specify
            the length at which reads which should be completely
            masked with Ns after adapter trimming (bcl2fastq
            --mask-short-adapter-reads option)
          adapter_sequence (str): optionally specify the adapter
            sequence to use for trimming (overrides sequence set
            in the sample sheet file)
          adapter_sequence_read2 (str): optionally specify the
            'read2' adapter sequence to use for trimming
            (overrides sequence set in the sample sheet file)
          icell8_atac_swap_i1_and_i2 (bool): if True then
            swap the I1 and I2 indexes when demultiplexing ICELL8
            ATAC data
          icell8_atac_reverse_complement (str): whether to reverse
            complement I1, I2, or both, when demultiplexing ICELL8
            ATAC data
          lanes (list): if set then specifies a list of lanes to
            include; all other lanes will be excluded
          trim_adapters (bool): if True (default) then perform
            adapter trimming as part of Fastq generation
          fastq_statistics (bool): if True (default) then generate
            statistics from Fastq files
          analyse_barcodes (bool): if True (default) then perform
            barcode analysis after Fastq generation
          lane_subsets (iterable): optionally define one or more
            subsets of lanes to process separately
        """
        # Initialise the pipeline superclass
        Pipeline.__init__(self,name="MakeFastqs")

        # Inputs
        self._run_dir = run_dir
        self._sample_sheet = os.path.abspath(sample_sheet)
        self._platform = platform

        # Limit set of lanes to use?
        if lanes:
            self._use_lanes = sorted([int(l) for l in lanes])
        else:
            self._use_lanes = None
        if self._use_lanes:
            self.report("Restrict to lanes: %s" % self._use_lanes)

        # Preflight checks
        #
        # Sample sheet
        if not self._sample_sheet_is_valid(self._sample_sheet):
            raise Exception("Problems detected in sample sheet")
        # Consistent lane definitions
        if lane_subsets:
            # Check each lane is only defined once
            defined_lanes = set()
            for s in lane_subsets:
                for l in s['lanes']:
                    ll = int(l)
                    if ll in defined_lanes:
                        # Lane is defined multiple times
                        raise Exception("Lane '%s' appears multiple "
                                        "times in lane subset "
                                        "definitions" % ll)
                    else:
                        defined_lanes.add(ll)
            # Check that subsets don't include lanes
            # that are excluded elsewhere
            if self._use_lanes:
                for l in defined_lanes:
                    if l not in self._use_lanes:
                        raise Exception("Lane '%s' is specified in a "
                                        "subset definition but is "
                                        "not included in supplied "
                                        "lane master list" % l)

        # Adapter sequences
        sample_sheet = SampleSheet(self._sample_sheet)
        if adapter_sequence is None:
            try:
                adapter_sequence = sample_sheet.settings['Adapter']
            except KeyError:
                adapter_sequence = ""
        if adapter_sequence_read2 is None:
            try:
                adapter_sequence_read2 = sample_sheet.settings['AdapterRead2']
            except KeyError:
                adapter_sequence_read2 = ""

        # ICELL8 well list
        if icell8_well_list:
            icell8_well_list = os.path.abspath(icell8_well_list)
        
        # Defaults
        self._bcl_converter = bcl_converter
        self._bases_mask = bases_mask
        self._adapter_sequence = adapter_sequence
        self._adapter_sequence_read2 = adapter_sequence_read2
        self._minimum_trimmed_read_length = minimum_trimmed_read_length
        self._mask_short_adapter_reads = mask_short_adapter_reads
        self._icell8_well_list = icell8_well_list
        self._icell8_atac_swap_i1_and_i2 = icell8_atac_swap_i1_and_i2
        self._icell8_atac_reverse_complement = \
                                icell8_atac_reverse_complement
        self._trim_adapters = bool(trim_adapters)
        self._fastq_statistics = bool(fastq_statistics)
        self._analyse_barcodes = bool(analyse_barcodes)

        # Define parameters
        self.add_param('data_dir',value=run_dir,type=str)
        self.add_param('sample_sheet',value=self._sample_sheet,type=str)
        self.add_param('analysis_dir',type=str)
        self.add_param('out_dir',type=str)
        self.add_param('primary_data_dir',type=str)
        self.add_param('barcode_analysis_dir',type=str)
        self.add_param('counts_dir',type=str)
        self.add_param('qc_report',type=str)
        self.add_param('force_copy_of_primary_data',value=False,type=bool)
        self.add_param('no_lane_splitting',value=False,type=bool)
        self.add_param('create_fastq_for_index_read',value=False,type=bool)
        self.add_param('find_adapters_with_sliding_window',value=False,
                       type=bool)
        self.add_param('create_empty_fastqs',value=False,type=bool)
        self.add_param('name',type=str)
        self.add_param('stats_file',type=str)
        self.add_param('stats_full',type=str)
        self.add_param('per_lane_stats',type=str)
        self.add_param('per_lane_sample_stats',type=str)
        self.add_param('nprocessors',type=int)
        self.add_param('require_bcl2fastq',type=str)
        self.add_param('require_bclconvert',type=str)
        self.add_param('cellranger_jobmode',value='local',type=str)
        self.add_param('cellranger_mempercore',type=int)
        self.add_param('cellranger_maxjobs',type=int)
        self.add_param('cellranger_jobinterval',type=int)
        self.add_param('cellranger_localcores',type=int)
        self.add_param('cellranger_localmem',type=int)

        # Internal parameters
        self.add_param('_platform')
        self.add_param('_bcl2fastq_info')
        self.add_param('_bclconvert_info')
        self.add_param('_cellranger_info')
        self.add_param('_cellranger_atac_info')
        self.add_param('_cellranger_arc_info')
        self.add_param('_spaceranger_info')
        self.add_param('_missing_fastqs',type=list)

        # Define runners
        self.add_runner('rsync_runner')
        self.add_runner('bcl2fastq_runner')
        self.add_runner('bclconvert_runner')
        self.add_runner('demultiplex_icell8_atac_runner')
        self.add_runner('cellranger_runner')
        self.add_runner('cellranger_atac_runner')
        self.add_runner('cellranger_arc_runner')
        self.add_runner('spaceranger_runner')
        self.add_runner('stats_runner')

        # Define module environment modules
        self.add_envmodules('bcl2fastq')
        self.add_envmodules('bclconvert')
        self.add_envmodules('cellranger_mkfastq')
        self.add_envmodules('cellranger_atac_mkfastq')
        self.add_envmodules('cellranger_arc_mkfastq')
        self.add_envmodules('spaceranger_mkfastq')

        # Pipeline outputs
        self.add_output('platform',self.params._platform)
        self.add_output('primary_data_dir',self.params.primary_data_dir)
        self.add_output('acquired_primary_data',Param())
        self.add_output('bcl2fastq_info',self.params._bcl2fastq_info)
        self.add_output('bclconvert_info',self.params._bclconvert_info)
        self.add_output('cellranger_info',self.params._cellranger_info)
        self.add_output('cellranger_atac_info',
                        self.params._cellranger_atac_info)
        self.add_output('spaceranger_info',self.params._spaceranger_info)
        self.add_output('cellranger_arc_info',
                        self.params._cellranger_arc_info)
        self.add_output('stats_file',Param())
        self.add_output('stats_full',Param())
        self.add_output('per_lane_stats',Param())
        self.add_output('per_lane_sample_stats',Param())
        self.add_output('missing_fastqs',self.params._missing_fastqs)

        # Lane subsets
        self._subsets = []

        # Create default lane subsets based on the lanes and
        # index sequences in the sample sheet
        sample_sheet_data = SampleSheet(self._sample_sheet)
        if sample_sheet_data.has_lanes:
            # Sample sheet has lanes
            # Lanes with the same "masked index sequence"
            # will be grouped together
            masked_indexes = dict()
            for data in sample_sheet_data:
                # Get index sequence
                index_sequence = samplesheet_index_sequence(data)
                # Filter out lanes not in restricted lane set
                lane = int(data['Lane'])
                if self._use_lanes:
                    if lane not in self._use_lanes:
                        # Lane is excluded so drop this line
                        continue
                # Store lanes against unique masked index sequences
                masked_index = self._mask_sequence(index_sequence)
                try:
                    masked_indexes[masked_index].add(lane)
                except KeyError:
                    masked_indexes[masked_index] = set((lane,))
            # Build subsets grouped by masked index
            for s in masked_indexes:
                lanes = sorted(list(masked_indexes[s]))
                self._add_subset(lanes=lanes,
                                 protocol=protocol,
                                 masked_index=s)
        else:
            # No lanes in sample sheet
            index_sequence = samplesheet_index_sequence(
                sample_sheet_data[0])
            masked_index = self._mask_sequence(index_sequence)
            self._add_subset(lanes=[],
                             protocol=protocol,
                             masked_index=masked_index)

        # Set pipeline defaults for automatically generated
        # lane subsets
        for s in self.subsets:
            self._update_subset(
                s,
                bcl_converter=self._bcl_converter,
                bases_mask=self._bases_mask,
                trim_adapters=self._trim_adapters,
                minimum_trimmed_read_length=\
                self._minimum_trimmed_read_length,
                mask_short_adapter_reads=\
                self._mask_short_adapter_reads,
                adapter_sequence=self._adapter_sequence,
                adapter_sequence2=self._adapter_sequence_read2,
                no_lane_splitting=self.params.no_lane_splitting,
                create_fastq_for_index_read=\
                self.params.create_fastq_for_index_read,
                find_adapters_with_sliding_window=\
                self.params.find_adapters_with_sliding_window,
                icell8_well_list=self._icell8_well_list,
                icell8_atac_swap_i1_and_i2=\
                self._icell8_atac_swap_i1_and_i2,
                icell8_atac_reverse_complement=\
                self._icell8_atac_reverse_complement,
                analyse_barcodes=self._analyse_barcodes
            )

        # Add user-defined subsets
        # The pipeline defaults will be inherited
        if lane_subsets:
            for user_subset in lane_subsets:
                # Look for an existing subset
                lanes = user_subset['lanes']
                assigned_subset = None
                for s in self.subsets:
                    if lanes == s['lanes']:
                        # Exact match
                        assigned_subset = s
                        break
                    else:
                        # Is it a subset of a subset?
                        if all([l in s['lanes'] for l in lanes]):
                            # Part of an existing subset
                            # which will be split
                            updated_lanes = [l for l in s['lanes']
                                             if l not in lanes]
                            # Make two new subsets copied from
                            # the superset, with lanes updated
                            self._copy_subset(updated_lanes,s)
                            assigned_subset = self._copy_subset(lanes,s)
                            # Remove the original superset
                            self._remove_subset(s['lanes'])
                            break
                # Check a subset was located or created
                if assigned_subset is None:
                    raise Exception("Failed to assign a lane subset for "
                                    "%s" % lanes)
                # Assign protocol to user subset
                try:
                    self._update_subset(assigned_subset,
                                        protocol=user_subset['protocol'])
                except KeyError:
                    pass

        # Check that subsets don't split projects
        if self._sample_sheet and len(self.subsets) > 1:
            if self._subsets_split_project():
                raise Exception("Subsets would split a project")

        # Check that barcodes are consistent with protocols
        for s in self.subsets:
            protocol = s['protocol']
            if protocol not in PROTOCOLS:
                raise Exception("Protocol '%s': not recognised" %
                                protocol)
            masked_index = s['masked_index']
            if masked_index == '__10X__':
                # 10xGenomics barcodes
                if protocol not in PROTOCOLS_10X:
                    raise Exception("Protocol '%s': can't handle "
                                    "10xGenomics barcodes" % protocol)
            else:
                # Standard barcodes
                if protocol in PROTOCOLS_10X:
                    raise Exception("Protocol '%s': needs 10xGenomics "
                                    "barcodes ('%s' not valid)" %
                                    (protocol,masked_index))
            
        # Update parameters on each subset according to the
        # assigned protocol (overriding pipeline defaults)
        for s in self.subsets:
            protocol = s['protocol']
            if protocol == 'mirna':
                # miRNA-seq protocol
                # Set minimum trimmed read length and turn off masking
                self._update_subset(s,
                                    minimum_trimmed_read_length=10,
                                    mask_short_adapter_reads=0)
            elif protocol == 'icell8':
                # ICELL8 protocol
                # Set minimum trimmed read length and turn off masking
                self._update_subset(s,
                                    minimum_trimmed_read_length=21,
                                    mask_short_adapter_reads=0)
            elif protocol == 'icell8_atac':
                # ICELL8 single-cell ATAC-seq
                self._update_subset(s,
                                    create_fastq_for_index_read=True)
            elif protocol == '10x_chromium_sc':
                # 10xGenomics Chromium SC
                # Disable adapter trimming
                self._update_subset(s,
                                    trim_adapters=False)
            elif protocol == '10x_atac':
                # 10xGenomics ATAC-seq
                # Disable adapter trimming
                self._update_subset(s,
                                    trim_adapters=False)
            elif protocol == '10x_visium':
                # 10xGenomics Visium
                # Disable adapter trimming
                self._update_subset(s,
                                    trim_adapters=False)
            elif protocol == '10x_multiome':
                # 10xGenomics multiome
                # Disable adapter trimming
                self._update_subset(s,
                                    trim_adapters=False)
            
        # Finally update parameters for user-defined
        # lane subsets (overriding both pipeline and
        # protocol defaults)
        if lane_subsets:
            for s in lane_subsets:
                lanes = s['lanes']
                assigned_subset = self._fetch_subset(lanes)
                self._update_subset(assigned_subset,
                                    **{ kw: s[kw] for kw in s
                                        if (kw != 'lanes' and
                                            kw != 'protocol') })

        # Perform checks for subsets
        for s in self.subsets:
            if s['protocol'] == 'icell8_atac':
                # ICELL8 ATAC
                # Check well list file is defined
                if not s['icell8_well_list']:
                    raise Exception("No ICELL8 well list assigned for "
                                    "lanes %s" % s['lanes'])

        # Reset lanes for single subset which implicitly
        # includes all lanes in the run
        if len(self.subsets) == 1 and not self._use_lanes:
            self.subsets[0]['lanes'] = []

        # Build the pipeline
        self._build_pipeline()

    @property
    def subsets(self):
        """
        Return list of lane subsets defined in pipeline
        """
        #return self._subsets
        return sorted(self._subsets,
                      key=lambda x: x['lanes'][0]
                      if len(x['lanes']) else 0)

    def _mask_sequence(self,s):
        """
        Internal: return 'masked' version of index sequence

        For standard index sequences (e.g. 'ATGGAT',
        'ATTGGGTA-TTATCCCA'), returns the sequence with
        A,C,G and T's replaced by N (e.g. 'NNNNNN',
        'NNNNNNNN-NNNNNNNN').

        For sequences matching 10xGenomics indices, returns
        the string '__10X__'.

        For empty index sequences, returns the string
        '__NO_SEQUENCE__'.

        Arguments:
          s (str): index sequence

        Returns:
          String: masked version of supplied index sequence.
        """
        if s is None or s == "":
            # No sequence
            return "__NO_SEQUENCE__"
        if s.startswith('SI-'):
            # Possible 10xGenomics index sequence
            try:
                if barcode_is_10xgenomics('-'.join(s.split('-')[0:3])):
                    return "__10X__"
            except IndexError:
                pass
        for c in "ACGT":
            # Mask bases with 'N's
            s = s.replace(c,'N')
        return s

    def _sample_sheet_is_valid(self,sample_sheet):
        """
        Internal: checks that sample sheet is valid

        Arguments:
          sample_sheet (str): path to sample sheet file

        Returns:
          Boolean: True if sample sheet has valid barcodes
            and doesn't contain non-ASCII characters; False
            otherwise.
        """
        # Check for invalid barcodes
        invalid_barcodes = SampleSheetLinter(
            sample_sheet_file=self._sample_sheet).has_invalid_barcodes()
        if invalid_barcodes:
            logger.error("Invalid barcodes detected")
            for line in invalid_barcodes:
                logger.critical("%s" % line)
        # Check for invalid characters
        invalid_characters = SampleSheetLinter(
            sample_sheet_file=self._sample_sheet).has_invalid_characters()
        if invalid_characters:
            logger.critical("Invalid non-printing/non-ASCII characters "
                            "detected")
        # Return boolean indicating if there's a problem
        if invalid_barcodes or invalid_characters:
            return False
        else:
            return True

    def _subsets_split_project(self):
        """
        Internal: checks if subsets split any projects

        Returns:
          Boolean: True if a project from the sample sheet
            is split between two or more subsets; False
            otherwise.
        """
        # Need at least one subset
        if len(self.subsets) == 1:
            return False
        # Load sample sheet
        ss = SampleSheet(self._sample_sheet)
        # Sample sheet must contain lane information
        if not ss.has_lanes:
            return False
        sample_project = ss.sample_project_column
        # Check projects in subsets
        projects = set()
        for s in self.subsets:
            subset_projects = set()
            for lane in s['lanes']:
                for data in ss.data:
                    if lane == data['Lane']:
                        subset_projects.add(data[sample_project])
                        break
            for project in subset_projects:
                if project in projects:
                    logger.critical("Project '%s' appears in multiple "
                                    "lane subsets" % project)
                    return True
            projects.update(subset_projects)
        # No projects split
        return False

    def _add_subset(self,lanes,**kws):
        """
        Internal: creates a new lane subset

        Arguments:
          lanes (list): list of lane numbers in the
            subset
          kws (mapping): optional list of key-value
            pairs assigning values to parameters for
            the new subset

        Returns:
          Dictionary: new subset
        """
        lanes_ = sorted([int(l) for l in lanes])
        s = subset(lanes=lanes_,**kws)
        self._subsets.append(s)
        return s

    def _fetch_subset(self,lanes):
        """
        Internal: return existing lane subset

        Arguments:
          lanes (list): set of lane numbers to
            match; all specified lanes must be 
            present in a subset for it to be
            returned

        Returns:
          Dictionary: subset matching the lane
            specification

        Raises:
          KeyError (if no matching subset is
            located)
        """
        lanes_ = sorted([int(l) for l in lanes])
        for s in self._subsets:
            if s['lanes'] == lanes_:
                return s
        raise KeyError("No subset matching '%s'" % lanes)

    def _remove_subset(self,lanes):
        """
        Internal: remove a subset from the pipeline

        Arguments:
          lanes (list): set of lane numbers to
            match; all specified lanes must be 
            present in a subset for it to be
            removed

        Raises:
          KeyError (if a subset is not located for
            removal)
        """
        lanes_ = sorted([int(l) for l in lanes])
        for i,s in enumerate(self._subsets):
            if s['lanes'] == lanes_:
                del(self._subsets[i])
                return
        raise KeyError("No subset matching '%s'" % lanes)

    def _copy_subset(self,lanes,existing_subset):
        """
        Internal: duplicate an existing subset

        Creates a new subset grouping the specified
        lanes and with all other parameters copied
        from the supplied subset.

        Arguments:
          lanes (list): list of lane numbers in the
            new subset
          existing_subset (dictionary): existing
            subset to copy parameters from

        Returns:
          Dictionary: duplicated subset
        """
        data = { kw: existing_subset[kw]
                 for kw in existing_subset if kw != 'lanes' }
        return self._add_subset(lanes,**data)

    def _update_subset(self,s,**kws):
        """
        Internal: update parameters in a lane subset

        Arguments:
          s (dictionary): lane subset to update
          kws (mapping): optional set of key-value
            pairs assigning new values to parameters
            for the subset
        """
        for attr in kws:
            if attr not in LANE_SUBSET_ATTRS:
                raise KeyError("Unsupported subset attribute '%s'"
                               % attr)
            s[attr] = kws[attr]

    def _build_pipeline(self):
        """
        Internal: construct the pipeline

        Assembles the pipeline tasks according to the
        data provided at instantiation
        """
        ####################
        # Deal with platform
        ####################
        if self._platform is None:
            # Try to get the platform from the run name
            self._platform = get_sequencer_platform(self._run_dir)
        if self._platform is None:
            # Set a generic platform name if it can't be identified
            self._platform = "illumina"

        #################
        # Report subsets
        #################
        self.report("Building pipeline for lane subsets:")
        for s in self.subsets:
            self.report("- Lanes: %s" % ','.join([str(l)
                                                  for l in s['lanes']]))
            for attr in s:
                if attr != 'lanes' and s[attr] is not None:
                    self.report("- %s: %s" % (attr,s[attr]))

        #####################
        # Fetch primary data
        #####################
        fetch_primary_data = FetchPrimaryData(
            "Fetch primary data",
            self.params.data_dir,
            self.params.primary_data_dir,
            force_copy=self.params.force_copy_of_primary_data
        )
        self.add_task(fetch_primary_data,
                      runner=self.runners['rsync_runner'])
        identify_platform = IdentifyPlatform(
            "Identify sequencer platform",
            fetch_primary_data.output.run_dir,
            platform=self._platform
        )
        self.add_task(identify_platform)

        # Load sample sheet data
        sample_sheet = SampleSheet(self._sample_sheet)

        # Keep track of Fastq generation output directories
        fastq_out_dirs = []

        # Keep track of missing fastqs
        missing_fastqs = []

        # Keep track of lanes to analyse barcodes for
        lanes_for_barcode_analysis = []

        # Placeholders for tasks acquiring software versions
        get_bcl2fastq = None
        get_bclconvert = None
        get_bcl2fastq_for_10x = None
        get_bcl2fastq_for_10x_atac = None
        get_bcl2fastq_for_10x_visium = None
        get_bcl2fastq_for_10x_multiome = None
        get_cellranger = None
        get_cellranger_atac = None
        get_cellranger_arc = None
        get_spaceranger = None

        #########################
        # Merge Fastq directories
        #########################
        merge_fastq_dirs = MergeFastqDirs(
            "Merge BCL to Fastq output directories",
            fastq_out_dirs,
            self.params.out_dir
        )
        self.add_task(merge_fastq_dirs)

        ##########################
        # Statistics and reporting
        ##########################
        if self._fastq_statistics:
            # Generate statistics
            fastq_statistics = FastqStatistics(
                "Generate statistics for Fastqs",
                self.params.out_dir,
                self._sample_sheet,
                self.params.analysis_dir,
                stats_file=self.params.stats_file,
                stats_full_file=self.params.stats_full,
                per_lane_stats_file=self.params.per_lane_stats,
                per_lane_sample_stats_file=\
                self.params.per_lane_sample_stats,
                nprocessors=self.params.nprocessors)
            self.add_task(fastq_statistics,
                          runner=self.runners['stats_runner'],
                          requires=(merge_fastq_dirs,))
            # Processing QC report
            report_qc = ReportProcessingQC(
                "Report Processing QC",
                name=self.params.name,
                analysis_dir=self.params.analysis_dir,
                stats_file=fastq_statistics.output.stats_full,
                per_lane_stats_file=\
                fastq_statistics.output.per_lane_stats,
                per_lane_sample_stats_file=\
                fastq_statistics.output.per_lane_sample_stats,
                report_html=self.params.qc_report
            )
            self.add_task(report_qc)

        ##################
        # Barcode analysis
        ##################
        if not self._fastq_statistics:
            do_barcode_analysis = False
        elif len(self.subsets) == 1 and not self.subsets[0]['lanes']:
            try:
                do_barcode_analysis = self.subsets[0]['analyse_barcodes']
            except KeyError:
                do_barcode_analysis = self._analyse_barcodes
        else:
            # Determine lanes to perform barcode analysis for
            for subset in self.subsets:
                if subset['analyse_barcodes']:
                    if subset['lanes']:
                        lanes_for_barcode_analysis.extend(subset['lanes'])
            if lanes_for_barcode_analysis:
                do_barcode_analysis = True
                lanes_for_barcode_analysis = sorted(
                    list(set(lanes_for_barcode_analysis)))
            else:
                do_barcode_analysis = False
        if do_barcode_analysis:
            # Set up pipeline for barcode analysis
            self.report("Lanes for barcode analysis: %s" %
                        ','.join([str(l) for l in lanes_for_barcode_analysis]))
            analyse_barcodes = AnalyseBarcodes(
                sample_sheet=self._sample_sheet)
            self.add_pipeline(analyse_barcodes,
                              params={
                                  'bcl2fastq_dir': self.params.out_dir,
                                  'title': FunctionParam(
                                      lambda run_dir,name:
                                      "Barcode analysis for %s%s" %
                                      (os.path.basename(run_dir),
                                       ' (%s)' % name if name else ''),
                                      self._run_dir,
                                      self.params.name),
                                  'lanes': Param(
                                      value=lanes_for_barcode_analysis),
                              },
                              requires=(merge_fastq_dirs,))

        # For each subset, add the appropriate set of
        # tasks for the protocol
        for subset in self.subsets:

            # Lanes in this subset
            lanes = subset['lanes']
            self.report("Adding tasks for %s" %
                        ("all lanes" if not lanes
                         else "lanes %s" % ','.join([str(l) for l in lanes])))

            # Protocol
            protocol = subset['protocol']
            self.report("- Protocol: %s" % protocol)

            # BCL to Fastq converter
            converter = subset['bcl_converter']

            #############
            # Bases mask
            #############
            bases_mask = subset['bases_mask']

            ###################
            # Adapter trimming
            ###################
            trim_adapters = subset['trim_adapters']

            if trim_adapters:
                minimum_trimmed_read_length = \
                            subset['minimum_trimmed_read_length']
                mask_short_adapter_reads = \
                            subset['mask_short_adapter_reads']
                adapter_sequence = subset['adapter_sequence']
                adapter_sequence_read2 = subset['adapter_sequence2']
            else:
                # No adapter trimming
                minimum_trimmed_read_length = 0
                mask_short_adapter_reads = 0
                adapter_sequence = ""
                adapter_sequence_read2 = ""
            # Report adapter settings
            self.report("- Adapter sequence      : %s" % (adapter_sequence
                                                          if adapter_sequence
                                                          else "<none>"))
            if adapter_sequence_read2:
                self.report("  Adapter sequence read2: %s" %
                            adapter_sequence_read2)

            #################
            # Lane splitting
            #################
            no_lane_splitting = subset['no_lane_splitting']

            #####################
            # Create index reads
            #####################
            create_fastq_for_index_read = \
                subset['create_fastq_for_index_read']

            # Use sliding window for adapter trimming
            find_adapters_with_sliding_window = \
                subset['find_adapters_with_sliding_window']

            #########
            # ICELL8
            #########
            well_list = subset['icell8_well_list']
            swap_i1_and_i2 = subset['icell8_atac_swap_i1_and_i2']
            reverse_complement = subset['icell8_atac_reverse_complement']
            
            ###################
            # Make samplesheet
            ###################
            make_sample_sheet = MakeSampleSheet(
                "Make sample sheet%s" % (" for lanes %s" %
                                         ','.join([str(x)  for x in lanes])
                                         if lanes else ""),
                self._sample_sheet,
                lanes=lanes,
                adapter=adapter_sequence,
                adapter_read2=adapter_sequence_read2)
            self.add_task(make_sample_sheet)

            # Construct a name for the output directory
            if lanes and len(self.subsets) > 1:
                # If lanes were specified and there are multiple
                # subsets then this will be
                #
                # ANALYSIS_DIR/BASENAME(OUT_DIR).LANES
                lanes_id = ".L%s" % ''.join([str(l) for l in lanes])
            else:
                # If no lanes specified and/or only one
                # subsets it will be
                #
                # ANALYSIS_DIR/BASENAME(OUT_DIR)
                lanes_id = ""
            fastq_out_dir = FunctionParam(
                lambda analysis_dir,out_dir,lanes_id:
                os.path.join(analysis_dir,
                             "%s%s" % (os.path.basename(out_dir),
                                       lanes_id)),
                self.params.analysis_dir,
                self.params.out_dir,
                lanes_id)

            # Flag if the final output directory exists
            final_output_exists = PathExistsParam(
                PathJoinParam(self.params.out_dir))
                
            ###############
            # Set up tasks
            ###############

            # Attempt to restore backup for this lane set
            restore_backup = RestoreBackupDirectory(
                "Restore backup%s" % (" for lanes %s" %
                                      ','.join([str(x) for x in lanes])
                                      if lanes else ""),
                fastq_out_dir,
                skip_restore=final_output_exists)
            self.add_task(restore_backup)

            # Standard protocols
            if protocol in ("standard","mirna"):

                if converter == "bcl2fastq":
                    # Get bcl2fastq information
                    if get_bcl2fastq is None:
                        # Create a new task only if one doesn't already
                        # exist
                        get_bcl2fastq = GetBcl2Fastq(
                            "Get information on bcl2fastq",
                            require_bcl2fastq=self.params.require_bcl2fastq)
                        self.add_task(get_bcl2fastq,
                                      envmodules=self.envmodules['bcl2fastq'])
                    # Run standard bcl2fastq
                    make_fastqs = RunBcl2Fastq(
                        "Run bcl2fastq%s" %
                        (" for lanes %s" % ','.join([str(x) for x in lanes])
                         if lanes else ""),
                        fetch_primary_data.output.run_dir,
                        fastq_out_dir,
                        make_sample_sheet.output.custom_sample_sheet,
                        bases_mask=bases_mask,
                        minimum_trimmed_read_length=\
                        minimum_trimmed_read_length,
                        mask_short_adapter_reads=\
                        mask_short_adapter_reads,
                        nprocessors=self.params.nprocessors,
                        no_lane_splitting=self.params.no_lane_splitting,
                        create_fastq_for_index_read=\
                        create_fastq_for_index_read,
                        find_adapters_with_sliding_window=\
                        find_adapters_with_sliding_window,
                        create_empty_fastqs=self.params.create_empty_fastqs,
                        platform=identify_platform.output.platform,
                        bcl2fastq_exe=get_bcl2fastq.output.bcl2fastq_exe,
                        bcl2fastq_version=\
                        get_bcl2fastq.output.bcl2fastq_version,
                        skip_bcl2fastq=final_output_exists)
                    # Add the Fastq generation task to the pipeline
                    self.add_task(make_fastqs,
                                  runner=self.runners['bcl2fastq_runner'],
                                  envmodules=self.envmodules['bcl2fastq'],
                                  requires=(restore_backup,))
                elif converter == "bclconvert":
                    # Get BCL Convert information
                    if get_bclconvert is None:
                        # Create a new task only if one doesn't already
                        # exist
                        get_bclconvert = GetBclConvert(
                            "Get information on BCL Convert",
                            require_version=self.params.require_bclconvert)
                        self.add_task(get_bclconvert,
                                      envmodules=self.envmodules['bclconvert'])
                    # Parameter to store intermediate output
                    # directories for merging
                    fastq_lanes_out_dirs = ListParam()
                    if lanes:
                        # Subset with multiple lanes
                        #
                        # Divide Fastq generation into individual lanes
                        # then merge outputs at the end
                        make_fastqs = MergeFastqs(
                            "Merge Fastqs from lanes %s"
                            % ','.join([str(l) for l in lanes]),
                            fastq_lanes_out_dirs,
                            fastq_out_dir,
                            make_sample_sheet.output.custom_sample_sheet,
                            no_lane_splitting=self.params.no_lane_splitting,
                            create_empty_fastqs=self.params.create_empty_fastqs)
                        self.add_task(make_fastqs)
                        # Run BCL Convert for each lane in subset
                        for lane in lanes:
                            # Output dir for this lane
                            fastq_lane_out_dir = FunctionParam(
                                lambda working_dir,out_dir,lane:
                                os.path.join(
                                    working_dir,
                                    "%s.L%s" % (os.path.basename(out_dir),
                                                lane)),
                                self.params.WORKING_DIR,
                                self.params.out_dir,
                                lane)
                            fastq_lanes_out_dirs.append(fastq_lane_out_dir)
                            # Run BCL Convert
                            bcl_convert = RunBclConvert(
                                "Run BCL Convert for lane %s" % lane,
                                fetch_primary_data.output.run_dir,
                                fastq_lane_out_dir,
                                make_sample_sheet.output.custom_sample_sheet,
                                lane=lane,
                                bases_mask=bases_mask,
                                minimum_trimmed_read_length=\
                                minimum_trimmed_read_length,
                                mask_short_adapter_reads=\
                                mask_short_adapter_reads,
                                nprocessors=self.params.nprocessors,
                                no_lane_splitting=False,
                                create_fastq_for_index_read=\
                                create_fastq_for_index_read,
                                create_empty_fastqs=False,
                                ignore_missing_fastqs=True,
                                platform=identify_platform.output.platform,
                                bclconvert_exe=\
                                get_bclconvert.output.bclconvert_exe,
                                bclconvert_version=\
                                get_bclconvert.output.bclconvert_version,
                                skip_bclconvert=final_output_exists)
                            # Add the Fastq generation task to the pipeline
                            self.add_task(
                                bcl_convert,
                                runner=self.runners['bclconvert_runner'],
                                envmodules=self.envmodules['bclconvert'],
                                requires=(restore_backup,))
                            # Merging task depends on BCL Convert
                            make_fastqs.requires(bcl_convert)
                    else:
                        # Subset with no lanes specified
                        #
                        # Handle no lane splitting
                        if sample_sheet.has_lanes:
                            # BCL Convert cannot handle no lane splitting
                            # directly for runs with "asymmetric" lanes
                            # (e.g. HISeq) i.e. platforms where lanes are
                            # explicitly specified in the sample sheet
                            #
                            # Delegate no lane splitting to final Fastq
                            # merging task
                            no_lane_splitting = False
                        else:
                            # Let BCL convert handle no lane splitting
                            # (final Fastq merging task will do nothing)
                            no_lane_splitting = \
                                self.params.no_lane_splitting
                        # Intermediate output Fastq dir
                        tmp_fastq_out_dir = FunctionParam(
                            lambda working_dir,out_dir:
                            os.path.join(
                                working_dir,
                                "__%s.tmp" % os.path.basename(out_dir)),
                            self.params.WORKING_DIR,
                            self.params.out_dir)
                        # Run BCL Convert for all lanes
                        bcl_convert = RunBclConvert(
                            "Run BCL Convert",
                            fetch_primary_data.output.run_dir,
                            tmp_fastq_out_dir,
                            make_sample_sheet.output.custom_sample_sheet,
                            bases_mask=bases_mask,
                            minimum_trimmed_read_length=\
                            minimum_trimmed_read_length,
                            mask_short_adapter_reads=\
                            mask_short_adapter_reads,
                            nprocessors=self.params.nprocessors,
                            no_lane_splitting=no_lane_splitting,
                            create_fastq_for_index_read=\
                            create_fastq_for_index_read,
                            create_empty_fastqs=\
                            self.params.create_empty_fastqs,
                            platform=identify_platform.output.platform,
                            bclconvert_exe=\
                            get_bclconvert.output.bclconvert_exe,
                            bclconvert_version=\
                            get_bclconvert.output.bclconvert_version,
                            skip_bclconvert=final_output_exists)
                        # Add the Fastq generation task to the pipeline
                        self.add_task(
                            bcl_convert,
                            runner=self.runners['bclconvert_runner'],
                            envmodules=self.envmodules['bclconvert'],
                            requires=(restore_backup,))
                        # Merge Fastqs (if requested) and move to
                        # final location
                        make_fastqs = MergeFastqs(
                            "Merge Fastqs across lanes",
                            ListParam((tmp_fastq_out_dir,)),
                            fastq_out_dir,
                            make_sample_sheet.output.custom_sample_sheet,
                            no_lane_splitting=\
                            self.params.no_lane_splitting,
                            create_empty_fastqs=\
                            self.params.create_empty_fastqs)
                        self.add_task(make_fastqs,
                                      requires=(bcl_convert,))
            # ICELL8 RNA-seq
            if protocol == "icell8":
                # Get bcl2fastq information
                if get_bcl2fastq is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_bcl2fastq = GetBcl2Fastq(
                        "Get information on bcl2fastq",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq,
                                  envmodules=self.envmodules['bcl2fastq'])
                # Get ICELL8 bases mask
                get_bases_mask = GetBasesMaskIcell8(
                    "Get bases mask for ICELL8",
                    fetch_primary_data.output.run_dir,
                    make_sample_sheet.output.custom_sample_sheet)
                self.add_task(get_bases_mask)
                # Run standard bcl2fastq
                make_fastqs = RunBcl2Fastq(
                    "Run bcl2fastq%s" % (" for lanes %s" %
                                         ','.join([str(x) for x in lanes])
                                         if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    fastq_out_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    bases_mask=get_bases_mask.output.bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    nprocessors=self.params.nprocessors,
                    no_lane_splitting=self.params.no_lane_splitting,
                    create_fastq_for_index_read=\
                    create_fastq_for_index_read,
                    find_adapters_with_sliding_window=\
                    find_adapters_with_sliding_window,
                    create_empty_fastqs=self.params.create_empty_fastqs,
                    platform=identify_platform.output.platform,
                    bcl2fastq_exe=get_bcl2fastq.output.bcl2fastq_exe,
                    bcl2fastq_version=get_bcl2fastq.output.bcl2fastq_version,
                    skip_bcl2fastq=final_output_exists)
                self.add_task(make_fastqs,
                              runner=self.runners['bcl2fastq_runner'],
                              envmodules=self.envmodules['bcl2fastq'],
                              requires=(restore_backup,))

            # ICELL8 ATAC-seq
            if protocol == "icell8_atac":
                # Parameter to indicate if bcl2fastq needs to run
                skip_bcl2fastq = FunctionParam(
                    lambda x,y: x or y,
                    final_output_exists,
                    PathExistsParam(fastq_out_dir))
                
                # Temporary dir for intermediate Fastqs
                tmp_bcl2fastq_dir = "bcl2fastq.icell8_atac%s" % lanes_id
                # Get bcl2fastq information
                if get_bcl2fastq is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_bcl2fastq = GetBcl2Fastq(
                        "Get information on bcl2fastq",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq,
                                  envmodules=self.envmodules['bcl2fastq'])
                # Get ICELL8 bases mask
                get_bases_mask = GetBasesMaskIcell8Atac(
                    "Get bases mask for ICELL8",
                    fetch_primary_data.output.run_dir)
                self.add_task(get_bases_mask)
                # Run standard bcl2fastq
                run_bcl2fastq = RunBcl2Fastq(
                    "Run bcl2fastq%s" % (" for lanes %s" %
                                         ','.join([str(x) for x in lanes])
                                         if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    tmp_bcl2fastq_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    bases_mask=get_bases_mask.output.bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    nprocessors=self.params.nprocessors,
                    no_lane_splitting=self.params.no_lane_splitting,
                    create_fastq_for_index_read=\
                    create_fastq_for_index_read,
                    find_adapters_with_sliding_window=\
                    find_adapters_with_sliding_window,
                    create_empty_fastqs=self.params.create_empty_fastqs,
                    bcl2fastq_exe=get_bcl2fastq.output.bcl2fastq_exe,
                    bcl2fastq_version=get_bcl2fastq.output.bcl2fastq_version,
                    skip_bcl2fastq=skip_bcl2fastq)
                self.add_task(run_bcl2fastq,
                              runner=self.runners['bcl2fastq_runner'],
                              envmodules=self.envmodules['bcl2fastq'],
                              requires=(restore_backup,))
                # Demultiplex ICELL8 ATAC Fastqs
                make_fastqs = DemultiplexIcell8Atac(
                    "Demultiplex ICELL8 ATAC fastqs%s" %
                    (" for lanes %s" % ','.join([str(x) for x in lanes])
                     if lanes else ""),
                    tmp_bcl2fastq_dir,
                    fastq_out_dir,
                    well_list,
                    nprocessors=self.params.nprocessors,
                    swap_i1_and_i2=swap_i1_and_i2,
                    reverse_complement=reverse_complement,
                    skip_demultiplex=final_output_exists)
                self.add_task(
                    make_fastqs,
                    runner=self.runners['demultiplex_icell8_atac_runner'],
                    requires=(run_bcl2fastq,))

            # 10x RNA-seq
            if protocol == "10x_chromium_sc":
                # Get bcl2fastq information
                if get_bcl2fastq_for_10x is None:
                    get_bcl2fastq_for_10x = GetBcl2Fastq(
                        "Get information on bcl2fastq for cellranger",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq_for_10x,
                                  envmodules=\
                                  self.envmodules['cellranger_mkfastq'])
                # Get cellranger information
                if get_cellranger is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_cellranger = Get10xPackage(
                        "Get information on cellranger",
                        require_package="cellranger")
                    self.add_task(get_cellranger,
                                  envmodules=\
                                  self.envmodules['cellranger_mkfastq'])
                # Run cellranger mkfastq
                make_fastqs = Run10xMkfastq(
                    "Run cellranger mkfastq%s" %
                    (" for lanes %s" % ','.join([str(x) for x in lanes])
                     if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    fastq_out_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    platform=identify_platform.output.platform,
                    bases_mask=bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    jobmode=self.params.cellranger_jobmode,
                    mempercore=self.params.cellranger_mempercore,
                    maxjobs=self.params.cellranger_maxjobs,
                    jobinterval=self.params.cellranger_jobinterval,
                    localcores=self.params.cellranger_localcores,
                    localmem=self.params.cellranger_localmem,
                    pkg_exe=get_cellranger.output.package_exe,
                    pkg_version=get_cellranger.output.package_version,
                    bcl2fastq_exe=get_bcl2fastq_for_10x.output.bcl2fastq_exe,
                    bcl2fastq_version=\
                    get_bcl2fastq_for_10x.output.bcl2fastq_version,
                    skip_mkfastq=final_output_exists
                )
                self.add_task(make_fastqs,
                              runner=self.runners['cellranger_runner'],
                              envmodules=self.envmodules['cellranger_mkfastq'],
                              requires=(restore_backup,))

            # 10x ATAC-seq
            if protocol == "10x_atac":
                # Get bcl2fastq information
                if get_bcl2fastq_for_10x_atac is None:
                    get_bcl2fastq_for_10x_atac = GetBcl2Fastq(
                        "Get information on bcl2fastq for cellranger-atac",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq_for_10x_atac,
                                  envmodules=\
                                  self.envmodules['cellranger_atac_mkfastq'])
                # Get cellranger-atac information
                if get_cellranger_atac is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_cellranger_atac = Get10xPackage(
                        "Get information on cellranger-atac",
                        require_package="cellranger-atac")
                    self.add_task(get_cellranger_atac,
                                  envmodules=\
                                  self.envmodules['cellranger_atac_mkfastq'])
                # Run cellranger mkfastq
                make_fastqs = Run10xMkfastq(
                    "Run cellranger-atac mkfastq%s" %
                    (" for lanes %s" % ','.join([str(x) for x in lanes])
                     if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    fastq_out_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    platform=identify_platform.output.platform,
                    bases_mask=bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    jobmode=self.params.cellranger_jobmode,
                    mempercore=self.params.cellranger_mempercore,
                    maxjobs=self.params.cellranger_maxjobs,
                    jobinterval=self.params.cellranger_jobinterval,
                    localcores=self.params.cellranger_localcores,
                    localmem=self.params.cellranger_localmem,
                    pkg_exe=get_cellranger_atac.output.package_exe,
                    pkg_version=\
                    get_cellranger_atac.output.package_version,
                    bcl2fastq_exe=\
                    get_bcl2fastq_for_10x_atac.output.bcl2fastq_exe,
                    bcl2fastq_version=\
                    get_bcl2fastq_for_10x_atac.output.bcl2fastq_version,
                    skip_mkfastq=final_output_exists
                )
                self.add_task(make_fastqs,
                              runner=self.runners['cellranger_runner'],
                              envmodules=\
                              self.envmodules['cellranger_atac_mkfastq'],
                              requires=(restore_backup,))

            # 10x Visium
            if protocol == "10x_visium":
                # Get bcl2fastq information
                if get_bcl2fastq_for_10x_visium is None:
                    get_bcl2fastq_for_10x_visium = GetBcl2Fastq(
                        "Get information on bcl2fastq for spaceranger",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq_for_10x_visium,
                                  envmodules=\
                                  self.envmodules['spaceranger_mkfastq'])
                # Get spaceranger information
                if get_spaceranger is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_spaceranger = Get10xPackage(
                        "Get information on spaceranger",
                        require_package="spaceranger")
                    self.add_task(get_spaceranger,
                                  envmodules=\
                                  self.envmodules['spaceranger_mkfastq'])
                # Run cellranger mkfastq
                make_fastqs = Run10xMkfastq(
                    "Run spaceranger mkfastq%s" %
                    (" for lanes %s" % ','.join([str(x) for x in lanes])
                     if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    fastq_out_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    platform=identify_platform.output.platform,
                    bases_mask=bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    jobmode=self.params.cellranger_jobmode,
                    mempercore=self.params.cellranger_mempercore,
                    maxjobs=self.params.cellranger_maxjobs,
                    jobinterval=self.params.cellranger_jobinterval,
                    localcores=self.params.cellranger_localcores,
                    localmem=self.params.cellranger_localmem,
                    pkg_exe=get_spaceranger.output.package_exe,
                    pkg_version=\
                    get_spaceranger.output.package_version,
                    bcl2fastq_exe=\
                    get_bcl2fastq_for_10x_visium.output.bcl2fastq_exe,
                    bcl2fastq_version=\
                    get_bcl2fastq_for_10x_visium.output.bcl2fastq_version,
                    skip_mkfastq=final_output_exists
                )
                self.add_task(make_fastqs,
                              runner=self.runners['spaceranger_runner'],
                              envmodules=\
                              self.envmodules['spaceranger_mkfastq'],
                              requires=(restore_backup,))

            # 10x multiome
            if protocol == "10x_multiome":
                # Get bcl2fastq information
                if get_bcl2fastq_for_10x_multiome is None:
                    get_bcl2fastq_for_10x_multiome = GetBcl2Fastq(
                        "Get information on bcl2fastq for cellranger-arc",
                        require_bcl2fastq=self.params.require_bcl2fastq)
                    self.add_task(get_bcl2fastq_for_10x_multiome,
                                  envmodules=\
                                  self.envmodules['cellranger_arc_mkfastq'])
                # Get cellranger-atac information
                if get_cellranger_arc is None:
                    # Create a new task only if one doesn't already
                    # exist
                    get_cellranger_arc = Get10xPackage(
                        "Get information on cellranger-arc",
                        require_package="cellranger-arc")
                    self.add_task(get_cellranger_arc,
                                  envmodules=\
                                  self.envmodules['cellranger_arc_mkfastq'])
                # Run cellranger mkfastq
                make_fastqs = Run10xMkfastq(
                    "Run cellranger-arc mkfastq%s" %
                    (" for lanes %s" % ','.join([str(x) for x in lanes])
                     if lanes else ""),
                    fetch_primary_data.output.run_dir,
                    fastq_out_dir,
                    make_sample_sheet.output.custom_sample_sheet,
                    platform=identify_platform.output.platform,
                    bases_mask=bases_mask,
                    minimum_trimmed_read_length=\
                    minimum_trimmed_read_length,
                    mask_short_adapter_reads=\
                    mask_short_adapter_reads,
                    jobmode=self.params.cellranger_jobmode,
                    mempercore=self.params.cellranger_mempercore,
                    maxjobs=self.params.cellranger_maxjobs,
                    jobinterval=self.params.cellranger_jobinterval,
                    localcores=self.params.cellranger_localcores,
                    localmem=self.params.cellranger_localmem,
                    pkg_exe=get_cellranger_arc.output.package_exe,
                    pkg_version=\
                    get_cellranger_arc.output.package_version,
                    bcl2fastq_exe=\
                    get_bcl2fastq_for_10x_multiome.output.bcl2fastq_exe,
                    bcl2fastq_version=\
                    get_bcl2fastq_for_10x_multiome.output.bcl2fastq_version,
                    skip_mkfastq=final_output_exists
                )
                self.add_task(make_fastqs,
                              runner=self.runners['cellranger_runner'],
                              envmodules=\
                              self.envmodules['cellranger_arc_mkfastq'],
                              requires=(restore_backup,))

            # Add task to list of tasks that downstream
            # tasks need to wait for
            merge_fastq_dirs.requires(make_fastqs)
            # Add output directory to list
            fastq_out_dirs.append(fastq_out_dir)
            # Add missing Fastqs
            missing_fastqs.append(make_fastqs.output.missing_fastqs)

        # Update outputs associated with primary data
        self.output.acquired_primary_data.set(
            PathExistsParam(
                fetch_primary_data.output.run_dir))
        self.params._platform.set(identify_platform.output.platform)

        # Update outputs with bcl2fastq information
        for get_bcl2fastq_task in (get_bcl2fastq,
                                   get_bcl2fastq_for_10x,
                                   get_bcl2fastq_for_10x_atac,
                                   get_bcl2fastq_for_10x_visium,
                                   get_bcl2fastq_for_10x_multiome,):
            if get_bcl2fastq_task:
                get_bcl2fastq = get_bcl2fastq_task
                break
        if get_bcl2fastq:
            self.params._bcl2fastq_info.set(
                get_bcl2fastq.output.bcl2fastq_info)

        # Update outputs with BCL Convert information
        if get_bclconvert:
            self.params._bclconvert_info.set(
                get_bclconvert.output.bclconvert_info)

        # Update outputs with 10x package information
        if get_cellranger:
            self.params._cellranger_info.set(
                get_cellranger.output.package_info)
        if get_cellranger_atac:
            self.params._cellranger_atac_info.set(
                get_cellranger_atac.output.package_info)
        if get_spaceranger:
            self.params._spaceranger_info.set(
                get_spaceranger.output.package_info)
        if get_cellranger_arc:
            self.params._cellranger_arc_info.set(
                get_cellranger_arc.output.package_info)

        # Update outputs associated with stats
        if self._fastq_statistics:
            self.output.stats_file.set(fastq_statistics.output.stats_file)
            self.output.stats_full.set(fastq_statistics.output.stats_full)
            self.output.per_lane_stats.set(
                fastq_statistics.output.per_lane_stats)
            self.output.per_lane_sample_stats.set(
                fastq_statistics.output.per_lane_sample_stats)

        # Update lists of missing Fastqs
        self.params._missing_fastqs.set(
            FunctionParam(self._merge_missing_fastqs,
                          missing_fastqs))

    def _merge_missing_fastqs(self,fastq_sets):
        """
        Internal: merge lists of "missing" Fastqs

        Given a list with lists of Fastqs, combine them
        into a single sorted list of Fastqs.
        """
        missing_fastqs = list()
        for fqs in fastq_sets:
            if fqs:
                missing_fastqs.extend(fqs)
        return sorted(missing_fastqs)

    def run(self,analysis_dir,out_dir=None,barcode_analysis_dir=None,
            primary_data_dir=None,force_copy_of_primary_data=False,
            no_lane_splitting=None,create_fastq_for_index_read=None,
            find_adapters_with_sliding_window=None,
            create_empty_fastqs=None,name=None,stats_file=None,
            stats_full=None,per_lane_stats=None,per_lane_sample_stats=None,
            nprocessors=None,require_bcl2fastq=None,require_bclconvert=None,
            cellranger_jobmode='local',cellranger_mempercore=None,
            cellranger_maxjobs=None,cellranger_jobinterval=None,
            cellranger_localcores=None,cellranger_localmem=None,
            working_dir=None,log_dir=None,log_file=None,
            batch_size=None,batch_limit=None,max_jobs=1,max_slots=None,
            poll_interval=5,runners=None,default_runner=None,
            envmodules=None,verbose=False):
        """
        Run the tasks in the pipeline

        Arguments:
          analysis_dir (str): directory to perform the processing
            and analyses in
          out_dir (str): (sub)directory for output from Fastq
            generation (defaults to 'bcl2fastq')
          barcode_analysis_dir (str): (sub)directory for barcode
            analysis (defaults to 'barcode_analysis')
          primary_data_dir (str): top-level directory
            holding the primary data
          force_copy_of_primary_data (bool): if True then force
            primary data to be copied (rsync'ed) even if it's on
            the local system (default is to link to primary data
            unless it's on a remote filesystem)
          no_lane_splitting (bool): if True then don't split
            output Fastqs across lanes (--no-lane-splitting)
          create_fastq_for_index_read (bool): if True then
            also output Fastqs for the index (I1 etc) reads
            (--create-fastq-for-index-read)
          find_adapters_with_sliding_window (bool): if True then
            use sliding window algorith to identify adapter
            sequences (--find-adapters-with-sliding-window)
          create_empty_fastqs (bool): if True then create empty
            "placeholder" Fastqs if not created by bcl2fastq
          name (str): optional identifier for output
            stats and report files
          stats_file (str): path to statistics output file
          stats_full (str): path to full statistics output file
          per_lane_stats (str): path to per-lane statistics
            output file
          per_lane_sample_stats(str): path to per-lane per-sample
            statistics output file
          nprocessors (int): number of threads to use for
            multithreaded applications (default is to take
            number of CPUs set in job runners)
          require_bcl2fastq (str): if set then specify bcl2fastq
            version requirement; should be a string of the form
           '1.8.4' or '>2.0'. The pipeline will fail if this
           requirement is not met
          require_bclconvert (str): if set then specify bcl-convert
            version requirement; should be a string of the form
           '3.7.5' or '>3.6'. The pipeline will fail if this
           requirement is not met
          cellranger_jobmode (str): job mode to run cellranger in
          cellranger_mempercore (int): memory assumed per core
          cellranger_maxjobs (int): maxiumum number of concurrent
            jobs to run
          cellranger_jobinterval (int): how often jobs are
            submitted (in ms)
          cellranger_localcores (int): maximum number of cores
            cellranger can request in jobmode 'local'
          cellranger_localmem (int): (optional) maximum memory
            cellranger can request in jobmode 'local'
          working_dir (str): optional path to a working
            directory (defaults to temporary directory in
            the current directory)
          log_dir (str): path of directory where log files
            will be written to
          batch_size (int): if set then run commands in
            each task in batches, with each batch running
            this many commands at a time (default is to run
            one command per job)
          batch_limit (int): if set then run commands in
            each task in batches, with the batch size set
            dyanmically so as not to exceed this limit
            (default is to use fixed batch sizes)
          max_jobs (int): optional maximum number of
            concurrent jobs in scheduler (defaults to 1)
          max_slots (int): optional maximum number of 'slots'
            (i.e. concurrent threads or maximum number of
            CPUs) available to the scheduler (defaults to
            no limit)
          poll_interval (float): optional polling interval
            (seconds) to set in scheduler (defaults to 5s)
          runners (dict): mapping of names to JobRunner
            instances; valid names are 'qc_runner',
            'report_runner','cellranger_runner',
            'verify_runner','default'
          envmodules (mapping): mapping of names to
            environment module file lists; valid names are
            'bcl2fastq','cellranger_mkfastq',
            'cellranger_atac_mkfastq'
          default_runner (JobRunner): optional default
            job runner to use
          verbose (bool): if True then report additional
            information for diagnostics
        """
        # Working directory
        clean_up_on_completion = False
        if working_dir is None:
            working_dir = tempfile.mkdtemp(prefix="__make_fastqs.",
                                           suffix=".tmp",
                                           dir=os.getcwd())
            clean_up_on_completion = True
        working_dir = os.path.abspath(working_dir)
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)

        # Output directory
        if out_dir is None:
            out_dir = "bcl2fastq"
        if not os.path.isabs(out_dir):
            out_dir = os.path.join(analysis_dir,out_dir)

        # Primary data directory
        if primary_data_dir is None:
            primary_data_dir = "primary_data"
        if not os.path.isabs(primary_data_dir):
            primary_data_dir = os.path.join(analysis_dir,
                                            primary_data_dir)

        # Statistics files
        # Basic stats
        if stats_file is None:
            if name:
                stats_file = "statistics.%s.info" % name
            else:
                stats_file = "statistics.info"
        if not os.path.isabs(stats_file):
            stats_file = os.path.join(analysis_dir,stats_file)
        # Full stats
        if stats_full is None:
            if name:
                stats_full = "statistics_full.%s.info" % name
            else:
                stats_full = "statistics_full.info"
        if not os.path.isabs(stats_full):
            stats_full = os.path.join(analysis_dir,stats_full)
        # Per-lane stats
        if per_lane_stats is None:
            if name:
                per_lane_stats = "per_lane_statistics.%s.info" % name
            else:
                per_lane_stats = "per_lane_statistics.info"
        if not os.path.isabs(per_lane_stats):
            per_lane_stats = os.path.join(analysis_dir,per_lane_stats)
        # Per-lane per-sample stats
        if per_lane_sample_stats is None:
            if name:
                per_lane_sample_stats = "per_lane_sample_stats.%s.info" % name
            else:
                per_lane_sample_stats = "per_lane_sample_stats.info"
        if not os.path.isabs(per_lane_sample_stats):
            per_lane_sample_stats = os.path.join(analysis_dir,
                                                 per_lane_sample_stats)

        # QC report
        if name:
            qc_report = "processing_qc_%s.html" % name
        else:
            qc_report = "processing_qc.html"
        qc_report = os.path.join(analysis_dir,qc_report)

        # Barcode analysis directory
        if barcode_analysis_dir is None:
            if name:
                barcode_analysis_dir = "barcode_analysis_%s" % name
            else:
                barcode_analysis_dir = "barcode_analysis"
        if not os.path.isabs(barcode_analysis_dir):
            barcode_analysis_dir = os.path.join(analysis_dir,
                                                barcode_analysis_dir)

        # Log and script directories
        if log_dir is None:
            log_dir = os.path.join(working_dir,"logs")
        scripts_dir = os.path.join(working_dir,"scripts")

        # Runners
        if runners is None:
            runners = dict()

        # Parameters
        analysis_dir = os.path.abspath(analysis_dir)
        params = {
            'analysis_dir': analysis_dir,
            'out_dir': out_dir,
            'primary_data_dir': primary_data_dir,
            'barcode_analysis_dir': barcode_analysis_dir,
            'counts_dir': os.path.join(barcode_analysis_dir,"counts"),
            'qc_report': qc_report,
            'force_copy_of_primary_data': force_copy_of_primary_data,
            'no_lane_splitting': no_lane_splitting,
            'create_fastq_for_index_read': create_fastq_for_index_read,
            'find_adapters_with_sliding_window':
            find_adapters_with_sliding_window,
            'create_empty_fastqs': create_empty_fastqs,
            'name': name,
            'stats_file': stats_file,
            'stats_full': stats_full,
            'per_lane_stats': per_lane_stats,
            'per_lane_sample_stats': per_lane_sample_stats,
            'nprocessors': nprocessors,
            'require_bcl2fastq': require_bcl2fastq,
            'require_bclconvert': require_bclconvert,
            'cellranger_jobmode': cellranger_jobmode,
            'cellranger_mempercore': cellranger_mempercore,
            'cellranger_maxjobs': cellranger_maxjobs,
            'cellranger_jobinterval': cellranger_jobinterval,
            'cellranger_localcores': cellranger_localcores,
            'cellranger_localmem': cellranger_localmem,
        }

        # Execute the pipeline
        status = Pipeline.run(self,
                              working_dir=working_dir,
                              log_dir=log_dir,
                              scripts_dir=scripts_dir,
                              log_file=log_file,
                              batch_size=batch_size,
                              batch_limit=batch_limit,
                              exit_on_failure=PipelineFailure.DEFERRED,
                              params=params,
                              poll_interval=poll_interval,
                              max_jobs=max_jobs,
                              max_slots=max_slots,
                              runners=runners,
                              default_runner=default_runner,
                              envmodules=envmodules,
                              verbose=verbose)

        # Clean up working dir
        if status == 0 and clean_up_on_completion:
            shutil.rmtree(working_dir)

        # Return pipeline status
        return status

######################################################################
# Pipeline task classes
######################################################################

class FetchPrimaryData(PipelineTask):
    """
    Fetch the primary data for processing
    """
    def init(self,data_dir,primary_data_dir,force_copy=False):
        """
        Initialise the FetchPrimaryData task

        Arguments:
          data_dir (str): location of the source sequencing data
          primary_data_dir (str): directory to copy data to (if
            source is a remote location) or link data from (if
            source is on the local system)
          force_copy (bool): if True then force primary data to
            be copied even if it's on the local system

        Outputs:
          run_dir: path to the local copy of the primary data
        """
        self.tmp_dir = None
        self.add_output('run_dir',Param(type='str'))
    def setup(self):
        # Check if primary data already exists
        final_run_dir = os.path.join(
            self.args.primary_data_dir,
            os.path.basename(self.args.data_dir))
        if os.path.exists(final_run_dir):
            print("Primary data already acquired")
            return
        # Create top level primary data directory
        if not os.path.exists(self.args.primary_data_dir):
            print("Making top-level primary data dir: %s" %
                  self.args.primary_data_dir)
            os.mkdir(self.args.primary_data_dir)
        # Check if source data are local or remote
        if not Location(self.args.data_dir).is_remote:
            # Local data
            print("Source data are on the local system")
            if not self.args.force_copy:
                # Make a symlink
                print("Making a symlink to source data")
                os.symlink(self.args.data_dir,
                           os.path.join
                           (self.args.primary_data_dir,
                            os.path.basename(self.args.data_dir)))
                return
        # Copy remote data (and local data when 'force copy'
        # was specified)
        print("Copying source data")
        # Make a temporary directory for rsyncing
        self.tmp_dir = os.path.join(self.args.primary_data_dir,
                                    "__%s.part" %
                                    os.path.basename(self.args.data_dir))
        # Construct the rsync command
        rsync_cmd = general_apps.rsync(
            "%s/" % self.args.data_dir,
            self.tmp_dir,
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
        print("Running %s" % rsync_cmd)
        self.add_cmd(PipelineCommandWrapper(
            "Rsync primary data from %s" % self.args.data_dir,
            *rsync_cmd.command_line))
    def finish(self):
        final_run_dir = os.path.join(
            self.args.primary_data_dir,
            os.path.basename(self.args.data_dir))
        # If copying then move the tmp copy to the
        # final location
        if self.tmp_dir is not None:
            print("Moving temporary copy of primary data to final "
                  "location")
            os.rename(self.tmp_dir,final_run_dir)
        if os.path.exists(final_run_dir):
            self.output.run_dir.set(final_run_dir)

class IdentifyPlatform(PipelineTask):
    """
    Identify the sequencer platform from the primary data
    """
    def init(self,run_dir,platform=None):
        """
        Initialise the IdentifyPlatform task

        Arguments:
          run_dir (str): path to the sequencer run data
          platform (str): optional, specify the platform
        Outputs:
          platform: sequencer platform
        """
        self.add_output('platform',Param(type='str'))
    def setup(self):
        # Load input data and acquire the platform
        illumina_run = IlluminaRun(self.args.run_dir,
                                   platform=self.args.platform)
        self.output.platform.set(illumina_run.platform)
        print("Platform identified as '%s'" % illumina_run.platform)

class MakeSampleSheet(PipelineTask):
    """
    Creates a custom sample sheet
    """
    def init(self,sample_sheet_file,lanes=(),adapter=None,
             adapter_read2=None):
        """
        Initialise the MakeSampleSheet task

        Arguments:
          sample_sheet_file (str): name and path of the base
            sample file to generate the new file from
          lanes (list): (optional) list of lane numbers to
            keep in the output sample sheet; if empty then all
            lanes will be kept
          adapter (str): (optional) if set then write to the
            `Adapter` setting
          adapter_read2 (str): (optional) if set then write to
            the `AdapterRead2` setting

        Outputs:
          custom_sample_sheet (PipelineParam): pipeline
            parameter instance that resolves to a string
            with the path to the output sample sheet file.
        """
        self.add_output('custom_sample_sheet',Param(type=str))
    def setup(self):
        print("Generating custom sample sheet from %s" %
              self.args.sample_sheet_file)
        # Construct a name for the output file
        if self.args.lanes:
            lanes_id = ".L%s" % ''.join([str(l) for l in self.args.lanes])
        else:
            lanes_id = ""
        output_sample_sheet = os.path.abspath("SampleSheet%s.%s.csv" %
                                              (lanes_id,
                                               time.strftime("%Y%m%d%H%M%S")))
        # Generate the custom sample sheet
        make_custom_sample_sheet(self.args.sample_sheet_file,
                                 output_sample_sheet,
                                 lanes=(None if not self.args.lanes
                                        else self.args.lanes),
                                 adapter=self.args.adapter,
                                 adapter_read2=self.args.adapter_read2)
        # Check the temporary sample sheet
        print("Checking custom sample sheet %s" % output_sample_sheet)
        linter = SampleSheetLinter(sample_sheet_file=output_sample_sheet)
        # Invalid barcodes
        invalid_barcodes = linter.has_invalid_barcodes()
        if invalid_barcodes:
            print("Invalid barcodes detected:")
            for line in invalid_barcodes:
                print("!!! %s" % line)
        # Invalid characters
        invalid_characters = linter.has_invalid_characters()
        if invalid_characters:
            print("Invalid non-printing/non-ASCII characters "
                  "detected")
        if invalid_barcodes or invalid_characters:
            raise Exception("Errors detected in generated sample sheet")
        # Sample sheet ok
        self.output.custom_sample_sheet.set(output_sample_sheet)

class GetBcl2Fastq(PipelineFunctionTask):
    """
    Get information on the bcl2fastq executable
    """
    def init(self,require_bcl2fastq=None):
        """
        Initialise the GetBcl2Fastq task

        Arguments:
          require_bcl2fastq (str): if set then should be a
            string of the form '1.8.4' or '>2.0', explicitly
            specifying the version of bcl2fastq to use. If
            not set then no version check will be made

        Outputs:
          bcl2fastq_exe (str): path to the bcl2fastq executable
          bcl2fastq_package (str): name of the bcl2fastq package
          bcl2fastq_version (str): the bcl2fastq version
          bcl2fastq_info (tuple): tuple consisting of
            (exe,package,version)
        """
        self.add_output('bcl2fastq_exe',Param(type=str))
        self.add_output('bcl2fastq_package',Param(type=str))
        self.add_output('bcl2fastq_version',Param(type=str))
        self.add_output('bcl2fastq_info',Param())
    def setup(self):
        if self.args.require_bcl2fastq:
            print("Requires bcl2fastq version %s" %
                  self.args.require_bcl2fastq)
        self.add_call("Check bcl2fastq version",
                      self.get_bcl2fastq,
                      self.args.require_bcl2fastq)
    def get_bcl2fastq(self,require_bcl2fastq=None):
        # Get bcl2fastq
        bcl2fastq = available_bcl2fastq_versions(require_bcl2fastq)
        if bcl2fastq:
            bcl2fastq_exe = bcl2fastq[0]
            bcl2fastq_info = bcl_to_fastq_info(bcl2fastq_exe)
        else:
            msg = "No appropriate bcl2fastq software located"
            if not require_bcl2fastq:
                msg += " matching requirement '%s'" % \
                       require_bcl2fastq
            raise Exception(msg)
        # Return the information on bcl2fastq
        return (bcl2fastq_exe,bcl2fastq_info[1],bcl2fastq_info[2])
    def finish(self):
        bcl2fastq = self.result()[0]
        bcl2fastq_exe = bcl2fastq[0]
        bcl2fastq_package = bcl2fastq[1]
        bcl2fastq_version = bcl2fastq[2]
        # Set outputs with info on bcl2fastq executable
        self.output.bcl2fastq_exe.set(bcl2fastq_exe)
        self.output.bcl2fastq_package.set(bcl2fastq_package)
        self.output.bcl2fastq_version.set(bcl2fastq_version)
        self.output.bcl2fastq_info.set((bcl2fastq_exe,
                                        bcl2fastq_package,
                                        bcl2fastq_version))
        # Report what was found
        print("Bcl2fastq exe    : %s" % bcl2fastq_exe)
        print("Bcl2fastq package: %s" % bcl2fastq_package)
        print("Bcl2fastq version: %s" % bcl2fastq_version)

class GetBclConvert(PipelineFunctionTask):
    """
    Get information on the bcl-convert executable
    """
    def init(self,require_version=None):
        """
        Initialise the GetBcl2Fastq task

        Arguments:
          require_version (str): if set then should be a
            string of the form '1.8.4' or '>2.0', explicitly
            specifying the version of bcl-convert to use. If
            not set then no version check will be made

        Outputs:
          bclconvert_exe (str): path to the bcl-convert executable
          bclconvert_package (str): name of the bcl-convert package
          bclconvert_version (str): the bcl-convert version
          bclconvert_info (tuple): tuple consisting of
            (exe,package,version)
        """
        self.add_output('bclconvert_exe',Param(type=str))
        self.add_output('bclconvert_package',Param(type=str))
        self.add_output('bclconvert_version',Param(type=str))
        self.add_output('bclconvert_info',Param())
    def setup(self):
        if self.args.require_version:
            print("Requires bcl-convert version %s" %
                  self.args.require_version)
        self.add_call("Check bcl2fastq version",
                      self.get_bclconvert,
                      self.args.require_version)
    def get_bclconvert(self,require_version=None):
        # Get bcl-convert
        bclconvert_exe = find_program("bcl-convert")
        if bclconvert_exe:
            # Get information on the version etc
            bclconvert_info_ = bclconvert_info(bclconvert_exe)
        else:
            msg = "No appropriate BCL Convert software located"
            if require_version:
                msg += " matching requirement '%s'" % require_version
            raise Exception(msg)
        # Return the information on bcl2fastq
        return (bclconvert_exe,bclconvert_info_[1],bclconvert_info_[2])
    def finish(self):
        bclconvert = self.result()[0]
        bclconvert_exe = bclconvert[0]
        bclconvert_package = bclconvert[1]
        bclconvert_version = bclconvert[2]
        # Set outputs with info on bcl-convert executable
        self.output.bclconvert_exe.set(bclconvert_exe)
        self.output.bclconvert_package.set(bclconvert_package)
        self.output.bclconvert_version.set(bclconvert_version)
        self.output.bclconvert_info.set((bclconvert_exe,
                                         bclconvert_package,
                                         bclconvert_version))
        # Report what was found
        print("BCL Convert exe    : %s" % bclconvert_exe)
        print("BCL Convert package: %s" % bclconvert_package)
        print("BCL Convert version: %s" % bclconvert_version)

class RestoreBackupDirectory(PipelineTask):
    """
    Check for and restore saved copy of directory

    Looks for a backup version of a directory, and
    restores it by renaming it back to the original
    name if found.

    Back up for directory `/path/to/dir` will be
    called `/path/to/save.dir`.
    """
    def init(self,dirn,skip_restore=False):
        """
        Initialise the RestoreBackupDirectory task

        Arguments:
          dirn (str): path to the original directory
            to look for backup of
          skip_restore (bool): if True then check for
            the backup but don't restore it if found
        """
        pass
    def setup(self):
        # Get path to for backup copy
        backup_dirn = os.path.join(
            os.path.dirname(self.args.dirn),
            "save.%s" % os.path.basename(self.args.dirn)
        )
        print("Checking for backup directory: '%s'" % backup_dirn)
        if os.path.exists(backup_dirn):
            if not self.args.skip_restore:
                print("Restoring to '%s'" % self.args.dirn)
                os.rename(backup_dirn,self.args.dirn)
            else:
                print("Backup located but restore skipped")
        else:
            print("Not found, skipped restore")

class RunBcl2Fastq(PipelineTask):
    """
    Run bcl2fastq to generate Fastqs from sequencing data
    """
    def init(self,run_dir,out_dir,sample_sheet,bases_mask='auto',
             ignore_missing_bcl=False,no_lane_splitting=False,
             minimum_trimmed_read_length=None,
             mask_short_adapter_reads=None,
             create_fastq_for_index_read=False,
             find_adapters_with_sliding_window=False,nprocessors=None,
             create_empty_fastqs=False,
             platform=None,bcl2fastq_exe=None,bcl2fastq_version=None,
             skip_bcl2fastq=False):
        """
        Initialise the RunBcl2Fastq task

        Arguments:
          run_dir (str): path to the source sequencing data
          out_dir (str): output directory for bcl2fastq
          sample_sheet (str): path to input samplesheet file
          bases_mask (str): if set then use this as an
            alternative bases mask setting
          ignore_missing_bcl (bool): if True then run
            bcl2fastq with --ignore-missing-bcl
          no_lane_splitting (bool): if True then run bcl2fastq
            with --no-lane-splitting
          minimum_trimmed_read_length (int): if set then supply
            to bcl2fastq via --minimum-trimmed-read-length
          mask_short_adapter_reads (int): if set then supply to
            bcl2fastq via --mask-short-adapter-reads
          create_fastq_for_index_read (boolean): if True then
            also create Fastq files for index reads (default,
            don't create index read Fastqs)
          find_adapters_with_sliding_window (bool): if True
            then use sliding window algorith for identifying
            adapter sequences (default is to use string
            matching algorithm)
          nprocessors (int): number of processors to use
            (taken from job runner by default)
          create_empty_fastqs (bool): if True then create empty
            placeholder Fastq files for any that are missing
            on successful completion of bcl2fastq
          platform (str): optional, sequencing platform that
            generated the data
          bcl2fastq_exe (str): the path to the bcl2fastq
            executable to use
          bcl2fastq_version (str): the version string for the
            bcl2fastq package
          skip_bcl2fastq (bool): if True then sets the output
            parameters but finishes before actually running
            bcl2fastq

        Outputs:
          bases_mask: actual bases mask used
          mismatches: number of mismatches allowed
          missing_fastqs: list of Fastqs missing after
            Fastq generation
        """
        # Internal variables
        self.supported_versions = ('1.8','2.17','2.20',)
        self.tmp_out_dir = None
        # Outputs
        self.add_output('bases_mask',Param(type='str'))
        self.add_output('mismatches',Param(type='int'))
        self.add_output('missing_fastqs',list())
    def setup(self):
        # Load input data
        illumina_run = IlluminaRun(self.args.run_dir,
                                   platform=self.args.platform)
        # Set bases mask
        if self.args.bases_mask == "auto":
            print("Setting bases mask from RunInfo.xml")
            bases_mask = get_bases_mask(illumina_run.runinfo_xml,
                                        self.args.sample_sheet)
        else:
            bases_mask = self.args.bases_mask
        if not bases_mask_is_valid(bases_mask):
            raise Exception("Invalid bases mask: '%s'" %
                            bases_mask)
        self.output.bases_mask.set(bases_mask)
        # Check sample sheet for collisions and set mismatches
        mismatches = get_nmismatches(bases_mask)
        while mismatches >= 0:
            if check_barcode_collisions(self.args.sample_sheet,
                                        mismatches):
                mismatches -= 1
            else:
                break
        if mismatches < 0:
            raise Exception("Barcode collisions detected even with "
                            "zero mismatches (indicates duplicated "
                            "indexes?)")
        self.output.mismatches.set(mismatches)
        # Check if Fastq generation should be skipped
        if self.args.skip_bcl2fastq:
            print("Skipping bcl2fastq run")
            return
        # Check if outputs already exist
        if os.path.exists(self.args.out_dir):
            print("Output directory %s already exists" %
                  self.args.out_dir)
            # Verify outputs
            illumina_data = IlluminaData(os.path.dirname(self.args.out_dir),
                                         os.path.basename(self.args.out_dir))
            if verify_run_against_sample_sheet(illumina_data,
                                               self.args.sample_sheet):
                print("Verified existing outputs against samplesheet")
                return
            else:
                raise Exception("Failed to verify existing outputs "
                                "against samplesheet")
        # Check that bcl2fastq version is supported
        supported_version = None
        for version in self.supported_versions:
            if self.args.bcl2fastq_version.startswith("%s." % version):
                supported_version = version
                break
        if supported_version:
            print("Bcl2fastq software matches supported version '%s'" %
                  supported_version)
        else:
            raise Exception("Don't know how to run bcl2fastq "
                            "version %s" % self.args.bcl2fastq_version)
        # Set up parameters
        params = {
            'mismatches': mismatches,
            'bases_mask': bases_mask,
            'ignore_missing_bcl': self.args.ignore_missing_bcl,
            'no_lane_splitting': self.args.no_lane_splitting,
            'minimum_trimmed_read_length':
            self.args.minimum_trimmed_read_length,
            'mask_short_adapter_reads':
            self.args.mask_short_adapter_reads,
            'create_fastq_for_index_reads':
            self.args.create_fastq_for_index_read,
            'find_adapters_with_sliding_window':
            self.args.find_adapters_with_sliding_window,
            'loading_threads': None,
            'demultiplexing_threads': None,
            'processing_threads': None,
            'writing_threads': None,
            'bcl2fastq_exe': self.args.bcl2fastq_exe,
        }
        # Update parameters based on bcl2fastq version
        if self.args.nprocessors:
            nprocessors = self.args.nprocessors
        else:
            nprocessors = self.runner_nslots
        if nprocessors is not None:
            # Explicitly set number of threads for each stage
            print("Nprocessors: %s" % nprocessors)
            if supported_version in ('2.17',):
                # bcl2fastq 2.17.*
                params['loading_threads'] = min(4,nprocessors)
                params['writing_threads'] = min(4,nprocessors)
                params['demultiplexing_threads'] = max(
                    int(float(nprocessors)*0.2),
                    nprocessors)
                params['processing_threads'] = nprocessors
            elif supported_version in ('2.20',):
                # bcl2fastq 2.20.*
                params['loading_threads'] = min(4,nprocessors)
                params['writing_threads'] = min(4,nprocessors)
                params['processing_threads'] = nprocessors
        # Report settings
        print("%-22s: %s" % ("Run dir",self.args.run_dir))
        print("%-22s: %s" % ("Sample sheet",
                             os.path.basename(self.args.sample_sheet)))
        print("%-22s: %s" % ("Output dir",self.args.out_dir))
        print("%-22s: %s" % ("Create empty Fastqs",self.args.create_empty_fastqs))
        for item,desc in (('bases_mask',"Bases mask"),
                          ('mismatches',"Allowed mismatches",),
                          ('ignore_missing_bcl',"Ignore missing bcl"),
                          ('no_lane_splitting',"No lane splitting"),
                          ('minimum_trimmed_read_length',"Min trimmed read len"),
                          ('mask_short_adapter_reads',"Mask short adptr reads"),
                          ('create_fastq_for_index_read',"Create index Fastqs"),
                          ('find_adapters_with_sliding_window',
                           "Use sliding window")):
            if item in params:
                print("%-22s: %s" % (desc,params[item]))
        print("Threads for each stage:")
        for item,desc in (('nprocessors',"Nprocessors"),
                          ('loading_threads',"Loading (-r)"),
                          ('demultiplexing_threads',"Demultiplexing (-d)"),
                          ('processing_threads',"Processing (-p)"),
                          ('writing_threads',"Writing (-w)")):
            if item in params and params[item] is not None:
                print("- %-20s: %s" % (desc,params[item]))
        # Set up temporary output dir
        self.tmp_out_dir = os.path.abspath("__%s.work" %
                                           os.path.basename(
                                               self.args.out_dir))
        # Build command to run bcl2fastq
        bcl2fastq2_cmd = bcl2fastq_apps.bcl2fastq2(
            self.args.run_dir,
            self.args.sample_sheet,
            output_dir=self.tmp_out_dir,
            **params)
        print("Running %s" % bcl2fastq2_cmd)
        self.add_cmd(PipelineCommandWrapper(
            "Run bcl2fastq",
            *bcl2fastq2_cmd.command_line))
    def finish(self):
        if self.tmp_out_dir:
            # Verify outputs and get list of missing Fastqs (if any)
            missing_fastqs = verify_run(self.tmp_out_dir,
                                        self.args.sample_sheet)
            # Update output
            for fq in missing_fastqs:
                self.output.missing_fastqs.append(fq)
            # Handle missing fastqs
            if missing_fastqs:
                if self.args.create_empty_fastqs:
                    # Create empty placeholder Fastqs
                    create_placeholder_fastqs(missing_fastqs,
                                              base_dir=self.tmp_out_dir)
                else:
                    # Terminate with an exception
                    raise Exception("Failed to verify outputs against "
                                    "samplesheet")
            # Move to final location
            print("Moving output to final location: %s" % self.args.out_dir)
            os.rename(self.tmp_out_dir,self.args.out_dir)

class RunBclConvert(PipelineTask):
    """
    Run BCL Convert to generate Fastqs from sequencing data
    """
    def init(self,run_dir,out_dir,sample_sheet,lane=None,
             bases_mask='auto',ignore_missing_bcl=False,
             no_lane_splitting=False,
             minimum_trimmed_read_length=None,
             mask_short_adapter_reads=None,
             create_fastq_for_index_read=False,nprocessors=None,
             create_empty_fastqs=False,ignore_missing_fastqs=False,
             platform=None,bclconvert_exe=None,bclconvert_version=None,
             skip_bclconvert=False):
        """
        Initialise the RunBclConvert task

        Arguments:
          run_dir (str): path to the source sequencing data
          out_dir (str): output directory for bcl2fastq
          sample_sheet (str): path to input samplesheet file
          lane (int): optional, run bcl-convert on a single lane
            with --bcl-only-lane
          bases_mask (str): if set then use this as an
            alternative bases mask setting
          no_lane_splitting (bool): if True then run bcl-convert
            with --no-lane-splitting
          minimum_trimmed_read_length (int): if set then supply
            to bcl-convert via sample sheet settings
          mask_short_adapter_reads (int): if set then supply to
            bcl-convert via sample sheet settings
          create_fastq_for_index_read (boolean): if True then
            also create Fastq files for index reads (default,
            don't create index read Fastqs)
          nprocessors (int): number of processors to use
            (taken from job runner by default)
          create_empty_fastqs (bool): if True then create empty
            placeholder Fastq files for any that are missing
            on successful completion of bcl-convert
          ignore_missing_fastqs (bool): if True then ignore
            missing Fastqs on successful completion of
            bcl-convert
          platform (str): optional, sequencing platform that
            generated the data
          bclconvert_exe (str): the path to the bcl-convert
            executable to use
          bclconvert_version (str): the version string for the
            bcl-convert package
          skip_bclconvert (bool): if True then sets the output
            parameters but finishes before actually running
            bcl-convert

        Outputs:
          bases_mask: actual bases mask used
          mismatches: number of mismatches allowed
          missing_fastqs: list of Fastqs missing after
            Fastq generation
        """
        # Internal variables
        self.supported_versions = ('3.7',)
        self.tmp_out_dir = None
        self.bclconvert_sample_sheet = None
        # Outputs
        self.add_output('bases_mask',Param(type='str'))
        self.add_output('mismatches',Param(type='int'))
        self.add_output('missing_fastqs',list())
    def setup(self):
        # Load input data
        illumina_run = IlluminaRun(self.args.run_dir,
                                   platform=self.args.platform)
        # Set bases mask
        if self.args.bases_mask == "auto":
            print("Setting bases mask from RunInfo.xml")
            bases_mask = get_bases_mask(illumina_run.runinfo_xml,
                                        self.args.sample_sheet)
        else:
            bases_mask = self.args.bases_mask
        if not bases_mask_is_valid(bases_mask):
            raise Exception("Invalid bases mask: '%s'" %
                            bases_mask)
        self.output.bases_mask.set(bases_mask)
        # Check sample sheet for collisions and set mismatches
        # NB BCL Convert sets one mismatch value per index
        mismatches = get_nmismatches(bases_mask,multi_index=True)
        updated_mismatches = []
        for ix,nmismatches in enumerate(mismatches,start=1):
            while nmismatches >= 0:
                if check_barcode_collisions(self.args.sample_sheet,
                                            nmismatches,use_index=ix):
                    nmismatches -= 1
                else:
                    break
            if nmismatches < 0:
                print("Warning: barcode collisions detected even with "
                      "zero mismatches for index %s" % ix)
            updated_mismatches.append(nmismatches)
        if all([(m < 0) for m in updated_mismatches]):
            # All index sequences are identical
            raise Exception("Barcode collisions detected across all "
                            "index sequences even with zero mismatches")
        else:
            # Reset the 'bad' mismatches to 1
            updated_mismatches = [(m if m >= 0 else 1)
                                  for m in updated_mismatches]
        mismatches = updated_mismatches
        self.output.mismatches.set(mismatches)
        # Check if Fastq generation should be skipped
        if self.args.skip_bclconvert:
            print("Skipping BCL Convert run")
            return
        # Check that BCL Convert version is supported
        supported_version = None
        for version in self.supported_versions:
            if self.args.bclconvert_version.startswith("%s." % version):
                supported_version = version
                break
        if supported_version:
            print("BCL Convert software matches supported version '%s'" %
                  supported_version)
        else:
            raise Exception("Don't know how to run BCL Convert "
                            "version %s" % self.args.bclconvert_version)
        # Write an updated sample sheet with the required settings
        sample_sheet = SampleSheet(self.args.sample_sheet)
        if bases_mask:
            # Convert bases mask to OverrideCycles format
            # (convert delimiters from comma to semi-colon)
            override_cycles = bases_mask.replace(',',';')
            sample_sheet.settings['OverrideCycles'] = override_cycles
        if mismatches:
            # Mismatches are set independently for index1 and index2
            for ix,nmismatches in enumerate(mismatches):
                sample_sheet.settings['BarcodeMismatchesIndex%d' % (ix+1)] = \
                    mismatches[ix]
        if self.args.minimum_trimmed_read_length is not None:
            sample_sheet.settings['MinimumTrimmedReadLength'] = \
                self.args.minimum_trimmed_read_length
        if self.args.mask_short_adapter_reads is not None:
            sample_sheet.settings['MaskShortReads'] = \
                self.args.mask_short_adapter_reads
        if self.args.create_fastq_for_index_read:
            sample_sheet.settings['CreateFastqForIndexReads'] = 1
        else:
            sample_sheet.settings['CreateFastqForIndexReads'] = 0
        # Correct adapter setting settings for BCL Convert
        if 'Adapter' in sample_sheet.settings:
            # Update to 'AdapterRead1'
            adapter_seq = sample_sheet.settings['Adapter']
            sample_sheet.settings['AdapterRead1'] = adapter_seq
            del(sample_sheet.settings['Adapter'])
        # Remove unwanted lanes
        if self.args.lane is not None:
            i = 0
            while i < len(sample_sheet):
                if sample_sheet[i]['Lane'] != self.args.lane:
                    # Drop this sample/lane
                    del(sample_sheet[i])
                else:
                    # Keep this sample/lane
                    i += 1
        self.bclconvert_sample_sheet = os.path.join(
            os.path.dirname(self.args.sample_sheet),
            "SampleSheet.bclconvert%s.%s.csv" %
            ('L%s' % self.args.lane if self.args.lane is not None else '',
            time.strftime("%Y%m%d%H%M%S")))
        sample_sheet.write(self.bclconvert_sample_sheet)
        # Check if outputs already exist
        if os.path.exists(self.args.out_dir):
            print("Output directory %s already exists" %
                  self.args.out_dir)
            # Verify outputs
            illumina_data = IlluminaData(os.path.dirname(self.args.out_dir),
                                         os.path.basename(self.args.out_dir))
            if verify_run_against_sample_sheet(illumina_data,
                                               self.bclconvert_sample_sheet):
                print("Verified existing outputs against samplesheet")
                return
            else:
                raise Exception("Failed to verify existing outputs "
                                "against samplesheet")
        # Set up parameters
        params = {
            'no_lane_splitting': self.args.no_lane_splitting,
            'sampleproject_subdirectories': True,
            'num_parallel_tiles': None,
            'num_conversion_threads': None,
            'num_compression_threads': None,
            'num_decompression_threads': None,
            'bclconvert_exe': self.args.bclconvert_exe,
        }
        # Add lane
        if self.args.lane:
            params['lane'] = self.args.lane
        # Update parameters based on bcl2fastq version
        if self.args.nprocessors:
            nprocessors = self.args.nprocessors
        else:
            nprocessors = self.runner_nslots
        if nprocessors is not None:
            # Explicitly set number of threads for each stage
            print("Nprocessors: %s" % nprocessors)
            params['num_parallel_tiles'] = 1
            params['num_conversion_threads'] = min(8,nprocessors)
            params['num_compression_threads'] = min(8,nprocessors)
            params['num_decompression_threads'] = \
                max(int(float(nprocessors)/2.0),1)
        # Report settings
        print("%-22s: %s" % ("Input Bcl dir",self.args.run_dir))
        print("%-22s: %s" % ("Sample sheet",
                             os.path.basename(self.bclconvert_sample_sheet)))
        print("%-22s: %s" % ("Output dir",self.args.out_dir))
        if self.args.lane:
            print("%-22s: %s" % ("Lane",self.args.lane))
        print("%-22s: %s" % ("Create empty Fastqs",
                             self.args.create_empty_fastqs))
        for item,desc in (('no_lane_splitting',"No lane splitting"),
                          ('sampleproject_subdirectories',
                           "Create SampleProject subdirs in output"),):
            if item in params:
                print("%-22s: %s" % (desc,params[item]))
        print("Threads for each stage:")
        for item,desc in (('nprocessors',"Nprocessors"),
                          ('num_parallel_tiles',
                           "Number of parallel tiles"),
                          ('num_conversion_threads',
                           "Threads per tile for Fastq conversion"),
                          ('num_compression_threads',
                           "Threads for compressing output Fastqs"),
                          ('num_decompression_threads',
                           "Threads for decompressing input BCLs")):
            if item in params and params[item] is not None:
                print("- %-38s: %s" % (desc,params[item]))
        print("SampleSheet settings:")
        for item in sample_sheet.settings_items:
            print("%-22s: %s" % (item,sample_sheet.settings[item]))
        # Set up temporary output dir
        self.tmp_out_dir = os.path.abspath("__%s.work" %
                                           os.path.basename(
                                               self.args.out_dir))
        # Build command to run bcl-convert
        bclconvert_cmd = bcl2fastq_apps.bclconvert(
            self.args.run_dir,
            self.tmp_out_dir,
            sample_sheet=self.bclconvert_sample_sheet,
            **params)
        print("Running %s" % bclconvert_cmd)
        self.add_cmd(PipelineCommandWrapper(
            "Run BCL Convert",
            *bclconvert_cmd.command_line))
    def finish(self):
        if self.tmp_out_dir:
            # Verify outputs and get list of missing Fastqs (if any)
            missing_fastqs = verify_run(self.tmp_out_dir,
                                        self.bclconvert_sample_sheet)
            # Update output
            for fq in missing_fastqs:
                self.output.missing_fastqs.append(fq)
            # Handle missing fastqs
            if missing_fastqs:
                if self.args.create_empty_fastqs:
                    # Create empty placeholder Fastqs
                    create_placeholder_fastqs(missing_fastqs,
                                              base_dir=self.tmp_out_dir)
                elif not self.args.ignore_missing_fastqs:
                    # Terminate with an exception
                    raise Exception("Failed to verify outputs against "
                                    "samplesheet")
            # Move to final location
            print("Moving output to final location: %s" % self.args.out_dir)
            os.rename(self.tmp_out_dir,self.args.out_dir)

class GetBasesMaskIcell8(PipelineTask):
    """
    Set the bases mask for ICELL8 RNA-seq data
    """
    def init(self,run_dir,sample_sheet):
        """
        Initialise the GetBasesMaskIcell8 task

        Arguments:
          run_dir (str): path to the directory with
            data from the sequencer run
          sample_sheet (str): path to the sample sheet
            file to be used for processing these data

        Outputs:
          bases_mask (str): bases mask to use in
            bcl2fastq for processing these data
        """
        self.add_output("bases_mask",Param(type='str'))
    def setup(self):
        # Reset the default bases mask
        bases_mask = IlluminaRunInfo(
            IlluminaRun(self.args.run_dir).runinfo_xml).bases_mask
        print("Initial bases_mask  : %s" % bases_mask)
        bases_mask = get_bases_mask_icell8(bases_mask,
                                           sample_sheet=self.args.sample_sheet)
        print("Corrected for ICELL8: %s" % bases_mask)
        self.output.bases_mask.set(bases_mask)

class GetBasesMaskIcell8Atac(PipelineTask):
    """
    Set the bases mask for ICELL8 ATAC-seq data
    """
    def init(self,run_dir):
        """
        Initialise the GetBasesMaskIcell8Atac task

        Arguments:
          run_dir (str): path to the directory with
            data from the sequencer run
          sample_sheet (str): path to the sample sheet
            file to be used for processing these data

        Outputs:
          bases_mask (str): bases mask to use in
            bcl2fastq for processing these data
        """
        self.add_output("bases_mask",Param(type='str'))
    def setup(self):
        # Get the bases mask
        bases_mask = get_bases_mask_icell8_atac(
            IlluminaRun(self.args.run_dir).runinfo_xml)
        print("Bases mask for ICELL8 ATAC: %s" % bases_mask)
        self.output.bases_mask.set(bases_mask)

class Get10xPackage(PipelineFunctionTask):
    """
    Get information on 10xGenomics software package
    """
    def init(self,require_package):
        """
        Initialise the Get10xPackage task

        If no matching package is located then the outputs
        are all set to 'None'.

        Arguments:
          require_package (str): name of the 10xGenomics
            package executable that is required (e.g. 'cellranger',
            'cellranger-atac')

        Outputs:
          package_name (str): name of the package
          package_exe (str): path to the package executable
          package_version (str): the package version
          package_info (tuple): tuple consisting of
            (exe,package,version)
        """
        self.add_output('package_name',Param(type=str))
        self.add_output('package_exe',Param(type=str))
        self.add_output('package_version',Param(type=str))
        self.add_output('package_info',Param())
    def setup(self):
        print("Look for '%s'" % self.args.require_package)
        self.add_call("Check %s" % self.args.require_package,
                      self.get_10x_package,
                      self.args.require_package)
    def get_10x_package(self,require_package):
        # Look for 10xGenomics package
        package_exe = find_program(require_package)
        package_name = os.path.basename(require_package)
        if package_exe:
            # Get information on the version etc
            if package_name in ('cellranger',
                                'cellranger-atac',
                                'cellranger-arc',):
                package_info = cellranger_info(package_exe)
            elif package_name == 'spaceranger':
                package_info = spaceranger_info(package_exe)
            else:
                raise Exception("%s: unknown 10xGenomics package" %
                                require_package)
        else:
            # Unable to locate appropriate software
            print("No appropriate %s software located" %
                  package_name)
            package_exe = None
            package_info = (None,None,None)
        # Return the information on the package
        return (package_exe,package_info[1],package_info[2])
    def finish(self):
        package = self.result()[0]
        package_exe = package[0]
        package_name = package[1]
        package_version = package[2]
        # Set outputs with info on package executable
        self.output.package_exe.set(package_exe)
        self.output.package_name.set(package_name)
        self.output.package_version.set(package_version)
        self.output.package_info.set((package_exe,
                                      package_name,
                                      package_version))
        # Report what was found
        print("Package exe    : %s" % package_exe)
        print("Package name   : %s" % package_name)
        print("Package version: %s" % package_version)

class DemultiplexIcell8Atac(PipelineTask):
    """
    Runs 'demultiplex_icell8_atac.py' to generate Fastqs
    """
    def init(self,fastq_dir,out_dir,well_list,
             nprocessors=None,swap_i1_and_i2=False,
             reverse_complement=None,
             skip_demultiplex=False):
        """
        Initialise the DemultiplexIcell8Atac task

        Arguments:
          fastq_dir (str): path to directory with
            Fastq files to demultiplex
          out_dir (str): path to output directory
          well_list (str): path to well list file
            to use for demultiplexing samples
          swap_i1_and_i2 (bool): if True then
            swap the I1 and I2 indexes when
            demultiplexing
          reverse_complement (str): whether to
            reverse complement I1, I2, or both,
            when demultiplexing
          skip_demultiplex (bool): if True then
            skip running the demultiplexing
        """
        # Internal variables
        self.tmp_out_dirs = {}
        self.illumina_data = None
        # Outputs
        # NB 'missing_fastqs' is a dummy output
        # included for compatibility with other
        # tasks in the pipeline
        self.add_output('missing_fastqs',list())
    def setup(self):
        # Check if demultiplexing should be skipped
        if self.args.skip_demultiplex:
            print("Skipping demultiplexing")
            return
        # Check if outputs already exist
        if os.path.exists(self.args.out_dir):
            print("Output directory %s already exists" %
                  self.args.out_dir)
            return
        # Collect Fastqs from bcl2fastq
        self.illumina_data = IlluminaData(
            os.path.dirname(self.args.fastq_dir),
            os.path.basename(self.args.fastq_dir))
        fastqs = []
        for project in self.illumina_data.projects:
            for sample in project.samples:
                for fq in sample.fastq:
                    fastqs.append(os.path.join(sample.dirn,fq))
        if not fastqs:
            raise Exception("No Fastqs found")
        # Check well list
        if not self.args.well_list:
            raise Exception("No well list file supplied")
        # Number of processors
        if self.args.nprocessors:
            nprocessors = self.args.nprocessors
        else:
            nprocessors = self.runner_nslots
        # Report settings
        for desc,param in (("Fastq dir",self.args.fastq_dir),
                           ("Well list file",self.args.well_list),
                           ("Output dir",self.args.out_dir),
                           ("Nprocessors",nprocessors),
                           ("Swap I1/I2",self.args.swap_i1_and_i2),
                           ("Reverse complement",
                            self.args.reverse_complement)):
            print("%-22s: %s" % (desc,param))
        # Do demultiplexing into samples based on well list
        for fqs in group_fastqs_by_name(fastqs):
            # Get sample name and lane from Fastq name
            illumina_fq = IlluminaFastq(fqs[0])
            name = "%s%s" % (illumina_fq.sample_name,
                             (".L%03d" % illumina_fq.lane_number
                              if illumina_fq.lane_number else ""))
            # Set up temporary output dir
            tmp_out_dir = os.path.abspath("demultiplex_icell8_atac.%s.work"
                                          % name)
            # Build demultiplexer command
            demultiplex_cmd = Command('demultiplex_icell8_atac.py',
                                      '--mode=samples',
                                      '--output-dir',tmp_out_dir,
                                      '--update-read-headers',
                                      '-n',nprocessors)
            if self.args.swap_i1_and_i2:
                demultiplex_cmd.add_args('--swap-i1-i2')
            if self.args.reverse_complement:
                demultiplex_cmd.add_args('--reverse-complement=%s' %
                                         self.args.reverse_complement)
            demultiplex_cmd.add_args('--unassigned',
                                     "Unassigned-%s" %
                                     illumina_fq.sample_name)
            demultiplex_cmd.add_args(self.args.well_list)
            demultiplex_cmd.add_args(*fqs)
            print("Running %s" % demultiplex_cmd)
            self.add_cmd(PipelineCommandWrapper(
                "Demultiplex ICELL8 ATAC Fastqs",
                *demultiplex_cmd.command_line))
            # Store the temporary output dir
            self.tmp_out_dirs[name] = tmp_out_dir
    def finish(self):
        if self.tmp_out_dirs:
            # Build bcl2fastq-style output directory
            bcl2fastq_dir = os.path.abspath("%s.work" %
                                            os.path.basename(self.args.out_dir))
            print("Building output dir: %s" % bcl2fastq_dir)
            for d in ('Stats','Reports',):
                mkdirs(os.path.join(bcl2fastq_dir,d))
            project = self.illumina_data.projects[0]
            mkdirs(os.path.join(bcl2fastq_dir,project.name))
            # Loop over samples
            for name in self.tmp_out_dirs:
                # Directory with demultiplexed Fastqs
                tmp_out_dir = self.tmp_out_dirs[name]
                print("Collecting Fastqs from %s:" % tmp_out_dir)
                # Copy (hard link) fastqs
                for fq in [f for f in os.listdir(tmp_out_dir)
                           if f.endswith(".fastq.gz")]:
                    print("-- %s" % fq)
                    os.link(os.path.join(tmp_out_dir,fq),
                            os.path.join(bcl2fastq_dir,
                                         project.name,fq))
                # Copy the reports and JSON file
                print("Collecting reports from %s:" % tmp_out_dir)
                for f in ('icell8_atac_stats.xlsx',
                          'icell8_atac_stats.json'):
                    print("-- %s" % f)
                    if len(self.tmp_out_dirs) > 1:
                        ff = "%s.%s%s" % (os.path.splitext(f)[0],
                                          name,
                                          os.path.splitext(f)[1])
                    else:
                        ff = f
                    os.link(os.path.join(tmp_out_dir,f),
                            os.path.join(bcl2fastq_dir,'Reports',ff))
            # Copy "undetermined" Fastqs to top level
            print("Collecting 'undetermined' Fastqs from %s:" %
                  self.args.fastq_dir)
            for fq in [f for f in os.listdir(self.args.fastq_dir)
                       if f.startswith("Undetermined_S0_")]:
                print("-- %s" % fq)
                os.link(os.path.join(self.args.fastq_dir,fq),
                        os.path.join(bcl2fastq_dir,fq))
            # Move outputs to final location
            print("Moving output to final location: %s" % self.args.out_dir)
            os.rename(bcl2fastq_dir,self.args.out_dir)

class Run10xMkfastq(PipelineTask):
    """
    Runs 10xGenomics 'mkfastq' to generate Fastqs
    """
    def init(self,run_dir,out_dir,sample_sheet,bases_mask='auto',
             minimum_trimmed_read_length=None,
             mask_short_adapter_reads=None,
             jobmode='local',maxjobs=None,
             mempercore=None,jobinterval=None,
             localcores=None,localmem=None,
             create_empty_fastqs=False,platform=None,
             pkg_exe=None,pkg_version=None,
             bcl2fastq_exe=None,bcl2fastq_version=None,
             skip_mkfastq=False):
        """
        Initialise the Run10xMkfastq task

        Arguments:
          run_dir (str): path to the directory with
            data from the sequencer run
          out_dir (str): output directory for cellranger
          sample_sheet (str): path to input samplesheet file
          bases_mask (str): if set then use this as an
            alternative bases mask setting
          minimum_trimmed_read_length (int): if set then supply
            to cellranger via --minimum-trimmed-read-length
          mask_short_adapter_reads (int): if set then supply to
            cellranger via --mask-short-adapter-reads
          jobmode (str): jobmode to use for
            running cellranger
          maxjobs (int): maximum number of concurrent
            jobs for 10xGenomics mkfastq to run
          mempercore (int): amount of memory available
            per core (for jobmode other than 'local')
          jobinterval (int): time to pause inbetween
            starting 10xGenomics mkfastq jobs
          localcores (int): number of cores available
            to 10xGenomics mkfastq in jobmode 'local'
          localmem (int): amount of memory available
            to 10xGenomics mkfastq in jobmode 'local'
          create_empty_fastqs (bool): if True then create empty
            placeholder Fastq files for any that are missing
            on successful completion of 10xGenomics mkfastq
          platform (str): optional, sequencing platform that
            generated the data
          pkg_exe (str): the path to the 10xGenomics
            software package to use (e.g. 'cellranger',
            'cellranger-atac', 'spaceranger')
          pkg_version (str): the version string for the
            10xGenomics package
          bcl2fastq_exe (str): the path to the bcl2fastq
            executable to use
          bcl2fastq_version (str): the version string for the
            bcl2fastq package
          skip_mkfastq (bool): if True then skip running
            the 'mkfastq' step within the task

        Outputs:
            missing_fastqs: list of Fastqs missing after
              Fastq generation

        """
        # Internal variables
        self.pkg = None
        self.tmp_out_dir = None
        self.lanes = None
        self.mkfastq_out_dir = None
        self.mro_file = None
        self.expect_qc_summary_json = True
        # Outputs
        self.add_output('missing_fastqs',list())
    def setup(self):
        # 10xGenomics software package
        if not self.args.pkg_exe:
            raise Exception("No 10xGenomics executable provided")
        self.pkg = os.path.basename(self.args.pkg_exe)
        # Check if mkfastq should be skipped
        if self.args.skip_mkfastq:
            print("Skipping %s mkfastq" % self.pkg)
            return
        # Check if outputs already exist
        if os.path.exists(self.args.out_dir):
            print("Output directory %s already exists" %
                  self.args.out_dir)
            return
        # Load input data
        illumina_run = IlluminaRun(self.args.run_dir,
                                   platform=self.args.platform)
        # Deal with lanes
        sample_sheet = SampleSheet(self.args.sample_sheet)
        if sample_sheet.has_lanes:
            self.lanes = sorted(list(set([line['Lane']
                                          for line in sample_sheet])))
        else:
            self.lanes = None
        # Determine expected names for cellranger outputs
        if self.lanes is not None:
            lanes_suffix = "_%s" % ''.join([str(l) for l in self.lanes])
        else:
            lanes_suffix = ""
        self.mkfastq_out_dir = "%s%s" % (illumina_run.runinfo.flowcell,
                                         lanes_suffix)
        self.mro_file = "__%s.mro" % self.mkfastq_out_dir
        # Set bases mask
        if self.pkg == "cellranger-atac":
            # scATAC-seq
            bases_mask = self.args.bases_mask
            if bases_mask is None:
                bases_mask = 'auto'
            if bases_mask == 'auto':
                # Update bases mask to only use first 8 bases from
                # first index e.g. I8nnnnnnnn and convert second index
                # to read e.g. Y16
                print("Determining bases mask from RunInfo.xml")
                bases_mask = get_bases_mask_10x_atac(illumina_run.runinfo_xml)
                print("Bases mask: %s (updated for 10x scATAC-seq)" %
                      bases_mask)
                if not bases_mask_is_valid(bases_mask):
                    raise Exception("Invalid bases mask: '%s'" %
                                    bases_mask)
        elif self.pkg in ('cellranger',
                          'cellranger-arc',
                          'spaceranger'):
            # scRNA-seq
            if self.args.bases_mask == "auto":
                bases_mask = None
            else:
                bases_mask = self.args.bases_mask
        else:
            raise Exception("%s: unsupported 10xGenomics package" %
                            self.pkg)
        # Check if outputs already exist
        if os.path.exists(self.args.out_dir):
            print("Output directory %s already exists" %
                  self.args.out_dir)
            return
        # Set up parameters
        params = {
            'bases_mask': bases_mask,
            'lanes': (','.join([str(l) for l in self.lanes])
                      if self.lanes else 'all'),
            'minimum_trimmed_read_length':
            self.args.minimum_trimmed_read_length,
            'mask_short_adapter_reads':
            self.args.mask_short_adapter_reads,
            'cellranger_jobmode': self.args.jobmode,
            'cellranger_maxjobs': self.args.maxjobs,
            'cellranger_mempercore': self.args.mempercore,
            'cellranger_jobinterval': self.args.jobinterval,
            'cellranger_localcores': self.args.localcores,
            'cellranger_localmem': self.args.localmem,
            'cellranger_working_dir': self.args.out_dir,
            'cellranger_mro_file': self.mro_file
        }
        # Report settings
        print("%-22s: %s" % ("Run dir",self.args.run_dir))
        print("%-22s: %s" % ("Sample sheet",
                             os.path.basename(self.args.sample_sheet)))
        print("%-22s: %s" % ("Output dir",self.args.out_dir))
        for item,desc in (('bases_mask',"Bases mask"),
                          ('lanes',"Lanes"),
                          ('minimum_trimmed_read_length',"Min trimmed read len"),
                          ('mask_short_adapter_reads',"Mask short adptr reads"),
                          ('cellranger_jobmode',"Cellranger jobmode"),
                          ('cellranger_maxjobs',"Cellranger maxjobs"),
                          ('cellranger_mempercore',"Cellranger mempercore"),
                          ('cellranger_jobinterval',"Cellranger jobinterval"),
                          ('cellranger_localcores',"Cellranger localcores"),
                          ('cellranger_localmem',"Cellranger localmem"),
                          ('cellranger_working_dir',"Cellranger working dir"),
                          ('cellranger_mro_file',"Cellranger .mro file")):
            if item in params and params[item] is not None:
                print("%-22s: %s" % (desc,params[item]))
        # Set up temporary output dir
        self.tmp_out_dir = os.path.abspath("__%s.work" %
                                           os.path.basename(
                                               self.args.out_dir))
        # Build command to run 'mkfastq'
        mkfastq_cmd = Command(self.args.pkg_exe,
                              "mkfastq",
                              "--run",self.args.run_dir,
                              "--samplesheet",self.args.sample_sheet,
                              "--output-dir",self.tmp_out_dir)
        include_qc_argument = True
        pkg_major_version = int(self.args.pkg_version.split('.')[0])
        if self.pkg == "cellranger":
            if pkg_major_version >= 6:
                # --qc removed in cellranger 6.0.0
                include_qc_argument = False
                self.expect_qc_summary_json = False
        elif self.pkg == "cellranger-atac":
            if pkg_major_version >= 2:
                # --qc removed in cellranger-atac 2.0.0
                include_qc_argument = False
                self.expect_qc_summary_json = False
        elif self.pkg == "cellranger-arc":
            if pkg_major_version >= 2:
                # --qc removed in cellranger-atac 2.0.0
                include_qc_argument = False
                self.expect_qc_summary_json = False
        if include_qc_argument:
            mkfastq_cmd.add_args("--qc")
        if self.lanes:
            mkfastq_cmd.add_args("--lanes",
                                 ','.join([str(l) for l in self.lanes]))
        if bases_mask:
            mkfastq_cmd.add_args("--use-bases-mask=%s" % bases_mask)
        if self.args.minimum_trimmed_read_length:
            mkfastq_cmd.add_args('--minimum-trimmed-read-length',
                                 self.args.minimum_trimmed_read_length)
        if self.args.mask_short_adapter_reads:
            mkfastq_cmd.add_args('--mask-short-adapter-reads',
                                 self.args.mask_short_adapter_reads)
        add_cellranger_args(
            mkfastq_cmd,
            jobmode=self.args.jobmode,
            mempercore=self.args.mempercore,
            maxjobs=self.args.maxjobs,
            jobinterval=self.args.jobinterval,
            localcores=self.args.localcores,
            localmem=self.args.localmem
        )
        print("Running %s" % mkfastq_cmd)
        self.add_cmd(PipelineCommandWrapper("Run %s" % self.pkg,
                                            *mkfastq_cmd.command_line))
    def finish(self):
        if self.tmp_out_dir:
            # Verify outputs
            #
            # Note: some 10x packages force the addition of a 'sample
            # name' subdirectory level (e.g. cellranger, cellranger-atac),
            # some don't (e.g. spaceranger), and some do both (e.g.
            # cellranger-arc) depending on the type of data
            # So we check both and see if either matches the prediction
            illumina_data = IlluminaData(os.path.dirname(self.tmp_out_dir),
                                         os.path.basename(self.tmp_out_dir))
            for include_sample_dir in (True,False):
                verified = verify_run_against_sample_sheet(
                    illumina_data,
                    self.args.sample_sheet,
                    include_sample_dir=include_sample_dir)
                if verified:
                    break
            if verified:
                print("Verified outputs against samplesheet")
            else:
                # Verification failed
                print("Failed to verify outputs against samplesheet")
                # List the missing Fastq files
                missing_fastqs = list_missing_fastqs(illumina_data,
                                                     self.args.sample_sheet)
                print("Missing Fastqs:")
                for fq in missing_fastqs:
                    print("- %s" % fq)
                    self.output.missing_fastqs.append(fq)
                if self.args.create_empty_fastqs:
                    # Create empty placeholder Fastqs
                    print("Making empty placeholder Fastqs")
                    for fq in missing_fastqs:
                        fastq = os.path.join(self.tmp_out_dir,fq)
                        # Make intermediate directory if required
                        if not os.path.exists(os.path.dirname(fastq)):
                            os.mkdir(os.path.dirname(fastq))
                        # Make empty file
                        with gzip.GzipFile(filename=fastq,mode='wb') as fp:
                            fp.write(''.encode())
                else:
                    # Terminate with an exception
                    raise Exception("Failed to verify outputs against "
                                    "samplesheet")
            # Check outputs and QC summary report
            if not os.path.isdir(self.mkfastq_out_dir):
                raise Exception("No output directory '%s'" %
                                self.mkfastq_out_dir)
            if self.expect_qc_summary_json:
                json_file = os.path.join(self.mkfastq_out_dir,
                                         "outs",
                                         "qc_summary.json")
                if not os.path.exists(json_file):
                    raise Exception("%s mkfastq failed to make JSON "
                                    "QC summary file (%s not found)"
                                    % (self.pkg,json_file))
                # Make HTML QC summary
                html_file = "%s_qc_summary%s.html" % \
                            (self.pkg,
                             "_%s" %
                             ''.join([str(l) for l in self.lanes])
                             if self.lanes is not None else "")
                make_qc_summary_html(json_file,html_file)
                print("Moving QC report '%s'" % html_file)
                os.rename(html_file,
                          os.path.join(
                              os.path.dirname(self.args.out_dir),
                              os.path.basename(html_file)))
            # Move outputs to final location
            print("Moving output to final location")
            os.rename(self.tmp_out_dir,self.args.out_dir)

class MergeFastqs(PipelineTask):
    """
    Merges Fastqs across multiple lanes
    """
    def init(self,fastq_dirs,out_dir,sample_sheet=None,
             no_lane_splitting=False,create_empty_fastqs=False):
        """
        Initialise the MergeFastqs task

        Arguments:
          fastq_dirs (list): set of directories with
            Fastqs in bcl2fastq-like structure, to
            merge together
          out_dir (str): path to output directory
          sample_sheet (str): optional sample sheet file to
            verify the merged files against
          no_lane_splitting (bool): if True then merge
            Fastqs across lanes
          create_empty_fastqs (bool): if True then create empty
            placeholder Fastq files for any that are missing
            on successful completion of Fastq merging

        Outputs:
            missing_fastqs: list of Fastqs missing after
              Fastq merging
        """
        # Internal variables
        self.tmp_merge_dir = None
        # Outputs
        self.add_output('missing_fastqs',list())
    def setup(self):
        # Sort out the Fastqs in the input project directories
        projects = dict()
        for fastq_dir in self.args.fastq_dirs:
            print("Examining %s" % fastq_dir)
            illumina_data = IlluminaData(os.path.dirname(fastq_dir),
                                         os.path.basename(fastq_dir))
            for project in illumina_data.projects:
                # Collect Fastqs in each sample for each project
                if project.name not in projects:
                    projects[project.name] = dict()
                for sample in project.samples:
                    if sample.name not in projects[project.name]:
                        projects[project.name][sample.name] = list()
                    projects[project.name][sample.name].extend(
                        [os.path.join(sample.dirn,fq)
                         for fq in sample.fastq])
            # Add the 'undetermined' Fastqs
            if illumina_data.undetermined:
                if '__undetermined__' not in projects:
                    projects['__undetermined__'] = dict(Undetermined=list())
                for sample in illumina_data.undetermined.samples:
                    projects['__undetermined__']['Undetermined'].extend(
                        [os.path.join(sample.dirn,fq)
                         for fq in sample.fastq])
        # Assign indices to samples based on sample sheet
        s_indices = {}
        if self.args.sample_sheet:
            sample_sheet = SampleSheet(self.args.sample_sheet)
            s_index = 0
            for line in sample_sheet:
                sample_name = line[sample_sheet.sample_name_column]
                if sample_name not in s_indices:
                    s_index += 1
                    s_indices[sample_name] = s_index
        # Make a new directory for the merging
        self.tmp_merge_dir = "__mergefastqs.%s.tmp" % \
                             time.strftime("%Y%m%d%H%M%S")
        os.mkdir(self.tmp_merge_dir)
        print("Made temporary directory for merging: %s" %
              self.tmp_merge_dir)
        # Merge each project
        for project in projects:
            print("Merging Fastqs in project '%s'" % project)
            if project != '__undetermined__':
                project_dir = os.path.join(self.tmp_merge_dir,project)
                os.mkdir(project_dir)
            else:
                project_dir = self.tmp_merge_dir
            # Handle each sample
            for sample in projects[project]:
                print("* %s" % sample)
                # Group the Fastqs within the sample
                fastq_groups = group_fastqs_by_name(
                    projects[project][sample])
                for idx in range(len(fastq_groups[0])):
                    fastqs = [grp[idx] for grp in fastq_groups]
                    if self.args.no_lane_splitting:
                        # Merge Fastqs across lanes
                        fastq_out = AnalysisFastq(fastqs[0])
                        try:
                            # Reset the Fastq sample index number
                            fastq_out.sample_number = \
                                s_indices[fastq_out.sample_name]
                        except KeyError:
                            pass
                        fastq_out.lane_number = None
                        fastq_out = os.path.join(project_dir,
                                                 "%s.fastq.gz" %
                                                 str(fastq_out))
                        if len(fastqs) > 1:
                            # Concatenate multiple Fastqs
                            concat_cmd = Command('zcat')
                            concat_cmd.add_args(*fastqs)
                            concat_cmd.add_args('|','gzip','-c', #'--best',
                                                '>',fastq_out)
                            self.add_cmd(PipelineCommandWrapper(
                                "Concatenate Fastqs for %s/%s" % (project,
                                                                  sample),
                                *concat_cmd.command_line))
                        else:
                            # Hard link to single Fastq
                            os.link(fastqs[0],os.path.join(fastq_out))
                    else:
                        # No merging across lanes
                        for fq in fastqs:
                            # Construct final Fastq name
                            fastq_out = AnalysisFastq(fq)
                            try:
                                # Reset the Fastq sample index number
                                fastq_out.sample_number = \
                                    s_indices[fastq_out.sample_name]
                            except KeyError:
                                pass
                            # Hard link to individual Fastqs
                            os.link(fq,os.path.join(project_dir,
                                                    "%s.fastq.gz" %
                                                    str(fastq_out)))
        # Add 'Stats' and 'Reports' directories
        for d in ("Stats","Reports"):
            os.mkdir(os.path.join(self.tmp_merge_dir,d))
            print("Made '%s' subdirectory" % d)
    def finish(self):
        if self.tmp_merge_dir:
            if self.args.sample_sheet:
                # Verify outputs and get list of missing Fastqs (if any)
                missing_fastqs = verify_run(self.tmp_merge_dir,
                                            self.args.sample_sheet)
                # Update output
                for fq in missing_fastqs:
                    self.output.missing_fastqs.append(fq)
                # Handle missing fastqs
                if missing_fastqs:
                    if self.args.create_empty_fastqs:
                        # Create empty placeholder Fastqs
                        create_placeholder_fastqs(missing_fastqs,
                                                  base_dir=self.tmp_merge_dir)
                    else:
                        # Terminate with an exception
                        raise Exception("Failed to verify outputs against "
                                        "samplesheet")
            print("Moving merged Fastqs to final location")
            os.rename(self.tmp_merge_dir,self.args.out_dir)

class MergeFastqDirs(PipelineFunctionTask):
    """
    Merges directories with subsets of Fastqs
    """
    def init(self,fastq_dirs,merged_fastq_dir):
        """
        Initialise the MergeFastqDirs task

        Arguments:
          fastq_dirs (list): set of directories with
            Fastqs in bcl2fastq-like structure, to
            merge together
          merged_fastq_dir (str): path to the output
            directory where all the Fastqs will be
            put together
        """
        self.tmp_merge_dir = None
        self.fastq_dirs = None
    def setup(self):
        # Check if output already exists
        if os.path.exists(self.args.merged_fastq_dir):
            print("%s already exists, nothing to do" %
                  self.args.merged_fastq_dir)
            return
        # Explicitly resolve the input parameters
        self.fastq_dirs = [resolve_parameter(d)
                           for d in self.args.fastq_dirs]
        # Collect input directories
        projects = []
        undetermined_fastqs = []
        for d in self.fastq_dirs:
            # Check the projects
            print("Examining %s" % d)
            illumina_data = IlluminaData(os.path.dirname(d),
                                         os.path.basename(d))
            for project in illumina_data.projects:
                if not [p for p in projects if p.name == project.name]:
                    print("- %s: will be merged" % project.name)
                    projects.append(project)
                else:
                    print("- %s: another project with the same name "
                          "already found")
                    raise Exception("collision: another project "
                                    "called '%s' was already found" %
                                    project.name)
            # Collect Fastqs with undetermined reads
            if illumina_data.undetermined:
                for sample in illumina_data.undetermined.samples:
                    for fq in sample.fastq:
                        undetermined_fastqs.append(
                            os.path.join(sample.dirn,fq))
        # Make a new directory for the merging
        self.tmp_merge_dir = "__mergefastqdirs.tmp"
        os.mkdir(self.tmp_merge_dir)
        print("Made temporary directory for merging: %s" %
              self.tmp_merge_dir)
        # Handle each project
        for project in projects:
            print("- Importing project '%s'" % project.name)
            self.add_call("Copy project %s" % project.name,
                          self.copy_project,
                          project.dirn,
                          os.path.join(self.tmp_merge_dir,
                                       os.path.basename(project.dirn)))
        # Handle the undetermined Fastqs
        undetermined = {}
        for fq in undetermined_fastqs:
            illuminafq = IlluminaFastq(fq)
            idx = "%s%s%d" % (('L%03d_' % illuminafq.lane_number
                               if illuminafq.lane_number else ''),
                              ('I' if illuminafq.is_index_read
                               else 'R'),
                              illuminafq.read_number)
            try:
                undetermined[idx].append(fq)
            except KeyError:
                undetermined[idx] = [fq]
        for idx in undetermined:
            fastqs_in = sorted(undetermined[idx])
            fastq_out = os.path.join(self.tmp_merge_dir,
                                     "Undetermined_S0_%s_001.fastq.gz"
                                     % idx)
            if len(fastqs_in) == 1:
                # Only one Fastq in list, copy it
                fq = fastqs_in[0]
                print("Copying %s" % fq)
                self.add_call("Linking to %s" % fq,
                              os.link,
                              fq,
                              fastq_out)
            else:
                # Multiple Fastqs, concat them
                print("Concatenting '%s' Fastqs:" % idx)
                for fq in fastqs_in:
                    print("- %s" % fq)
                concat_cmd = Command('zcat')
                concat_cmd.add_args(*fastqs_in)
                concat_cmd.add_args('|','gzip','-c',
                                    '>',fastq_out)
                self.add_cmd(PipelineCommandWrapper(
                    "Concatenate %s undetermined Fastqs" % idx,
                    *concat_cmd.command_line))
        # Add 'Stats' and 'Reports' directories
        for d in ("Stats","Reports"):
            os.mkdir(os.path.join(self.tmp_merge_dir,d))
            print("Made '%s' subdirectory" % d)
    def copy_project(self,project_dir,dest):
        # Copy contents of project to destination dir
        if not os.path.exists(dest):
            print("Making directory %s" % dest)
            os.mkdir(dest)
        # Walk through the project directory
        for f in walk(project_dir):
            if os.path.isdir(f):
                # Make an equivalent subdirectory
                d = os.path.join(dest,os.path.relpath(f,project_dir))
                if not os.path.exists(d):
                    print("Making subdirectory %s" % d)
                    os.mkdir(d)
            elif os.path.isfile(f):
                # Make a hard link to the file
                src = os.path.join(project_dir,f)
                tgt = os.path.join(dest,os.path.relpath(f,project_dir))
                print("Making link to file %s" % src)
                os.link(src,tgt)
            else:
                raise Exception("Don't know how to handle %s" %
                                os.path.join(project_dir,f))
    def finish(self):
        if self.tmp_merge_dir:
            print("Moving merged Fastq dir files to final location")
            os.rename(self.tmp_merge_dir,self.args.merged_fastq_dir)
            for fastq_dir in self.fastq_dirs:
                fastq_dir_backup = os.path.join(
                    os.path.dirname(fastq_dir),
                    "save.%s" % os.path.basename(fastq_dir))
                print("Moving source directory '%s' out of the way" %
                      fastq_dir)
                os.rename(fastq_dir,fastq_dir_backup)

class FastqStatistics(PipelineTask):
    """
    Generates statistics for Fastq files
    """
    def init(self,bcl2fastq_dir,sample_sheet,out_dir,
             stats_file=None,stats_full_file=None,
             per_lane_stats_file=None,
             per_lane_sample_stats_file=None,
             add_data=False,force=False,nprocessors=None):
        """
        Initialise the FastqStatistics task

        Arguments:
          bcl2fastq_dir (str): path to directory with
            Fastqs from bcl2fastq
          sample_sheet (str): path to sample sheet file
          out_dir (str): path to directory to write the
            output stats files to
          stats_file (str): path to statistics output file
          stats_full_file (str): path to full statistics
            output file
          per_lane_stats_file (str): path to per-lane
            statistics output file
          per_lane_sample_stats_file (str): path to
            per-lane per-sample statistics output file
          add_data (bool): if True then add stats to the
            existing stats files (default is to overwrite
            existing stats files)
          force (bool): if True then force update of the
            stats files even if they are newer than the
            Fastq files (by default stats are only updated
            if they are older than the Fastqs)
          nprocessors (int): number of cores to use when
            running 'fastq_statistics.py'

        Outputs:
          stats_file: path to basic stats file
          stats_full: path to full stats file
          per_lane_stats: path to per-lane stats file
          per_lane_sample_stats: path to per-lane sample
            stats file
        """
        # Flag to indicate if statistics should be (re)generated
        self.generate_stats = False
        # Names for intermediate output stats files
        self.stats_file = "statistics.info"
        self.stats_full = "statistics_full.info"
        self.per_lane_stats = "per_lane_statistics.info"
        self.per_lane_sample_stats = "per_lane_sample_stats.info"
        # Outputs
        self.add_output("stats_file",Param(type=str))
        self.add_output("stats_full",Param(type=str))
        self.add_output("per_lane_stats",Param(type=str))
        self.add_output("per_lane_sample_stats",Param(type=str))
    def setup(self):
        # Sort out final output file names
        # Basic statistics
        if self.args.stats_file:
            self.final_stats = self.args.stats_file
        else:
            self.final_stats = os.path.join(self.args.out_dir,
                                            self.stats_file)
        # Per-lane statistics
        if self.args.per_lane_stats_file:
            self.final_per_lane_stats = \
                                self.args.per_lane_stats_file
        else:
            self.final_per_lane_stats = os.path.join(
                self.args.out_dir,self.per_lane_stats)
        # Full statistics
        if self.args.stats_full_file:
            self.final_stats_full = self.args.stats_full_file
        else:
            self.final_stats_full = os.path.join(self.args.out_dir,
                                                 self.stats_full)
        # Per-lane per-sample statistics
        if self.args.per_lane_sample_stats_file:
            self.final_per_lane_sample_stats = \
                self.args.per_lane_sample_stats_file
        else:
            self.final_per_lane_sample_stats = os.path.join(
                self.args.out_dir,self.per_lane_sample_stats)
        # Get most recent timestamp on existing files
        newest_mtime = 0
        for f in (self.final_stats,
                  self.final_per_lane_stats,
                  self.final_stats_full,
                  self.final_per_lane_sample_stats):
            try:
                # Update newest timestamp
                newest_mtime = max(newest_mtime,
                                   os.path.getmtime(f))
            except OSError:
                # File doesn't exist
                newest_mtime = 0
                self.generate_stats = True
        # Load data from bcl2fastq outputs
        try:
            illumina_data = IlluminaData(
                os.path.dirname(self.args.bcl2fastq_dir),
                os.path.basename(self.args.bcl2fastq_dir))
        except Exception as ex:
            raise Exception("Failed to load bcl2fastq data from %s: "
                            "%s" % (self.args.bcl2fastq_dir,ex))
        # Check if any Fastqs are newer than stats
        if newest_mtime > 0:
            for project in illumina_data.projects:
                for sample in project.samples:
                    for fq in sample.fastq:
                        fastq = os.path.join(sample.dirn,fq)
                        if (os.path.getmtime(fastq) > newest_mtime):
                            self.generate_stats = True
                            break
            if self.generate_stats:
                print("Fastqs are newer than existing stats files")
            else:
                print("Stats files are newer than Fastqs")
                if self.args.force:
                    # Force regenerate the statistics
                    self.generate_stats = True
        # Don't continue if nothing to do
        if not self.generate_stats:
            print("Nothing to do")
            return
        # Number of processors
        if self.args.nprocessors:
            nprocessors = self.args.nprocessors
        else:
            nprocessors = self.runner_nslots
        # Run the fastq_statistics.py utility
        fastq_statistics_cmd = Command(
            'fastq_statistics.py',
            '--unaligned',os.path.basename(self.args.bcl2fastq_dir),
            '--sample-sheet',self.args.sample_sheet,
            '--output',self.stats_file,
            '--full-stats',self.stats_full,
            '--per-lane-stats',self.per_lane_stats,
            '--per-lane-sample-stats',self.per_lane_sample_stats,
            '--nprocessors',nprocessors)
        if self.args.add_data:
            fastq_statistics_cmd.add_args('--update')
        fastq_statistics_cmd.add_args(
            os.path.dirname(self.args.bcl2fastq_dir))
        print("Running %s" % fastq_statistics_cmd)
        self.add_cmd(PipelineCommandWrapper(
            "Run fastq_statistics",
            *fastq_statistics_cmd.command_line))
    def finish(self):
        if self.generate_stats:
            print("Moving stats files to final locations")
            for f,ff in ((self.stats_file,self.final_stats),
                         (self.per_lane_stats,self.final_per_lane_stats),
                         (self.stats_full,self.final_stats_full),
                         (self.per_lane_sample_stats,
                          self.final_per_lane_sample_stats)):
                print("- %s -> %s" % (f,ff))
                if os.path.abspath(f) == os.path.abspath(ff):
                    raise Exception("'%s' and '%s' are the same file?"
                                    % (f,ff))
                elif not os.path.exists(f):
                    raise Exception("'%s' not found in %s" % (f,os.getcwd()))
                elif not os.path.isdir(os.path.dirname(ff)):
                    raise Exception("No path to '%s'" % ff)
                os.rename(f,ff)
        # Assign outputs
        self.output.stats_file.set(self.final_stats)
        self.output.stats_full.set(self.final_stats_full)
        self.output.per_lane_stats.set(self.final_per_lane_stats)
        self.output.per_lane_sample_stats.set(
            self.final_per_lane_sample_stats)

class ReportProcessingQC(PipelineTask):
    """
    Generate HTML report on the processing QC
    """
    def init(self,name,analysis_dir,stats_file,per_lane_stats_file,
             per_lane_sample_stats_file,report_html):
        """
        Initialise the ReportProcessingQC task

        Arguments:
          name (str): identifier for report title
          analysis_dir (str): directory with the
            statistics files
          stats_file (str): path to full statistics
            file
          per_lane_stats_file (str): path to the
            per-lane statistics file
          per_lane_sample_stats_file (str): path to
            the per-lane per-sample statistics file
          report_html (str): path to the output
            HTML QC report
        """
        self.tmp_report = None
    def setup(self):
        print("Generating processing QC report")
        self.tmp_report = os.path.join(self.args.analysis_dir,
                                       "processing_qc_report.html.tmp")
        ProcessingQCReport(
            self.args.analysis_dir,
            self.args.stats_file,
            self.args.per_lane_stats_file,
            self.args.per_lane_sample_stats_file,
            name=self.args.name).\
            write(self.tmp_report)
    def finish(self):
        print("Moving processing QC report to final location")
        print("- %s" % self.args.report_html)
        os.rename(self.tmp_report,self.args.report_html)

def subset(lanes,**kws):
    """
    Create a dictionary representing a set of lanes

    Returns a dictionary which holds information
    about a set of lanes grouped together for
    processing, along with values of parameters
    that should be used for this set of lanes.

    Keys must be one of the parameter names
    listed in the LANE_SET_ATTRIBUTES constant;
    specifying an unrecognised key will result
    in a KeyError exception.

    Arguments:
      lanes (list): lanes that comprise the set
      kws (mapping): set of key-value
        pairs assigning values to parameters
        for the group of lanes

    Raises:
      KeyError: if a supplied key is not a valid
        attribute.
    """
    s = dict(lanes=sorted([int(l) for l in lanes]))
    for attr in kws:
        if attr not in LANE_SUBSET_ATTRS:
            raise KeyError("Unsupported subset attribute '%s'"
                           % attr)
        s[attr] = kws[attr]
    return s

def verify_run(fastq_dir,sample_sheet):
    """
    Verify Fastq dir contents against sample sheet

    Check the contents of a Bcl-to-Fastq output directory
    against a sample sheet, and return a list of missing
    Fastqs (or an empty list if all expected Fastqs are
    present).

    Arguments:
      fastq_dir (str): path to Bcl-to-Fastq output directory
      sample_sheet (str): path to sample sheet file

    Returns:
      List: list of missing Fastqs, or an empty list if
        all expected Fastqs are present.
    """
    illumina_data = IlluminaData(os.path.dirname(fastq_dir),
                                 os.path.basename(fastq_dir))
    if verify_run_against_sample_sheet(illumina_data,sample_sheet):
        # Verification ok
        print("Verified Fastqs against samplesheet")
        return []
    else:
        # Verification failed
        print("Failed to verify outputs against samplesheet")
        # List the missing Fastq files
        missing_fastqs = list_missing_fastqs(illumina_data,
                                             sample_sheet)
        print("Missing Fastqs:")
        for fq in missing_fastqs:
            print("- %s" % fq)
        return missing_fastqs

def create_placeholder_fastqs(fastqs,base_dir=None):
    """
    Create empty 'placeholder' Fastq files

    Arguments:
      fastqs (list): paths to Fastq file names to
        create
      base_dir (str): if supplied then used as the base
        directory; Fastqs will be created relative to
        this dir
    """
    # Create empty placeholder Fastqs
    print("Making empty placeholder Fastqs")
    for fastq in fastqs:
        # Prepend base dir
        if base_dir:
            fastq = os.path.join(base_dir,fastq)
        # Make intermediate directory if required
        if not os.path.exists(os.path.dirname(fastq)):
            os.mkdir(os.path.dirname(fastq))
        # Make empty file
        with gzip.GzipFile(filename=fastq,mode='wb') as fp:
            fp.write(''.encode())
