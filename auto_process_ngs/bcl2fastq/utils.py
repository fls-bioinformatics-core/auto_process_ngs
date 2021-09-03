#!/usr/bin/env python
#
#     bcl2fastq/utils.py: utility functions for bcl2fastq conversion
#     Copyright (C) University of Manchester 2013-2021 Peter Briggs
#
########################################################################
#
# bcl2fastq/utils.py
#
#########################################################################

"""
bcl2fastq/utils.py

Utility functions for bcl to fastq conversion operations:

- get_sequencer_platform: get sequencing instrument platform
- available_bcl2fastq_versions: list available bcl2fastq converters
- bcl_to_fastq_info: retrieve information on the bcl2fastq software
- bclconvert_info: retrieve information on the BCL Convert software
- make_custom_sample_sheet: create a corrected copy of a sample sheet file
- get_required_samplesheet_format: fetch format required by bcl2fastq version
- get_bases_mask: get a bases mask string
- bases_mask_is_valid: check if bases mask string is valid
- get_nmismatches: determine number of mismatches from bases mask
- convert_bases_mask_to_override_cycles: convert bases mask for BCL convert
- check_barcode_collisions: look for too-similiar pairs of barcode sequences

"""

#######################################################################
# Imports
#######################################################################
import os
import re
import logging
from ..command import Command
from ..utils import find_executables
from ..utils import parse_version
from ..samplesheet_utils import barcode_is_10xgenomics
import bcftbx.IlluminaData as IlluminaData
import bcftbx.platforms as platforms
import bcftbx.utils as bcf_utils
from bcftbx.JobRunner import SimpleJobRunner

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def get_sequencer_platform(dirn,instrument=None,settings=None):
    """
    Return the platform for the sequencing instrument

    Attempts to identify the platform (e.g. 'hiseq', 'miseq' etc)
    for a sequencing run.

    If 'settings' is supplied then the platform is looked up
    based on the instrument names and platforms listed in the
    'sequencers' section of the configuration. If 'instrument'
    is also supplied then this is used; otherwise the instrument
    name is extracted from the supplied directory name.

    If no match can be found then there is a final attempt to
    determine the platform from the hard-coded names in the
    'bcftbx.platforms' module.

    Arguments:
      dirn (str): path to the data or analysis directory
      instrument (str): (optional) the instrument name
      settings (Settings):  (optional) a Settings instance
        with the configuration loaded

    Returns:
      String: either the platform or None, if the platform
        cannot be determined.
    """
    # Attempt to look up the instrument name
    platform = None
    if instrument is None:
        print("Extracting instrument name from directory name")
        try:
            datestamp,instrument,run_number,\
                flow_cell_prefix,flow_cell_id = \
                    IlluminaData.split_run_name_full(dirn)
        except Exception as ex:
            logger.warning("Unable to extract instrument name: "
                           "%s" % ex)
    if instrument and settings:
        print("Identifying platform from instrument name")
        try:
            return settings.sequencers[instrument].platform
        except KeyError:
            # Instrument not listed in the settings
            logger.warning("Instrument name '%s' not found in "
                           "configuration file" % instrument)
    # Fall back to old method
    print("Identifying platform from data directory name")
    platform = platforms.get_sequencer_platform(dirn)
    if platform is None:
        logger.warning("Unable to identify platform from "
                       "directory name")
    return platform

def available_bcl2fastq_versions(reqs=None,paths=None):
    """
    List available bcl2fastq converters

    By default searches the PATH for likely bcl2fastq
    converters and returns a list of executables with the
    full path.

    The 'reqs' argument allows a specific version or range
    of versions to be requested; in this case the returned
    list will only contain those packages which satisfy
    the requested versions.

    A range of version specifications can be requested by
    separating multiple specifiers with a comma - for
    example '>1.8.3,<2.16'.

    The full set of operators is:

    - ==, >, >=, <=, <

    If no versions are requested then the packages will
    be returned in PATH order; otherwise they will be
    returned in version order (highest to lowest).

    Arguments:
      reqs (str): optional version requirement expression
        (for example '>=1.8.4'). If supplied then only
        executables fulfilling the requirement will be
        returned. If no operator is supplied then '=='
        is implied.
      paths (list): optional set of directory paths to
        search when looking for bcl2fastq software. If
        not supplied then the set of paths specified in
        the PATH environment variable will be searched.

    Returns:
      List: full paths to bcl2fastq converter executables.

    """
    return find_executables(('bcl2fastq',
                             'configureBclToFastq.pl'),
                            info_func=bcl_to_fastq_info,
                            reqs=reqs,
                            paths=paths)

def bcl_to_fastq_info(path=None):
    """
    Retrieve information on the bcl2fastq software

    If called without any arguments this will locate the first
    bcl-to-fastq conversion package executable (either
    'configureBclToFastq.pl' or 'bcl2fastq') that is available on
    the user's PATH (as returned by 'available_bcl2fastq_versions')
    and attempts to guess the package name (either `bcl2fastq` or
    `CASAVA`) and the version that it belongs to.

    Alternatively if the path to an executable is supplied then
    the package name and version will be determined from that
    instead.

    If no package is identified then the script path is still
    returned, but without any version info.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the bcl2fastq program or
        configureBclToFastq.pl script and PACKAGE and VERSION are
        guesses for the package/version that it belongs to. If any
        value can't be determined then it will be returned as an
        empty string.

    """
    # Initialise
    bcl2fastq_path = ''
    package_name = ''
    package_version = ''
    # Locate the core script
    if not path:
        exes = available_bcl2fastq_versions()
        if exes:
            bcl2fastq_path = exes[0]
    else:
        bcl2fastq_path = os.path.abspath(path)
    # Identify the version
    if os.path.basename(bcl2fastq_path) == 'configureBclToFastq.pl':
        # Found CASAVA or bcl2fastq 1.8.* version
        # Look for the top-level directory
        path = os.path.dirname(bcl2fastq_path)
        # Look for etc directory
        etc_dir = os.path.join(os.path.dirname(path),'etc')
        if os.path.isdir(etc_dir):
            for d in bcf_utils.list_dirs(etc_dir):
                m = re.match(r'^(bcl2fastq|CASAVA)-([0-9.]+)$',d)
                if m:
                    package_name = m.group(1)
                    package_version = m.group(2)
                    break
    elif os.path.basename(bcl2fastq_path) == 'bcl2fastq':
        # Found bcl2fastq v2.*
        # Run the program to get the version
        version_cmd = Command(bcl2fastq_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            if line.startswith('bcl2fastq'):
                # Extract version from line of the form
                # bcl2fastq v2.17.1.14
                package_name = 'bcl2fastq'
                try:
                    package_version = line.split()[1][1:]
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': %s" %
                                   (line,ex))
    else:
        # No package supplied or located
        logger.warning("Unable to identify bcl-to-fastq conversion package "
                       "from '%s'" % bcl2fastq_path)
    # Return what we found
    return (bcl2fastq_path,package_name,package_version)

def bclconvert_info(path=None):
    """
    Retrieve information on the bcl-convert software

    If called without any arguments this will locate the first
    bcl-concert executable that is available on the user's PATH.

    Alternatively if the path to an executable is supplied then
    the package name and version will be determined from that
    instead.

    If no package is identified then the script path is still
    returned, but without any version info.

    Returns:
      Tuple: tuple consisting of (PATH,PACKAGE,VERSION) where PATH
        is the full path for the bcl-convert program, and PACKAGE
        and VERSION the package/version that it belongs to (PACKAGE
        will be 'BCL Convert' if a matching executable is located).
        If any value can't be determined then it will be returned
        as an empty string.

    """
    # Initialise
    bclconvert_path = ''
    package_name = ''
    package_version = ''
    # Locate the bcl-convert program
    if not path:
        bclconvert_path = bcf_utils.find_program('bcl-convert')
    else:
        bclconvert_path = os.path.abspath(path)
    # Identify the version
    if bclconvert_path:
        # Run the program to get the version
        version_cmd = Command('bcl-convert','-V')
        output = version_cmd.subprocess_check_output()[1]
        print(output)
        for line in output.split('\n'):
            if line.startswith('bcl-convert'):
                # Extract version from line of the form
                # bcl-convert Version 00.000.000.3.7.5
                package_name = 'BCL Convert'
                try:
                    package_version = '.'.join(line.split('.')[-3:])
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': %s" %
                                   (line,ex))
    else:
        # No package supplied or located
        logger.warning("Unable to identify BCLConvert package from '%s'" %
                       bclconvert_path)
    # Return what we found
    return (bclconvert_path,package_name,package_version)

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None,
                             lanes=None,adapter=None,adapter_read2=None,
                             fmt=None):
    """
    Creates a corrected copy of a sample sheet file

    Creates and returns a SampleSheet object with a copy of the
    input sample sheet, with any illegal or duplicated names fixed.
    Optionally it can also: write the updated sample sheet data to a
    new file, switch the format, and include only a subset of lanes
    from the original file

    Arguments:
      input_sample_sheet (str): name and path of the original sample
        sheet file
      output_sample_sheet (str): (optional) name and path to write
        updated sample sheet to, or `None`
      lanes (list): (optional) list of lane numbers to keep in the
        output sample sheet; if `None` then all lanes will be kept
        (the default), otherwise lanes will be dropped if they don't
        appear in the supplied list
      adapter (str): (optional) if set then write to the `Adapter`
        setting
      adapter_read2 (str): (optional) if set then write to the
        `AdapterRead2` setting
      fmt (str): (optional) format for the output sample sheet,
        either 'CASAVA' or 'IEM'; if this is `None` then the format
        of the original file will be used

    Returns:
      SampleSheet object with the data for the corrected sample
      sheet.

    """
    # Load the sample sheet data
    sample_sheet = IlluminaData.SampleSheet(input_sample_sheet)
    # Determine the column names for this format
    if sample_sheet.format == 'CASAVA':
        sample_col = 'SampleID'
        project_col = 'SampleProject'
    elif sample_sheet.format == 'IEM':
        sample_col = 'Sample_ID'
        project_col = 'Sample_Project'
    else:
        raise Exception("Unknown sample sheet format: %s" %
                        sample_sheet.format)
    # Add project names if not supplied
    for line in sample_sheet:
        if not line[project_col]:
            line[project_col] = line[sample_col]
    # Fix other problems
    sample_sheet.fix_illegal_names()
    sample_sheet.fix_duplicated_names()
    # Put data into lane order (if lanes specified)
    if sample_sheet.has_lanes:
        sample_sheet.data.sort(lambda line: line['Lane'])
    # Select subset of lanes if requested
    if lanes is not None:
        logger.debug("Updating to include only specified lanes: %s" %
                     ','.join([str(l) for l in lanes]))
        i = 0
        while i < len(sample_sheet):
            line = sample_sheet[i]
            if line['Lane'] in lanes:
                logger.debug("Keeping %s" % line)
                i += 1
            else:
                del(sample_sheet[i])
    # Update adapter sequences
    if adapter is not None:
        if adapter:
            sample_sheet.settings['Adapter'] = str(adapter)
        else:
            try:
                del(sample_sheet.settings['Adapter'])
            except KeyError:
                pass
    if adapter_read2 is not None:
        if adapter_read2:
            sample_sheet.settings['AdapterRead2'] = str(adapter_read2)
        else:
            try:
                del(sample_sheet.settings['AdapterRead2'])
            except KeyError:
                pass
    # Write out new sample sheet
    if output_sample_sheet is not None:
        sample_sheet.write(output_sample_sheet,fmt=fmt)
    return sample_sheet

def get_required_samplesheet_format(bcl2fastq_version):
    """
    Returns sample sheet format required by bcl2fastq

    Given a bcl2fastq version, returns the format of the
    sample sheet that is required for that version.

    Arguments:
      bcl2fastq_version (str): version of bcl2fastq

    Returns:
      String: Sample sheet format (e.g. 'CASAVA', 'IEM'
        etc).

    """
    version = parse_version(bcl2fastq_version)
    major,minor = version[0:2]
    if (major,minor) == parse_version('1.8')[0:2]:
        # Version 1.8.*
        return 'CASAVA'
    elif major == parse_version('2')[0]:
        # Version 2.*
        return 'IEM'
    else:
        # Not a known version
        raise NotImplementedError('unknown version: %s' %
                                  bcl2fastq_version)

def get_bases_mask(run_info_xml,sample_sheet_file=None):
    """
    Get bases mask string

    Generates initial bases mask based on data in RunInfo.xml (which
    says how many reads there are, how many cycles in each read, and
    which are index reads), and optionally updates this using the
    barcode information in the sample sheet file.

    Arguments:
      run_info_xml: name and path of RunInfo.xml file from the
        sequencing run
      sample_sheet_file: (optional) path to sample sheet file

    Returns:
      Bases mask string e.g. 'y101,I6'. 
    """
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print("Bases mask: %s (from RunInfo.xml)" % bases_mask)
    if sample_sheet_file is not None:
        # Update bases mask from sample sheet
        example_barcode = IlluminaData.samplesheet_index_sequence(
            IlluminaData.SampleSheet(sample_sheet_file).data[0])
        if example_barcode is None:
            example_barcode = ""
        if barcode_is_10xgenomics(example_barcode):
            print("Bases mask: barcode is 10xGenomics sample set ID")
        else:
            bases_mask = IlluminaData.fix_bases_mask(bases_mask,
                                                     example_barcode)
        print("Bases mask: %s (updated for barcode sequence '%s')" %
              (bases_mask,example_barcode))
    return bases_mask

def bases_mask_is_valid(bases_mask):
    """
    Check if a bases mask is valid

    Arguments:
      bases_mask: bases mask string to check

    Returns:
      Boolean: True if the supplied bases mask is valid,
        False if not.
    """
    try:
        for read in bases_mask.upper().split(','):
            if not re.match(r'^([IY][0-9]+|[IY]*)(N[0-9]+|N*)$',read):
                return False
        return True
    except AttributeError:
        return False

def get_nmismatches(bases_mask,multi_index=False):
    """
    Determine number of mismatches from bases mask

    Automatically determines the maximum number of mismatches that
    should be allowed for a bcl to fastq conversion run, based on
    the tag length i.e. the length of the index barcode sequences.

    Tag lengths of 6 or more use 1 mismatch, otherwise use zero
    mismatches.

    The number of mismatches should be supplied to the bclToFastq
    conversion process.

    Raises an exception if the supplied bases mask is not valid.

    Arguments:
      bases_mask: bases mask string of the form e.g. 'y101,I6,y101'
      multi_index: boolean flag, if False (default) then use the
        total length of all indices and return a single integer
        number of allowed mismatches; if True then return a list
        with the number of mismatches for each index (so a dual
        index will be a pair of allowed mismatches)

    Returns:
      Integer value of number of mismatches. (If the bases mask doesn't
        contain any index reads then returns zero for single-index mode,
        or an empty list for multi-read mode.)

    """
    # Check mask is valid
    if not bases_mask_is_valid(bases_mask):
        raise Exception("'%s': not a valid bases mask" % bases_mask)
    # Get the lengths of each index read
    index_lengths = []
    for read in bases_mask.upper().split(','):
        if read.startswith('I'):
            index_length = 0
            try:
                i = read.index('N')
                read = read[:i]
            except ValueError:
                pass
            try:
                index_length += int(read[1:])
            except ValueError:
                index_length += len(read)
            index_lengths.append(index_length)
    if multi_index:
        # Return list of mismatches
        return [1 if index_length >= 6 else 0
                for index_length in index_lengths]
    else:
        # Total the length of all index reads
        index_length = sum(index_lengths)
        # Return number of mismatches
        if index_length >= 6:
            return 1
        else:
            return 0

def convert_bases_mask_to_override_cycles(bases_mask):
    """
    Converts bcl2fastq-format bases mask to BCL Convert format

    Given a bases mask string (e.g. 'y76,I8,I8,y76'), returns
    the equivalent BCL Convert format for use with 'OverrideCycles'
    in a sample sheet (e.g. 'Y76;I8;I8;Y76').

    Arguments:
      bases_mask (str): bcl2fastq bases mask string

    Returns:
      String: the original bases mask converted to BCL Convert
        format
    """
    override_cycles = []
    for item in str(bases_mask).upper().split(','):
        value = ''
        count = 0
        for c in item:
            if not value:
                # First character
                value += c
            elif c.isdigit():
                # Digit character
                value += c
            else:
                # Non-digit character
                if c == value[-1]:
                    # Same as previous character
                    count += 1
                else:
                    # Different from previous character
                    if count:
                        # Increase count by 1 to include
                        # the very first character
                        value += str(count+1)
                        count = 0
                    value += c
        # Tidy up trailing count
        if count:
            # Increase count by 1 to include
            # the very first character
            value += str(count+1)
        # Add converted item
        override_cycles.append(value)
    # Assemble with appropriate delimiter and return
    return ';'.join(override_cycles)

def check_barcode_collisions(sample_sheet_file,nmismatches,
                             use_index='all'):
    """
    Check sample sheet for barcode collisions

    Check barcode index sequences within each lane (or across
    all samples, if no lane information is present) and find
    any which differ in fewer bases than a threshold number
    which is calculated as:

    less than 2 times the number of mismatches plus 1

    (as is stated in the output from bcl2fastq v2.)

    Pairs of barcodes which are too similar (i.e. which collide)
    are reported as a list of tuples, e.g.

    [('ATTCCT','ATTCCG'),...]

    Arguments:
      sample_sheet_file (str): path to a SampleSheet.csv file
        to analyse for barcode collisions
      nmismatches (int): maximum number of mismatches to allow
      use_index (str): flag indicating how to treat index
        sequences: 'all' (the default) combines indexes into a
        single sequence before checking for collisions, '1' only
        checks index 1 (i7), and '2' only checks index 2 (i5)

    Returns:
      List: list of pairs of colliding barcodes (with each pair
        wrapped in a tuple), or an empty list if no collisions
        were detected.

    """
    # Load the sample sheet data
    sample_sheet = IlluminaData.SampleSheet(sample_sheet_file)
    # Convert index flag to string
    use_index = str(use_index)
    # List of index sequences (barcodes)
    barcodes = {}
    has_lanes = sample_sheet.has_lanes
    for line in sample_sheet:
        # Lane
        if has_lanes:
            lane = line['Lane']
        else:
            lane = 1
        # Extract i7 index sequence
        indx_i7 = None
        try:
            # IEM4 format
            indx_i7 = line['index'].strip()
        except KeyError:
            # CASAVA format
            try:
                indx_i7 = line['Index'].strip()
            except KeyError:
                pass
        # Extract i5 index sequence
        indx_i5 = None
        try:
            # IEM4 format
            indx_i5 = line['index2'].strip()
        except KeyError:
            # No i5 for CASAVA
            pass
        # Assemble index sequence to check for mismatches
        if use_index == "all":
            # Combine i5 and i7 into a single sequence
            indx = "%s%s" % (indx_i7 if indx_i7 else '',
                             indx_i5 if indx_i5 else '')
        elif use_index == "1":
            # Only use i7
            indx = indx_i7
        elif use_index == "2":
            # Only use i5
            indx = indx_i5
        else:
            # Undefined index type
            raise Exception("Unrecognised index: '%s'" % use_index)
        # Explicitly set empty index to None
        if not indx:
            indx = None
        try:
            barcodes[lane].append(indx)
        except KeyError:
            barcodes[lane] = [indx,]
    # Mismatch threshold
    mismatch_threshold = 2*nmismatches + 1
    # Check for collisions
    collisions = []
    for lane in barcodes:
        for i,seq1 in enumerate(barcodes[lane][:-1]):
            for seq2 in barcodes[lane][i+1:]:
                ndiff = 0
                for c1,c2 in zip(seq1,seq2):
                    if c1 != c2:
                        ndiff += 1
                if ndiff < mismatch_threshold:
                    collisions.append((seq1,seq2))
    return collisions
