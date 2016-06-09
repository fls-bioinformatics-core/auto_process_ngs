#!/usr/bin/env python
#
#     bcl2fastq_utils.py: utility functions for bcl2fastq conversion
#     Copyright (C) University of Manchester 2013-2016 Peter Briggs
#
########################################################################
#
# bclToFastq.py
#
#########################################################################

"""
bcl2fastq_utils.py

Utility functions for bcl to fastq conversion operations:

- available_bcl2fastq_versions: list available bcl2fastq converters
- bcl_to_fastq_info: retrieve information on the bcl2fastq software
- make_custom_sample_sheet: create a fixed copy of a sample sheet file
- get_required_samplesheet_format: fetch format required by bcl2fastq version
- get_nmismatches: determine number of mismatches from bases mask
- check_barcode_collisions: look for too-similiar pairs of barcode sequences
- get_bases_mask: get a bases mask string
- run_bcl2fastq_1_8: run bcl-to-fastq conversion from CASAVA/bcl2fastq 1.8.*
- run_bcl2fastq_2_17: run bcl-to-fastq conversion from bcl2fastq 2.17.*

"""

#######################################################################
# Imports
#######################################################################
import os
import re
import operator
import logging
import auto_process_ngs.applications as applications
import bcftbx.IlluminaData as IlluminaData
import bcftbx.utils as bcf_utils
from pkg_resources import parse_version

#######################################################################
# Functions
#######################################################################

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
    # Search paths
    if paths is None:
        paths = os.environ['PATH'].split(os.pathsep)
    # Search for executables
    available_exes = []
    for path in paths:
        if os.path.isfile(path):
            path = os.path.dirname(path)
        for name in ('bcl2fastq','configureBclToFastq.pl',):
            prog_path = os.path.abspath(os.path.join(path,name))
            if bcf_utils.PathInfo(prog_path).is_executable:
                available_exes.append(prog_path)
    # Filter on requirement
    if reqs:
        # Loop over ranges
        for req in reqs.split(','):
            logging.debug("Filtering on expression: %s" % req)
            # Determine operator and version
            req_op = None
            req_version = None
            for op in ('==','>=','<=','>','<'):
                if req.startswith(op):
                    req_op = op
                    req_version = req[len(op):].strip()
                    break
            if req_version is None:
                req_op = '=='
                req_version = req.strip()
            logging.debug("Required version: %s %s" % (req_op,req_version))
            if req_op == '==':
                op = operator.eq
            elif req_op == '>=':
                op = operator.ge
            elif req_op == '>':
                op = operator.gt
            elif req_op == '<':
                op = operator.lt
            elif req_op == '<=':
                op = operator.le
            # Filter the available executables on version
            logging.debug("Pre filter: %s" % available_exes)
            logging.debug("Versions  : %s" % [bcl_to_fastq_info(x)[2]
                                              for x in available_exes])
            available_exes = filter(lambda x: op(
                parse_version(bcl_to_fastq_info(x)[2]),
                parse_version(req_version)),
                                    available_exes)
            logging.debug("Post filter: %s" % available_exes)
        # Sort into version order, highest to lowest
        available_exes.sort(
            key=lambda x: parse_version(bcl_to_fastq_info(x)[2]),
            reverse=True)
        logging.debug("Post sort: %s" % available_exes)
        logging.debug("Versions : %s" % [bcl_to_fastq_info(x)[2]
                                 for x in available_exes])
    return available_exes

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
        version_cmd = applications.Command(bcl2fastq_path,'--version')
        output = version_cmd.subprocess_check_output()[1]
        for line in output.split('\n'):
            if line.startswith('bcl2fastq'):
                # Extract version from line of the form
                # bcl2fastq v2.17.1.14
                package_name = 'bcl2fastq'
                try:
                    package_version = line.split()[1][1:]
                except ex:
                    logging.warning("Unable to get version from '%s': %s" %
                                    (line,ex))
    else:
        # No package supplied or located
        logging.warning("Unable to identify bcl-to-fastq conversion package "
                        "from '%s'" % bcl2fastq_path)
    # Return what we found
    return (bcl2fastq_path,package_name,package_version)

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None,
                             fmt=None):
    """
    Creates a corrected copy of a sample sheet file

    Creates and returns a SampleSheet object with a copy of the
    input sample sheet, with any illegal or duplicated names fixed.
    Optionally also writes the updated sample sheet data to a new
    file.

    Arguments:
      input_sample_sheet (str): name and path of the original sample
        sheet file
      output_sample_sheet (str): (optional) name and path to write
        updated sample sheet to, or `None`
      fmt (str): (optional) format for the output sample sheet,
        either 'CASAVA' or 'IEM'; if this `None` then the format of
        the original file will be used

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

def get_bases_mask(run_info_xml,sample_sheet_file):
    """
    Get bases mask string

    Generates initial bases mask based on data in RunInfo.xml (which
    says how many reads there are, how many cycles in each read, and
    which are index reads). Then updates this using the barcode
    information in the sample sheet file.

    Arguments:
      run_info_xml: name and path of RunInfo.xml file from the
        sequencing run
      sample_sheet_file: name and path of sample sheet file.

    Returns:
      Bases mask string e.g. 'y101,I6'. 

    """
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print "Bases mask: %s (from RunInfo.xml)" % bases_mask
    # Update bases mask from sample sheet
    example_barcode = IlluminaData.get_casava_sample_sheet(sample_sheet_file)[0]['Index']
    bases_mask = IlluminaData.fix_bases_mask(bases_mask,example_barcode)
    print "Bases mask: %s (updated for barcode sequence '%s')" % (bases_mask,
                                                                  example_barcode)
    return bases_mask

def get_nmismatches(bases_mask):
    """
    Determine number of mismatches from bases mask

    Automatically determines the maximum number of mismatches that shoud
    be allowed for a bcl to fastq conversion run, based on the tag
    length i.e. the length of the index barcode sequences.

    Tag lengths of 6 or more use 1 mismatch, otherwise use zero
    mismatches.

    The number of mismatches should be supplied to the bclToFastq
    conversion process.

    Arguments:
      bases_mask: bases mask string of the form e.g. 'y101,I6,y101'

    Returns:
      Integer value of number of mismatches. (If the bases mask doesn't
      contain any index reads then returns zero.)

    """
    for read in bases_mask.split(','):
        if read.startswith('I'):
            try:
                i = read.index('n')
                read = read[:i]
            except ValueError:
                pass
            index_length = int(read[1:].rstrip('n'))
            if index_length >= 6:
                return 1
            else:
                return 0
    # Failed to find any indexed reads
    return 0

def check_barcode_collisions(input_sample_sheet,nmismatches):
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
      input_sample_sheet (str): path to a SampleSheet.csv file
        to analyse for barcode collisions
      nmismatches (int): maximum number of mismatches to allow

    Returns:
      List: list of pairs of colliding barcodes (with each pair
        wrapped in a tuple), or an empty list if no collisions
        were detected.

    """
    # Load the sample sheet data
    sample_sheet = IlluminaData.SampleSheet(input_sample_sheet)
    # List of index sequences (barcodes)
    barcodes = {}
    has_lanes = sample_sheet.has_lanes
    for line in sample_sheet:
        # Lane
        if has_lanes:
            lane = line['Lane']
        else:
            lane = 1
        # Index sequence
        try:
            # Try dual-indexed IEM4 format
            indx = "%s%s" %(line['index'].strip(),
                            line['index2'].strip())
        except KeyError:
            # Try single indexed IEM4 (no index2)
            try:
                indx = line['index'].strip()
            except KeyError:
                # Try CASAVA format
                indx = line['Index'].strip()
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

def run_bcl2fastq_1_8(basecalls_dir,sample_sheet,
                      output_dir="Unaligned",
                      mismatches=None,
                      bases_mask=None,
                      nprocessors=None,
                      force=False,
                      ignore_missing_bcl=False,
                      ignore_missing_stats=False,
                      ignore_missing_control=False):
    """
    Wrapper for running the CASAVA/bcl2fastq 1.8.* pipeline

    Runs the CASAVA-style bcl to fastq pipeline, specifically:

    1. Executes the 'configureBclToFastq.pl' script to generate a
       Makefile to perform the conversion, then
    2. Executes 'make' using this Makefile to generate fastq files
       from bcl files.

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
      nprocessors: optional, number of processors to use when running
        'make' step
      force: optional, if True then force overwrite of an existing
        output directory (default is False).
      ignore_missing_bcl: optional, if True then interpret missing bcl
        files as no call (default is False)
      ignore_missing_stats: optional, if True then fill in with zeroes
        when *.stats files are missing (default is False)
      ignore_missing_control: optional, if True then interpret missing
        control files as not-set control bits (default is False)

    Returns:
      0 on success; if a problem is encountered then returns -1 for
      errors within the function (e.g. missing Makefile) or the exit
      code from the failed program.

    """
    # Set up and run configureBclToFastq
    configure_cmd = applications.bcl2fastq.configureBclToFastq(
        basecalls_dir,
        sample_sheet,
        output_dir=output_dir,
        mismatches=mismatches,
        bases_mask=bases_mask,
        force=force,
        ignore_missing_bcl=ignore_missing_bcl,
        ignore_missing_stats=ignore_missing_stats,
        ignore_missing_control=ignore_missing_control
    )
    # Check the executable exists
    if not configure_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % configure_cmd.command)
        return -1
    print "Running command: %s" % configure_cmd
    returncode = configure_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("configureToBclFastq.pl returned %s" % returncode)
        return returncode
    # Check outputs (directory and makefile)
    if not os.path.isdir(output_dir):
        logging.error("Output directory '%s' not found" % output_dir)
        return -1
    makefile = os.path.join(output_dir,'Makefile')
    if not os.path.isfile(makefile):
        logging.error("Makefile not found in %s" % output_dir)
        return -1
    # Set up and run make command
    make_cmd = applications.general.make(makefile=makefile,
                                         working_dir=output_dir,
                                         nprocessors=nprocessors)
    if not make_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % make_cmd.command)
        return -1
    print "Running command: %s" % make_cmd
    returncode = make_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("make returned %s" % returncode)
    return returncode

def run_bcl2fastq_2_17(basecalls_dir,sample_sheet,
                       output_dir="Unaligned",
                       mismatches=None,
                       bases_mask=None,
                       ignore_missing_bcl=False,
                       no_lane_splitting=False,
                       minimum_trimmed_read_length=None,
                       mask_short_adapter_reads=None,
                       loading_threads=None,
                       demultiplexing_threads=None,
                       processing_threads=None,
                       writing_threads=None):
    """
    Wrapper for running bcl2fastq 2.17.*

    Runs the bcl2fastq 2.17.* software to generate fastq files
       from bcl files.

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
      ignore_missing_bcl: optional, if True then interpret missing bcl
        files as no call (default is False)
      no_lane_splitting: optional, if True then output FASTQ files
        will not be split by lane (default is False, i.e. do split
        FASTQs by lane)
      minimum_trimmed_read_length: optional, specify minimum length
        for reads after adapter trimming (shorter reads will be padded
        with Ns to make them long enough)
      mask_short_adapter_reads: optional, specify the minimum length
        of ACGT bases that must be present in a read after adapter
        trimming for it not to be masked completely with Ns.
      loading_threads: optional, specify number of threads to use
        for loading bcl data (--loading-threads)
      demultiplexing_threads: optional, specify number of threads to
        use for demultiplexing (--demultiplexing-threads)
      processing_threads: optional, specify number of threads to use
        for processing (--processing-threads)
      writing_threads: optional, specify number of threads to use for
        writing FASTQ data (--writing-threads)

    Returns:
      0 on success; if a problem is encountered then returns -1 for
      errors within the function (e.g. missing Makefile) or the exit
      code from the failed program.

    """
    # Set up and run bcl2fastq2
    bcl2fastq2_cmd = applications.bcl2fastq.bcl2fastq2(
        basecalls_dir,
        sample_sheet,
        output_dir=output_dir,
        mismatches=mismatches,
        bases_mask=bases_mask,
        ignore_missing_bcl=ignore_missing_bcl,
        no_lane_splitting=no_lane_splitting,
        minimum_trimmed_read_length=minimum_trimmed_read_length,
        mask_short_adapter_reads=mask_short_adapter_reads,
        loading_threads=loading_threads,
        demultiplexing_threads=demultiplexing_threads,
        processing_threads=processing_threads,
        writing_threads=writing_threads
    )
    # Check the executable exists
    if not bcl2fastq2_cmd.has_exe:
        logging.error("'%s' missing, cannot run" % bcl2fastq2_cmd.command)
        return -1
    print "Running command: %s" % bcl2fastq2_cmd
    returncode = bcl2fastq2_cmd.run_subprocess()
    # Check returncode
    if returncode != 0:
        logging.error("bcl2fastq returned %s" % returncode)
        return returncode
    # Check outputs (directory and makefile)
    if not os.path.isdir(output_dir):
        logging.error("Output directory '%s' not found" % output_dir)
        return -1
    return returncode
