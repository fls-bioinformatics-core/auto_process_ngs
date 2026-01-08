#!/usr/bin/env python
#
#     tenx/multiome.py: utilities for handling 10xGenomics SC multiome
#     Copyright (C) University of Manchester 2023 Peter Briggs
#

"""
Utilities for working with 10x Genomics single cell multiome
pipelines:

- MultiomeLibraries
"""

#######################################################################
# Imports
#######################################################################

import os
from ..analysis import locate_project
from ..analysis import split_sample_reference

#######################################################################
# Classes
#######################################################################

class MultiomeLibraries:
    """
    Class to handle '10x_multiome_libraries.info' files

    These files link sample names in an analysis project
    with those in another project. They consist of
    tab-delimited lines of the form:

    <LOCAL_SAMPLE>   <REMOTE_SAMPLE>

    where LOCAL_SAMPLE should be the name of a sample in
    the local project directory, and REMOTE_SAMPLE is a
    compound ID for a sample in a different project, of
    the form:

    [RUN:]PROJECT/SAMPLE_NAME

    In turn, RUN can be a run reference ID, a run name, or
    path to an analysis directory.
    """
    def __init__(self,filen):
        """
        Create new MultiomeLibraries instance

        Arguments:
          filen (str): path to a 10x Multiome libraries
            info file
        """
        self._filen = os.path.abspath(filen)
        self._samples = dict()
        # Read data in from libraries file
        with open(self._filen,'rt') as fp:
            for line in fp:
                # Ignore comment lines
                if line.startswith('#'):
                    continue
                # Lines should be 'local_sample<TAB>remote_sample'
                local_sample,remote_sample = line.split()
                try:
                    self._samples[local_sample].append(remote_sample)
                except KeyError:
                    self._samples[local_sample] = [remote_sample,]

    def _library_description(self,library_type):
        """
        Internal: get cellranger-arc description of library type
        """
        if library_type in ('ATAC', 'snATAC'):
            return "Chromatin Accessibility"
        elif library_type in ('GEX', 'snGEX'):
            return "Gene Expression"
        else:
            raise Exception("Unsupported library: '%s'"
                            % library_type)

    @property
    def local_samples(self):
        """
        List the local sample names
        """
        return sorted(list(self._samples.keys()))

    def linked_samples(self,sample):
        """
        List the remote samples associated with a local sample
        """
        return self._samples[sample]

    def linked_projects(self):
        """
        List the projects linked to this one
        """
        projects = list()
        for sample in self.local_samples:
            for linked_sample in self.linked_samples(sample):
                # Extract and locate the location of the remote sample
                run,project,remote_sample = split_sample_reference(
                    linked_sample)
                linked_project = locate_project(
                    "%s%s" % ('%s:' % run if run else '',
                              project),
                    start_dir=os.path.dirname(self._filen),
                    ascend=True)
                if not linked_project:
                    raise Exception("Failed to locate project for "
                                    "linked sample: %s"
                                    % linked_sample)
                else:
                    # Check linked project not already in the list
                    if linked_project.dirn not in \
                       [p.dirn for p in projects]:
                        projects.append(linked_project)
        # Return list of projects
        return sorted(projects,key=lambda p: p.name)

    def write_libraries_csv(self,sample,fastq_dir,library_type,
                            filen=None):
        """
        Create a cellranger-arc libraries.csv file

        The format is a header line of the form:

        fastqs,sample,library_type

        See https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/using/count#libraries

        Arguments:
          sample (str): local sample name
          fastq_dir (str): path to directory with the
            Fastqs for the sample
          library_type (str): name of the library type
            for the sample
          filen (str): optional, path to the output CSV
            file (defaults to 'libraries.SAMPLE.csv')
        """
        # Determine output file
        if filen is None:
            filen = "libraries.%s.csv" % sample
        filen = os.path.abspath(filen)
        # Get linked samples
        linked_samples = self._samples[sample]
        # Write libraries.csv file for cellranger-arc
        with open(filen,'wt') as fp:
            fp.write("fastqs,sample,library_type\n")
            # Data for local sample
            fp.write("%s,%s,%s\n" %
                     (fastq_dir,
                      sample,
                      self._library_description(library_type)))
            # Data for linked sample(s)
            for sample_id in linked_samples:
                print("Locating linked sample: '%s'" % sample_id)
                run,project,linked_sample = \
                    split_sample_reference(sample_id)
                project = locate_project(
                    "%s%s" % ('%s:' % run if run else '',
                              project),
                    start_dir=os.path.dirname(self._filen),
                    ascend=True)
                if not project:
                    raise Exception("Failed to locate project for "
                                    "linked sample: %s"
                                    % sample_id)
                fp.write("%s,%s,%s\n" %
                         (project.fastq_dir,
                          linked_sample,
                          self._library_description(
                              project.info.library_type)))
            print("Generated %s" % filen)
