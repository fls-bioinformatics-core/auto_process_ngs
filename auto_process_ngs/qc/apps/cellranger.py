#!/usr/bin/env python
#
#     qc/apps/cellranger: handle outputs from Cellranger variants
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs

"""
Provides utility classes and functions for handline Cellranger outputs.

Provides the following classes:

- CellrangerCount: handle outputs from cellranger count
- CellrangerMulti: handle outputs from cellranger multi

Provides the following functions:

- cellranger_count_output: get names for cellranger count output
- cellranger_atac_count_output: get names for cellranger-atac count output
- cellranger_arc_count_output: get names for cellranger-arc count output
- cellranger_multi_output: get names for cellranger multi output
- fetch_cellranger_multi_output_dirs: get list of cellranger multi output dirs
"""

#######################################################################
# Imports
#######################################################################

import os
from bcftbx.utils import walk
from ...tenx.metrics import GexSummary
from ...tenx.metrics import AtacSummary
from ...tenx.metrics import MultiplexSummary
from ...tenx.metrics import MultiomeSummary
from ...tenx.cellplex import CellrangerMultiConfigCsv

#######################################################################
# Classes
#######################################################################

class CellrangerCount:
    """
    Wrapper class for handling outputs from cellranger count

    The ``CellrangerCount`` object gives access to various
    details of the outputs (such as sample name and file
    paths).
    """
    def __init__(self,cellranger_count_dir,cellranger_exe=None,
                 version=None,reference_data=None):
        """
        Create a new CellrangerCount instance

        Arguments:
          cellranger_count_dir (str): path to the cellranger
            count output directory for a sample (e.g.
            ``.../cellranger_count/SAMPLE``)
          cellranger_exe (str): the cellranger executable that
            generated the outputs
          version (str): the version of cellranger used
          reference_data (str): the reference dataset
        """
        # Store path to top-level directory and sample name
        self._cellranger_count_dir = os.path.abspath(
            cellranger_count_dir)
        try:
            self._sample_name = os.path.basename(self.dir)
        except OSError:
            self._sample_name = None
        # Command line file
        try:
            self._cmdline_file = os.path.join(self.dir,
                                              "_cmdline")
        except OSError:
            self._cmdline_file = None
        # Deal with additional data items
        if not cellranger_exe:
            try:
                cellranger_exe = self.cmdline.split()[0]
            except Exception:
                pass
        self._cellranger_exe = cellranger_exe
        if not reference_data:
            try:
                args = self.cmdline.split()
                for idx,arg in enumerate(args):
                    if arg in ('--transcriptome',
                               '--reference',):
                        reference_data = args[idx+1]
                        break
            except Exception:
                pass
        self._reference_data = reference_data
        self._version = version
        # Paths to metrics and web summary files
        try:
            if self.pipeline_name in ('cellranger',):
                metrics_csv = os.path.join(self.dir,
                                           "outs",
                                           "metrics_summary.csv")
            elif self.pipeline_name in ('cellranger-atac',
                                        'cellranger-arc'):
                metrics_csv = os.path.join(self.dir,
                                           "outs",
                                           "summary.csv")
            else:
                metrics_csv = None
            if metrics_csv:
                if not os.path.exists(metrics_csv):
                    metrics_csv = None
        except OSError:
            # OSError raise by self.dir if top-level is
            # missing
            metrics_csv = None
        self._metrics_csv = metrics_csv
        # Path to web summary HTML file
        try:
            self._web_summary = os.path.join(self.dir,
                                             "outs",
                                             "web_summary.html")
            if not os.path.exists(self._web_summary):
                self._web_summary = None
        except OSError:
            # OSError raise by self.dir if top-level is
            # missing
            self._web_summary = None

    @property
    def mode(self):
        return "count"

    @property
    def dir(self):
        """
        Path to the directory with the cellranger count outputs
        """
        if os.path.isdir(self._cellranger_count_dir):
            return self._cellranger_count_dir
        raise OSError("'%s': not a directory" %
                      self._cellranger_count_dir)

    @property
    def sample_name(self):
        """
        Sample name derived from the directory name
        """
        return self._sample_name

    @property
    def cellranger_exe(self):
        """
        Cellranger executable
        """
        return self._cellranger_exe

    @property
    def pipeline_name(self):
        """
        Pipeline name i.e. name of the software package
        """
        exe = self.cellranger_exe
        if exe:
            if exe.endswith("bin/count"):
                # Special case: 'count' was run directly
                for name in ("cellranger",
                             "cellranger-atac",
                             "cellranger-arc",):
                    if ("%s-cs" % name) in exe.split(os.sep):
                        return name
            # Default: take name of executable
            return os.path.basename(exe)
        else:
            # Couldn't get the executable
            return None

    @property
    def version(self):
        """
        Cellranger version
        """
        return self._version

    @property
    def reference_data(self):
        """
        Reference dataset
        """
        return self._reference_data

    @property
    def metrics_csv(self):
        """
        Path to the cellranger count 'metrics.csv' file
        """
        if self._metrics_csv and os.path.isfile(self._metrics_csv):
            return self._metrics_csv
        raise OSError("'%s': not a file" % self._metrics_csv)

    @property
    def web_summary(self):
        """
        Path to the cellranger count 'web_summary.html' file
        """
        if self._web_summary and os.path.isfile(self._web_summary):
            return self._web_summary
        raise OSError("'%s': not a file" % self._web_summary)

    @property
    def metrics(self):
        """
        Return the appropriate 'MetricsSummary' object
        """
        if self.pipeline_name == "cellranger":
            metrics = GexSummary
        elif self.pipeline_name == "cellranger-atac":
            metrics = AtacSummary
        elif self.pipeline_name == "cellranger-arc":
            metrics = MultiomeSummary
        else:
            metrics = None
        return metrics(self.metrics_csv)

    @property
    def cmdline_file(self):
        """
        Path to the 'cellranger* count' '_cmdline' file
        """
        return self._cmdline_file

    @property
    def cmdline(self):
        """
        Return the command line used to run 'cellranger* count'
        """
        if self.cmdline_file:
            with open(self.cmdline_file,'rt') as fp:
                return fp.read().strip()
        else:
            return None

class CellrangerMulti:
    """
    Wrapper class for handling outputs from cellranger multi

    The ``CellrangerMulti`` object gives access to various
    details of the outputs (such as sample name and file
    paths).
    """
    def __init__(self,cellranger_multi_dir,cellranger_exe=None,
                 version=None,reference_data=None,
                 config_csv=None):
        """
        Create a new CellrangerMulti instance

        Arguments:
          cellranger_multi_dir (str): path to the cellranger
            multi output directory for a sample
          cellranger_exe (str): the cellranger executable that
            generated the outputs
          version (str): the version of cellranger used
          reference_data (str): the reference dataset
          config_csv (str): the config.csv file used for
            running the 'multi' command
        """
        # Store path to top-level directory
        self._cellranger_multi_dir = os.path.abspath(
            cellranger_multi_dir)
        # Identify the samples
        self._samples = {}
        try:
            per_sample_outs_dir = os.path.join(self.dir,
                                               "outs",
                                               "per_sample_outs")
            for d in os.listdir(per_sample_outs_dir):
                self._samples[d] = {
                    'metrics_csv': None,
                    'web_summary': None,
                }
        except OSError:
            pass
        # Command line file
        self._cmdline_file = os.path.join(self.dir,
                                          "_cmdline")
        if not os.path.exists(self._cmdline_file):
            self._cmdline_file = None
        # Locate multi config.csv file
        if not config_csv:
            try:
                args = self.cmdline.split()
                for idx,arg in enumerate(args):
                    if arg in ('--csv',):
                        config_csv = args[idx+1]
                        break
            except Exception:
                pass
        self._config_csv = config_csv
        # Deal with additional data items
        if not cellranger_exe:
            try:
                cellranger_exe = self.cmdline.split()[0]
            except Exception:
                pass
        self._cellranger_exe = cellranger_exe
        if not reference_data:
            try:
                reference_data = self.config.reference_data_path
            except Exception:
                pass
        self._reference_data = reference_data
        try:
            self._probe_set = self.config.probe_set_path
        except Exception:
            self._probe_set = None
        self._version = version
        # Paths to metrics and web summary files
        for smpl in self.sample_names:
            smpl_dir = os.path.join(self.dir,
                                    "outs",
                                    "per_sample_outs",
                                    smpl)
            self._samples[smpl]['metrics_csv'] = os.path.join(
                smpl_dir,
                "metrics_summary.csv")
            self._samples[smpl]['web_summary'] = os.path.join(
                smpl_dir,
                "web_summary.html")

    @property
    def mode(self):
        return "multi"

    @property
    def dir(self):
        """
        Path to the directory with the 'cellranger multi' outputs
        """
        if os.path.isdir(self._cellranger_multi_dir):
            return self._cellranger_multi_dir
        raise OSError("'%s': not a directory" %
                      self._cellranger_multi_dir)

    @property
    def sample_names(self):
        """
        Sample names derived from the subdirectory names
        """
        return sorted(list(self._samples.keys()))

    @property
    def cellranger_exe(self):
        """
        Cellranger executable
        """
        return self._cellranger_exe

    @property
    def pipeline_name(self):
        """
        Pipeline name i.e. name of the software package
        """
        exe = self.cellranger_exe
        if exe:
            if exe.endswith("bin/count"):
                # Special case: 'count' was run directly
                for name in ("cellranger",
                             "cellranger-atac",
                             "cellranger-arc",):
                    if ("%s-cs" % name) in exe.split(os.sep):
                        return name
            # Default: take name of executable
            return os.path.basename(exe)
        else:
            # Couldn't get the executable
            return None

    @property
    def version(self):
        """
        Cellranger version
        """
        return self._version

    @property
    def reference_data(self):
        """
        Reference dataset
        """
        return self._reference_data

    @property
    def probe_set(self):
        """
        Probe set
        """
        return self._probe_set

    def metrics_csv(self,name):
        """
        Path to the cellranger multi 'metrics.csv' file for a sample
        """
        metrics_csv = self._samples[name]['metrics_csv']
        if os.path.isfile(metrics_csv):
            return metrics_csv
        raise OSError("'%s': not a file" % metrics_csv)

    def web_summary(self,name):
        """
        Path to the cellranger multi 'web_summary.html' file for a sample
        """
        web_summary = self._samples[name]['web_summary']
        if os.path.isfile(web_summary):
            return web_summary
        raise OSError("'%s': not a file" % web_summary)

    def metrics(self,name):
        """
        Return a 'MultiplexSummary' object for a sample
        """
        return MultiplexSummary(self.metrics_csv(name))

    @property
    def cmdline_file(self):
        """
        Path to the 'cellranger multi' '_cmdline' file
        """
        return self._cmdline_file

    @property
    def cmdline(self):
        """
        Return the command line used to run 'cellranger* count'
        """
        if self.cmdline_file:
            with open(self.cmdline_file,'rt') as fp:
                return fp.read().strip()
        else:
            return None

    @property
    def config(self):
        """
        Return CellrangerMultiConfigCsv instance
        """
        if self._config_csv:
            return CellrangerMultiConfigCsv(self._config_csv)
        else:
            return None

#######################################################################
# Functions
#######################################################################

def cellranger_count_output(project,sample_name=None,
                            prefix="cellranger_count"):
    """
    Generate list of 'cellranger count' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/metrics_summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_atac_count_output(project,sample_name=None,
                                 prefix="cellranger_count"):
    """
    Generate list of 'cellranger-atac count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-atac
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_arc_count_output(project,sample_name=None,
                                prefix="cellranger_count"):
    """
    Generate list of 'cellranger-arc count' outputs

    Given an AnalysisProject, the outputs from 'cellranger-arc
    count' will look like:

    - {PREFIX}/{SAMPLE_n}/outs/summary.csv
    - {PREFIX}/{SAMPLE_n}/outs/web_summary.html

    for each SAMPLE_n in the project.

    If a sample name is supplied then outputs are limited
    to those for that sample

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      sample_name (str): sample to limit outputs to
      prefix (str): directory for outputs (defaults
        to "cellranger_count")

    Returns:
       tuple: cellranger count outputs (without leading paths)
    """
    outputs = []
    # Metrics and web summary files
    for sample in project.samples:
        if sample_name and sample_name != sample.name:
            continue
        sample_count_dir = os.path.join(prefix,
                                        sample.name)
        for f in ("summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_count_dir,
                                        "outs",f))
    return tuple(outputs)

def cellranger_multi_output(project,config_csv,sample_name=None,
                            prefix="cellranger_multi"):
    """
    Generate list of 'cellranger multi' outputs

    Given an AnalysisProject, the outputs from 'cellranger
    multi' will look like:

    - {PREFIX}/outs/multi/multiplexing_analysis/tag_calls_summary.csv

    and

    - {PREFIX}/outs/per_sample_outs/{SAMPLE_n}/metrics_summary.csv
    - {PREFIX}/outs/per_sample_outs/{SAMPLE_n}/web_summary.html

    for each multiplexed SAMPLE_n defined in the config.csv file
    (nb these are not equivalent to the 'samples' defined by the
    Fastq files in the project).

    If a sample name is supplied then outputs are limited
    to those for that sample; if the supplied config.csv file isn't
    found then no outputs will be returned.

    Arguments:
      project (AnalysisProject): project to generate
        output names for
      config_csv (str): path to the cellranger multi
        config.csv file
      sample_name (str): multiplexed sample to limit outputs
        to (optional)
      prefix (str): directory for outputs (optional, defaults
        to "cellranger_multi")

    Returns:
       tuple: cellranger multi outputs (without leading paths)
    """
    outputs = []
    # Check that config.csv file exists
    if not os.path.isfile(config_csv):
        return outputs
    # Per-sample metrics and web summary files
    for sample in CellrangerMultiConfigCsv(config_csv).sample_names:
        if sample_name and sample_name != sample:
            continue
        sample_dir = os.path.join(prefix,
                                  "outs",
                                  "per_sample_outs",
                                  sample)
        for f in ("metrics_summary.csv",
                  "web_summary.html"):
            outputs.append(os.path.join(sample_dir,f))
    # Multiplexing outputs
    multi_analysis_dir = os.path.join(prefix,
                                      "outs",
                                      "multi",
                                      "multiplexing_analysis")
    for f in ("tag_calls_summary.csv",):
        outputs.append(os.path.join(multi_analysis_dir,f))
    return tuple(outputs)

def fetch_cellranger_multi_output_dirs(top_dir):
    """
    Locate output directories from cellranger multi

    Recursively searches the directory structure under
    the supplied top-level directory and returns a
    list of paths to each possible "cellranger multi"
    output directory.

    Putative output directories will contain at minimum
    subdirectories called "outs" and "per_sample_outs".

    Arguments:
      top_dir (str): path to directory to search under

    Returns:
       List: list of paths to putative "cellranger multi"
         output directories.
    """
    # Collect all paths ending with "outs"
    outs_dirs = []
    for f in walk(top_dir):
        if os.path.basename(f) == "outs":
            outs_dirs.append(f)
    # Reduce to those which look like cellranger multi outputs
    # i.e. have a "per_sample_outs" subdirectory, and discard
    # "outs" part of path
    multi_dirs = []
    for d in sorted([os.path.dirname(f) for f in outs_dirs
                     if os.path.isdir(os.path.join(f, "per_sample_outs"))]):
        # Prune directories which are under other multi outputs dirs
        prune_dir = False
        for dd in multi_dirs:
            if d.startswith(dd + os.sep):
                prune_dir = True
                break
        if not prune_dir:
            multi_dirs.append(d)
    return multi_dirs
