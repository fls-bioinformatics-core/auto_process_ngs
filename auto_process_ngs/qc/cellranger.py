#!/usr/bin/env python
#
# Handle outputs from cellranger
import os
from ..tenx_genomics_utils import GexSummary
from ..tenx_genomics_utils import AtacSummary
from ..tenx_genomics_utils import MultiplexSummary
from ..tenx_genomics_utils import MultiomeSummary
from ..tenx_genomics_utils import CellrangerMultiConfigCsv

class CellrangerCount(object):
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
        if self.pipeline_name in ('cellranger',):
            metrics_csv = "metrics_summary.csv"
        elif self.pipeline_name in ('cellranger-atac',
                                    'cellranger-arc'):
            metrics_csv = "summary.csv"
        else:
            metrics_csv = None
        try:
            self._metrics_csv = os.path.join(self.dir,
                                             "outs",
                                             metrics_csv)
            if not os.path.exists(self._metrics_csv):
                self._metrics_csv = None
        except OSError:
            # OSError raise by self.dir if top-level is
            # missing
            self._metrics_csv = None
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

class CellrangerMulti(object):
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
