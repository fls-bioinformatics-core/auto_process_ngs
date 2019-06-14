#!/usr/bin/env python
#
# Handle outputs from cellranger
import os

class CellrangerCount(object):
    """
    Wrapper class for handling outputs from cellranger count

    The ``CellrangerCount`` object gives access to various
    details of the outputs (such as sample name and file
    paths).
    """
    def __init__(self,cellranger_count_dir):
        """
        Create a new CellrangerCount instance

        Arguments:
          cellranger_count_dir (str): path to the cellranger
            count output directory for a sample (e.g.
            ``.../cellranger_count/SAMPLE``)
        """
        self._cellranger_count_dir = os.path.abspath(
            cellranger_count_dir)
        try:
            self._sample_name = os.path.basename(self.dir)
        except OSError:
            self._sample_name = None
        try:
            self._metrics_csv = os.path.join(self.dir,
                                             "outs",
                                             "metrics_summary.csv")
        except OSError:
            self._metrics_csv = None
        try:
            self._web_summary = os.path.join(self.dir,
                                             "outs",
                                             "web_summary.html")
        except OSError:
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
        
        
