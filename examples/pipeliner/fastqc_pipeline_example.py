#!/usr/bin/env python
#
# Example pipeline to run Fastqc on one or more Fastq files
# but ignoring any with zero reads

import os
import argparse
from bcftbx.FASTQFile import nreads
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineFunctionTask
from auto_process_ngs.pipeliner import PipelineCommandWrapper
from auto_process_ngs.pipeliner import Pipeline

class RunFastqc(PipelineTask):
    # Run Fastqc on multiple files
    def init(self,fastqs,out_dir):
        # Inputs:
        # - fastqs: list of input Fastq files
        # - out_dir: where to put the Fastqc outputs
        # Outputs:
        # - files: list of output Fastqc HTML files
        self.add_output('files',list())
    def setup(self):
        if not os.path.exists(self.args.out_dir):
            os.mkdir(self.args.out_dir)
        for fq in self.args.fastqs:
            self.add_cmd(
                PipelineCommandWrapper("Run FastQC",
                                       "fastqc",
                                       "-o",self.args.out_dir,
                                       fq))
    def finish(self):
        for fq in self.args.fastqs:
            if fq.endswith(".gz"):
                fq = os.path.splitext(fq)[0]
            out_file = os.path.join(
                self.args.out_dir,
                os.path.splitext(
                    os.path.basename(fq))[0]+"_fastqc.html")
            if not os.path.exists(out_file):
                self.fail(message="Missing output file: %s" % out_file)
            else:
                self.output.files.append(out_file)

class FilterEmptyFastqs(PipelineFunctionTask):
    # Filter Fastq files based on read count
    def init(self,fastqs):
        self.add_output('fastqs',list())
    def setup(self):
        for fq in self.args.fastqs:
            self.add_call("Filter out empty fastqs",
                          self.filter_empty_fastqs,fq)
    def filter_empty_fastqs(self,*fastqs):
        filtered_fastqs = list()
        for fq in fastqs:
            if nreads(fq) > 0:
                print("%s" % fq)
                filtered_fastqs.append(fq)
        return filtered_fastqs
    def finish(self):
        for result in self.result():
            for fq in result:
                self.output.fastqs.append(fq)

if __name__ == "__main__":
    # Command line
    p = argparse.ArgumentParser()
    p.add_argument("fastqs",nargs='+',metavar="FASTQ")
    args = p.parse_args()

    # Make and run a pipeline
    ppl = Pipeline()
    filter_empty_fastqs = FilterEmptyFastqs("Filter empty Fastqs",
                                            args.fastqs)
    run_fastqc = RunFastqc("Run Fastqc",
                           filter_empty_fastqs.output.fastqs,
                           os.getcwd())
    ppl.add_task(filter_empty_fastqs)
    ppl.add_task(run_fastqc,requires=(filter_empty_fastqs,))
    ppl.run()
    print(run_fastqc.output())
