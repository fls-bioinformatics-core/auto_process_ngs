#!/usr/bin/env python
#
# Example pipeline to run Fastqc on one or more Fastq files
# but ignoring any with zero reads

import os
import argparse
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineCommandWrapper
from auto_process_ngs.pipeliner import Pipeline

class RunFastqc(PipelineTask):
    # Run Fastqc on multiple files
    def init(self,fastqs,out_dir):
        # fastqs: list of input Fastq files
        # out_dir: where to put the Fastqc outputs
        self.out_files = list()
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
                self.out_files.append(out_file)
    def output(self):
        # Returns a list of Fastqc HTML files
        return self.out_files

class CountReads(PipelineTask):
    # Count reads in Fastq files
    def init(self,fastqs):
        # fastqs: list of input Fastq files
        self.counts = dict()
    def setup(self):
        for fq in self.args.fastqs:
            print fq
            if os.path.splitext(fq)[1] == ".gz":
                cat = "zcat"
            else:
                cat = "cat"
            self.add_cmd(
                PipelineCommandWrapper("Count reads",
                                       "echo","-n",fq,"' '","&&",
                                       cat,fq,"|",
                                       "wc","-l"))
    def finish(self):
        for line in self.stdout.split('\n'):
            if not line:
                continue
            fq = line.split()[0]
            read_count = int(line.split()[1])/4
            self.counts[fq] = read_count
    def output(self):
        # Returns a dictionary with file paths as
        # keys mapping to read counts
        return self.counts

class FilterEmptyFastqs(PipelineTask):
    # Filter Fastq files based on read count
    def init(self,read_counts):
        # read_counts: dictionary of read counts
        self.filtered_fastqs = list()
    def setup(self):
        for fq in self.args.read_counts:
            if self.args.read_counts[fq] > 0:
                self.filtered_fastqs.append(fq)
    def output(self):
        # Returns a list of Fastq files
        return self.filtered_fastqs

if __name__ == "__main__":
    # Command line
    p = argparse.ArgumentParser()
    p.add_argument("fastqs",nargs='+',metavar="FASTQ")
    args = p.parse_args()

    # Make and run a pipeline
    ppl = Pipeline()
    read_counts = CountReads("Count the reads",args.fastqs)
    print read_counts.output()
    filter_empty_fastqs = FilterEmptyFastqs("Filter empty Fastqs",
                                            read_counts.output())
    run_fastqc = RunFastqc("Run Fastqc",
                           filter_empty_fastqs.output(),
                           os.getcwd())
    ppl.add_task(read_counts)
    ppl.add_task(filter_empty_fastqs,requires=(read_counts,))
    ppl.add_task(run_fastqc,requires=(filter_empty_fastqs,))
    ppl.run()
    sched.stop()
    print run_fastqc.output()
