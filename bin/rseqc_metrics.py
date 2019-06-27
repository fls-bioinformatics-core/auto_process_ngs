#!/usr/bin/env python
#
# Prototype pipeline for RSeQC functions
#

__version__ = "0.0.1"

#######################################################################
# Imports
#######################################################################

import sys
import os
import argparse
import shutil
import random
import json
import logging
from collections import defaultdict
from bcftbx.utils import find_program
from bcftbx.utils import mkdir
from bcftbx.ngsutils import getreads
from bcftbx.ngsutils import getreads_subset
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineFunctionTask
from auto_process_ngs.pipeliner import PipelineParam
from auto_process_ngs.pipeliner import sanitize_name
from auto_process_ngs.applications import Command
from auto_process_ngs.fastq_utils import pair_fastqs

######################################################################
# Pipeline classes
######################################################################

class JoinDir(PipelineParam):
    """
    """
    def __init__(self,base_dir,subdir):
        PipelineParam.__init__(self,
                               type=self.join_dirs,
                               value=subdir)
        self._base_dir = base_dir
    @property
    def base_dir(self):
        try:
            return self._base_dir.value
        except AttributeError:
            return self._base_dir
    def join_dirs(self,s):
        return os.path.join(self.base_dir,s)

class RSeQCPipeline(Pipeline):
    """
    """
    def __init__(self):
        """
        """
        # Initialise the pipeline superclass
        Pipeline.__init__(self,name="RSeQC metrics")

        # Define parameters
        self.add_param('fastqs',type=list)
        self.add_param('ref_gene_model',type=str)
        self.add_param('genomedir',type=str)
        self.add_param('subset',type=int,value=10000)
        self.add_param('nthreads',type=int,value=1)
        self.add_param('project_name',type=str)
        self.add_param('outdir',type=str)
        self.add_param('working_dir',type=str)

        # Define runners
        self.add_runner('star')
        self.add_runner('rseqc')

        ####################
        # Build the pipeline
        ####################
        
        # Pair and sort Fastqs
        pair_and_sort_fastqs = PairAndSortFastqs(
            "Pair and sort Fastqs",
            self.params.fastqs)
        self.add_task(pair_and_sort_fastqs,
                      runner=self.runners['default'])

        # Generate subsetted Fastqs
        subset_fastqs = SubsetFastqs(
            "Subset Fastq pairs",
            self.params.subset,
            pair_and_sort_fastqs.output.fastq_pairs,
            JoinDir(self.params.working_dir,"_fastqs.subset")
        )
        self.add_task(subset_fastqs,
                      requires=(pair_and_sort_fastqs,),
                      runner=self.runners['default'])

        # Make BAM files
        make_bams = MakeBamFiles(
            "Generate and sort BAM files",
            subset_fastqs.output.fastq_pairs,
            self.params.genomedir,
            self.params.nthreads,
            JoinDir(self.params.working_dir,"_bam_files")
        )
        self.add_task(make_bams,
                      requires=(subset_fastqs,),
                      runner=self.runners['star'])

        # Run RSeQC gene body coverage
        rseqc_gene_body_coverage = RSeQCGeneBodyCoverage(
            "Run RSeQC gene body coverage",
            make_bams.output.bam_files,
            self.params.ref_gene_model,
            self.params.outdir,
            self.params.working_dir,
            self.params.project_name
        )
        self.add_task(rseqc_gene_body_coverage,
                      requires=(make_bams,),
                      runner=self.runners['rseqc'])

        # Run RSeQC inner distance
        rseqc_inner_distance = RSeQCInnerDistance(
            "Run RSeQC inner distance",
            make_bams.output.bam_files,
            self.params.ref_gene_model,
            self.params.outdir,
            self.params.working_dir
        )
        self.add_task(rseqc_inner_distance,
                      requires=(make_bams,),
                      runner=self.runners['rseqc'])

    def run(self,fastqs,ref_gene_model,genomedir,
            subset=10000,nthreads=1,outdir=None,
            project_name=None,working_dir=None,
            runners=None):
        """
        """
        # Top-level output dir
        if outdir is None:
            outdir = sanitize_name(self._name)
        outdir = os.path.abspath(outdir)

        # Working directory
        if working_dir is None:
            working_dir = "__%s.ppl.out" % os.path.basename(outdir)
        working_dir = os.path.abspath(working_dir)

        # Log and script directories
        log_dir = os.path.join(working_dir,"logs")
        scripts_dir = os.path.join(working_dir,"scripts")

        # Check parameters
        params = dict(fastqs=fastqs,
                      ref_gene_model=ref_gene_model,
                      genomedir=genomedir,
                      subset=subset,
                      nthreads=nthreads,
                      outdir=outdir,
                      project_name=project_name)
        params_file = os.path.join(working_dir,".params.json")
        if os.path.exists(params_file):
            with open(params_file,'r') as fp:
                stored_params = json.load(fp)
                if params != stored_params:
                    raise Exception("FATAL: trying to rerun pipeline "
                                    "with different parameters")
                else:
                    self.report("Rerunning pipeline: parameters "
                                "are consistent")

        # Check for required software
        if self.check_software('STAR',
                               'samtools',
                               'geneBody_coverage.py',
                               'inner_distance.py',
                               'Rscript'):
            raise Exception("Missing required software")

        # Report the runners
        for r in self.runners:
            self.report("Runner '%s': %s" % (r,self.runners[r]))

        # Make directories
        if not os.path.exists(working_dir):
            mkdir(working_dir)
        if not os.path.exists(outdir):
            mkdir(outdir)

        # Store parameters
        self.report("Storing pipeline parameters")
        with open(params_file,'w') as fp:
            json.dump(params,fp)

        # Execute the pipeline
        return Pipeline.run(self,
                            working_dir=working_dir,
                            log_dir=log_dir,
                            scripts_dir=scripts_dir,
                            params={
                                'fastqs': [fq for fq in fastqs],
                                'ref_gene_model': ref_gene_model,
                                'genomedir': genomedir,
                                'subset': subset,
                                'nthreads': nthreads,
                                'project_name': project_name,
                                'outdir': outdir,
                                'working_dir': working_dir,
                            },
                            runners=runners)

    def check_software(self,*progs):
        """
        """
        # Check for required software
        self.report("Checking for required software")
        missing = []
        for prog in progs:
            exe = find_program(prog)
            if exe:
                self.report("* %s ... %s" % (prog,exe))
            else:
                self.report("* %s ... not found" % prog)
                missing.append(prog)
        return missing

######################################################################
# Pipeline task classes
######################################################################

class PairAndSortFastqs(PipelineFunctionTask):
    """
    Arrange Fastqs into sorted R1/R2 pairs

    This is a utility task that can be used
    to arrange a list of Fastq files into
    R1/R2 pairs according to their contents.
    Essentially it is a wrapper for the
    'pair_fastqs' function.
    """
    def init(self,fastqs):
        """
        Initialise the PairAndSortFastqs task

        Arguments:
          fastqs (list): list of paths to
            the Fastq files to be paired

        Outputs:
          fastq_pairs: list of tuples with
            R1/R2 Fastq pairs
        """
        self.add_output('fastq_pairs',list())
    def setup(self):
        self.add_call("Pair and sort Fastq files",
                      self.pair_and_sort_fastqs,
                      self.args.fastqs)
    def pair_and_sort_fastqs(self,fqs):
        # Do initial pairing
        paired,unpaired = pair_fastqs(fqs)
        # Combine into a single list
        paired.extend([(fq,) for fq in unpaired])
        return sorted(paired,cmp=lambda x,y: cmp(x[0],y[0]))
    def finish(self):
        # 'pair_and_sort_fastqs' returns a list of
        # of Fastq 'pairs'
        for pair in self.result()[0]:
            if len(pair) == 1:
                print "-- %s" % pair[0]
            elif len(pair) == 2:
                print "-- R1: %s" % pair[0]
                print "   R2: %s" % pair[1]
            else:
                print "Bad Fastq 'pair': %s" % pair
                self.fail("Sorting/pairing failed")
            self.output.fastq_pairs.append(pair)

class SubsetFastqs(PipelineFunctionTask):
    """
    Create Fastqs with subset of reads
    """
    def init(self,subset_size,fastq_pairs,out_dir):
        """
        Initialise the SubsetFastqs task

        Arguments:
          r1 (str): path to the R1 Fastq
          r2 (str): path to the R2 Fastq
          subset (int): number of reads in subset
          out_dir (str): path to directory
            to write the final Fastqs to

        Outputs:
          fastqs_pairs: list of tuples with subsetted
            R1/R2 Fastq pairs
        """
        self.add_output('fastq_pairs',list())
    def setup(self):
        if os.path.exists(self.args.out_dir):
            # Outputs have already been generated
            print "%s already exists" % self.args.out_dir
            # Set the outputs
            for pair in self.args.fastq_pairs:
                fqs = []
                for fq in pair:
                    if fq is not None:
                        fq = os.path.join(self.args.out_dir,
                                          os.path.basename(fq))
                        if fq.endswith('.gz'):
                            fq = fq[:-3]
                    fqs.append(fq)
                self.output.fastq_pairs.append(tuple(fqs))
            return
        # Set up temporary directory
        self.tmp_out_dir = tmp_dir(self.args.out_dir)
        print "Made temp dir %s" % self.tmp_out_dir
        # Set up the subsetting jobs
        for pair in self.args.fastq_pairs:
            self.add_call("Subset Fastqs",
                          self.subset_fastqs,
                          self.args.subset_size,
                          pair[0],pair[1],
                          self.tmp_out_dir)
    def subset_fastqs(self,size,r1,r2,out_dir):
        """
        Make Fastqs with subset of input reads or read pairs
        """
        # Determine number of reads from input Fastqs
        nreads = sum(1 for i in getreads(os.path.abspath(r1)))
        if size > nreads:
            print "Number of reads smaller than requested subset"
            size = 0
        if size == 0:
            print "Using all reads/read pairs in Fastq file(s)"
            size = nreads
        else:
            print "Using random subset of %d reads/read pairs" \
                % size
        # Generate subset indices to extract
        if size == nreads:
            subset_indices = [i for i in xrange(nreads)]
        else:
            subset_indices = random.sample(xrange(nreads),size)
        # Do the subsetting
        fqs_in = filter(lambda fq: fq is not None,(r1,r2))
        fastqs = []
        for fq in fqs_in:
            fq_subset = os.path.join(out_dir,os.path.basename(fq))
            if fq_subset.endswith(".gz"):
                fq_subset = '.'.join(fq_subset.split('.')[:-1])
            with open(fq_subset,'w') as fp:
                for read in getreads_subset(os.path.abspath(fq),
                                            subset_indices):
                    fp.write('\n'.join(read) + '\n')
                fastqs.append(fq_subset)
        # Return the subsetted Fastq files
        return fastqs
    def finish(self):
        # Output dir already exists
        if os.path.exists(self.args.out_dir):
            return
        # Move the tmp dir to the final location
        os.rename(self.tmp_out_dir,self.args.out_dir)
        # 'subset_fastqs' returns pairs of subsetted Fastqs
        # in the temporary directory
        # Update the list so that Fastqs have the correct path
        for result in self.result():
            pair = []
            for fq in result:
                if fq is not None:
                    fq = os.path.join(self.args.out_dir,
                                      os.path.basename(fq))
                    if not os.path.exists(fq):
                        self.fail("Missing Fastq: %s not found" %
                                  fq)
                pair.append(fq)
            self.output.fastq_pairs.append(tuple(pair))
        print self.output

class MakeBamFiles(PipelineFunctionTask):
    """
    Make BAM files from Fastqs using STAR
    """
    def init(self,fastqs,genomedir,nthreads,out_dir):
        """
        Initialise the MakeBamFiles task

        Arguments:
          fastqs (list): either R1 Fastq or R1/R2
            Fastq pair
          genomedir (str): path to directory with
            STAR genome index to map against
          nthreads (int): number of cores for STAR
            to use
          out_dir (str): path to directory
            to write the final BAM files to

        Outputs:
          bam_files: list of sorted BAM files
        """
        self.add_output('bam_files',list())
    def setup(self):
        if os.path.exists(self.args.out_dir):
            # Outputs have already been generated
            print "%s already exists" % self.args.out_dir
            # Set the outputs
            for pair in self.args.fastqs:
                bam_file = os.path.join(self.args.out_dir,
                                        "%s.bam" %
                                        strip_extension(
                                            os.path.basename(pair[0])))
                self.output.bam_files.append(bam_file)
            return
        self.tmp_out_dir = tmp_dir(self.args.out_dir)
        print "Made temp dir %s" % self.tmp_out_dir
        for fastq_pair in self.args.fastqs:
            self.add_call("Make BAM files",
                          self.make_bam_files,
                          fastq_pair,
                          self.args.genomedir,
                          self.args.nthreads,
                          self.tmp_out_dir)
    def make_bam_files(self,fastqs,genomedir,nthreads,out_dir):
        # Generate and sort the BAM files
        # Basename and prefix for output files
        basename = strip_extension(os.path.basename(fastqs[0]))
        prefix = "%s_" % basename
        # Working directory for STAR outputs
        working_dir = os.path.join(out_dir,"STAR_%s" % basename)
        print "Creating directory for STAR outputs: %s" % working_dir
        mkdir(working_dir)
        # Output BAM file name will be "<prefix>Aligned.out.bam"
        star_bam_file = os.path.join(working_dir,
                                     "%sAligned.out.bam" % prefix)
        # Final BAM file will be "<FASTQ>.bam"
        bam_file = os.path.join(out_dir,"%s.bam" % basename)
        if os.path.exists(bam_file):
            print "BAM file '%s' already exists" % bam_file
            return bam_file
        # Check the Fastqs
        if fastqs[0].endswith('.gz'):
            print "Input Fastq(s) are gzipped, making uncompressed versions"
            uncompressed_fastqs = []
            for fq in fastqs:
                gunzip_cmd = Command('gzip','-d','-c',fq)
                fq_gunzip = os.path.join(out_dir,
                                         os.path.basename(fq[0:-3]))
                print "Running %s" % ' '.join(gunzip_cmd)
                status = gunzip_cmd.run_subprocess(working_dir=working_dir,
                                                   log=fq_gunzip)
                if status != 0:
                    raise Exception("Error uncompressing %s: exit code: %s"
                                    % (fq,status))
                uncompressed_fastqs.append(fq_gunzip)
            fastqs = uncompressed_fastqs
        # Build the STAR command line for mapping
        star_cmd = Command('STAR',
                           '--runMode','alignReads',
                           '--genomeLoad','NoSharedMemory',
                           '--genomeDir',os.path.abspath(genomedir),
                           '--readFilesIn',fastqs[0])
        if len(fastqs) > 1:
            star_cmd.add_args(fastqs[1])
        star_cmd.add_args('--outSAMtype','BAM','Unsorted',
                          '--outSAMstrandField','intronMotif',
                          '--outFileNamePrefix',prefix,
                          '--runThreadN',nthreads)
        print "Running %s" % ' '.join(star_cmd)
        status = star_cmd.run_subprocess(working_dir=working_dir)
        if status != 0:
            raise Exception("STAR returned non-zero exit code: %s"
                            % status)
        # Sort the BAM file
        sorted_bam_file_prefix = os.path.join(working_dir,
                                              "%s.sorted" %
                                              os.path.basename(star_bam_file)[:-4])
        # Run the sorting
        samtools_sort_cmd = Command('samtools',
                                    'sort',
                                    star_bam_file,
                                    sorted_bam_file_prefix)
        print "Running %s" % ' '.join(samtools_sort_cmd)
        status = samtools_sort_cmd.run_subprocess(working_dir=working_dir)
        if status != 0:
            raise Exception("samtools sort returned non-zero exit code: %s" %
                            status)
        sorted_bam_file = "%s.bam" % sorted_bam_file_prefix
        # Index the sorted BAM file
        samtools_index_cmd = Command('samtools',
                                     'index',
                                     sorted_bam_file)
        print "Running %s" % ' '.join(samtools_index_cmd)
        status = samtools_index_cmd.run_subprocess(working_dir=working_dir)
        if status != 0:
            raise Exception("samtools index returned non-zero exit code: %s" %
                            status)
        # Output BAM index file (.bai)
        bai_file = "%s.bai" % sorted_bam_file
        # Move the BAM and BAI files
        os.rename(sorted_bam_file,bam_file)
        os.rename(bai_file,os.path.join(out_dir,"%s.bai" % bam_file))
        # Remove the STAR working directory
        shutil.rmtree(working_dir)
        # Return the BAM file name
        return bam_file
    def finish(self):
        # Output dir already exists
        if os.path.exists(self.args.out_dir):
            return
        # Move the tmp dir to the final location
        os.rename(self.tmp_out_dir,self.args.out_dir)
        # 'make_bam_files' returns the path to a BAM file
        # in the temporary directory
        # Update the paths to point to the final directory
        for bam_file in self.result():
            bam_file = os.path.join(self.args.out_dir,
                                    os.path.basename(bam_file))
            if not os.path.exists(bam_file):
                self.fail("Missing BAM file: %s not found" %
                          bam_file)
            self.output.bam_files.append(bam_file)
        print self.output

class RSeQCGeneBodyCoverage(PipelineFunctionTask):
    """
    Run the geneBody_coverage utility from RSeQC
    """
    def init(self,bam_files,reference_gene_model,out_dir,
             working_dir=None,project_name=None):
        """
        Initialise the RSeQCGeneBodyCoverage task

        Arguments:
          bam_files (list): list of paths to the
            input BAM files
          reference_gene_model (str): path to
            BED file with reference gene model
          out_dir (str): path to directory
            to write the gene body coverage files
          working_dir (str): path to working
            directory (or 'None')
          project_name (str): either a base name
            to use for the output files, or
            'None'

        Outputs:
          out_files: list of output files from
            geneBody_coverage.py
        """
        self.add_output('out_files',list())
    def setup(self):
        if not self.args.bam_files:
            print "No BAM files supplied: nothing to do"
            return
        self.tmp_out_dir = tmp_dir("%s.gene_body_coverage"
                                   % self.args.out_dir)
        if self.args.project_name:
            self.project_name = str(self.args.project_name)
        else:
            self.project_name = self.basename(self.args.bam_files[0])
        print "Made temp dir %s" % self.tmp_out_dir
        for f in self.rseqc_genebody_coverage_files(self.project_name,
                                                    self.args.out_dir):
            if not os.path.exists(f):
                self.add_call("RSeQC gene body coverage",
                              self.rseqc_genebody_coverage,
                              self.args.bam_files,
                              self.args.reference_gene_model,
                              self.project_name,
                              self.tmp_out_dir)
                break
    def basename(self,f):
        return strip_extension(os.path.basename(f))
    def rseqc_genebody_coverage_files(self,prefix,dirn=None):
        # Outputs are:
        # -- <prefix>.geneBodyCoverage.txt
        # -- <prefix>.geneBodyCoverage.r
        # -- log.txt
        # -- <prefix>.geneBodyCoverage.curves.png
        # -- <prefix>.geneBodyCoverage.heatMap.png
        outputs = ["%s.geneBodyCoverage.txt" % prefix,
                   "%s.geneBodyCoverage.r" % prefix,
                   "%s.geneBodyCoverage.curves.png" % prefix]
        if len(self.args.bam_files) > 2:
            outputs.append("%s.geneBodyCoverage.heatMap.png" % prefix)
        if dirn is not None:
            outputs = [os.path.join(dirn,f) for f in outputs]
        return tuple(outputs)
    def rseqc_genebody_coverage(self,bam_files,reference_gene_model,
                                prefix,out_dir):
        """
        Run RSeQC 'genebody_coverage.py'
        """
        # Build the command line for genebody_coverage.py
        genebody_coverage_cmd = Command('geneBody_coverage.py',
                                        '-i',','.join(bam_files),
                                        '-r',reference_gene_model,
                                        '-f','png',
                                        '-o',prefix)
        print "Running %s" % genebody_coverage_cmd
        status = genebody_coverage_cmd.run_subprocess(working_dir=out_dir)
        if status != 0:
            raise Exception("geneBody_coverage.py returned non-zero exit "
                            "code: %s" % status)
        # Print the contents of the out_dir
        print "Contents of %s" % out_dir
        for f in os.listdir(out_dir):
            print "-- %s" % f
        # Return list of outputs
        return self.rseqc_genebody_coverage_files(prefix,out_dir)
    def finish(self):
        # Store the outputs
        for f in self.rseqc_genebody_coverage_files(self.project_name,
                                                    self.tmp_out_dir):
            ff = os.path.join(self.args.out_dir,os.path.basename(f))
            if os.path.exists(f):
                os.rename(f,ff)
            if not os.path.exists(ff):
                self.fail(message="Missing file: %s" % ff)
            self.output.out_files.append(ff)
        # Remove the temporary directory
        shutil.rmtree(self.tmp_out_dir)

class RSeQCInnerDistance(PipelineFunctionTask):
    """
    Run the inner_distance utility from RSeQC
    """
    def init(self,bam_files,reference_gene_model,out_dir,
             working_dir=None):
        """
        Initialise the RSeQCInnerDistance task

        Arguments:
          bam_files (list): list of paths to the
            input BAM files
          reference_gene_model (str): path to
            BED file with reference gene model
          out_dir (str): path to directory
            to write the gene body coverage files
          working_dir (str): path to working
            directory (or 'None')

        Outputs:
          out_files: dictionary with BAM files as
            keys and list of associated output files
            from inner_distance.py as values
        """
        self.add_output('out_files',dict())
    def setup(self):
        if not self.args.bam_files:
            print "No BAM files supplied: nothing to do"
            return
        self.tmp_out_dir = tmp_dir("%s.inner_distance"
                                   % self.args.out_dir)
        print "Made temp dir %s" % self.tmp_out_dir
        for bam_file in self.args.bam_files:
            for f in self.rseqc_inner_distance_files(bam_file,
                                                     self.args.out_dir):
                if not os.path.exists(f):
                    self.add_call("RSeQC inner distance",
                                  self.rseqc_inner_distance,
                                  bam_file,
                                  self.args.reference_gene_model,
                                  self.tmp_out_dir)
                    break
    def basename(self,f):
        return strip_extension(os.path.basename(f))
    def rseqc_inner_distance_files(self,bam_file,dirn=None):
        # Outputs are:
        # -- <prefix>.inner_distance.txt
        # -- <prefix>.inner_distance_plot.r
        # -- <prefix>.inner_distance_freq.txt
        # -- <prefix>.inner_distance_plot.pdf
        prefix = self.basename(bam_file)
        outputs = ["%s.inner_distance.txt" % prefix,
                   "%s.inner_distance_plot.r" % prefix,
                   "%s.inner_distance_plot_png.r" % prefix,
                   "%s.inner_distance_freq.txt" % prefix,
                   "%s.inner_distance_plot.pdf" % prefix,
                   "%s.inner_distance_plot.png" % prefix]
        if dirn is not None:
            outputs = [os.path.join(dirn,f) for f in outputs]
        return tuple(outputs)
    def rseqc_inner_distance(self,bam_file,reference_gene_model,
                             out_dir):
        """
        Run RSeQC 'inner_distance.py'
        """
        # Check if BAM file is paired-end
        # See https://www.biostars.org/p/178730/#178732
        samtools_view_cmd = Command('samtools',
                                    'view',
                                    '-c','-f',1,
                                    bam_file)
        print "Running %s" % samtools_view_cmd
        logfile = os.path.join(out_dir,
                               "%s.pe_count" %
                               os.path.basename(bam_file))
        status = samtools_view_cmd.run_subprocess(working_dir=out_dir,
                                                  log=logfile)
        with open(logfile,'r') as fp:
            npairs = int(fp.read().strip())
            if not npairs:
                # Single-ended so don't run inner_distance.py
                print "%s: single-ended reads" % bam_file
                return (bam_file,)
        # Output prefix
        prefix = os.path.join(out_dir,self.basename(bam_file))
        # Build the command line for inner_distance.py
        inner_distance_cmd = Command('inner_distance.py',
                                     '-i',bam_file,
                                     '-r',reference_gene_model,
                                     '-o',prefix)
        print "Running %s" % inner_distance_cmd
        status = inner_distance_cmd.run_subprocess(working_dir=out_dir)
        if status != 0:
            raise Exception("inner_distance.py returned non-zero exit "
                            "code: %s" % status)
        # Make modified version of R script for PNG
        rscript = "%s.inner_distance_plot.r" % prefix
        rscript_png = "%s.inner_distance_plot_png.r" % prefix
        with open(rscript,'r') as frpdf:
            with open(rscript_png,'w') as frpng:
                frpng.write(frpdf.read().replace("pdf","png"))
        status = Command('Rscript',
                         '--vanilla',
                         rscript_png).run_subprocess(working_dir=out_dir)
        if status != 0:
            raise Exception("Failed to generate PNG of inner distance plot: "
                            "%s" % status)
        return (bam_file,self.rseqc_inner_distance_files(bam_file,out_dir))
    def finish(self):
        # Store the outputs
        for bam_file in self.args.bam_files:
            self.output.out_files[bam_file] = list()
            for f in self.rseqc_inner_distance_files(bam_file,
                                                     self.tmp_out_dir):
                ff = os.path.join(self.args.out_dir,
                                  os.path.basename(f))
                if os.path.exists(f):
                    os.rename(f,ff)
                if not os.path.exists(ff):
                    self.fail(message="Missing file: %s" % f)
                self.output.out_files[bam_file].append(ff)
        # Remove the temporary directory
        shutil.rmtree(self.tmp_out_dir)

######################################################################
# Functions
######################################################################

def tmp_dir(d):
    """
    Create a temp dir for directory 'd'

    **Copied from icell8/pipeline.py**
    """
    # Make temp directory for outputs
    tmp = "%s.tmp" % d
    if os.path.exists(tmp):
        print "Removing existing tmp dir '%s'" % tmp
        shutil.rmtree(tmp)
    print "Creating tmp dir '%s'" % tmp
    mkdir(tmp)
    return tmp

def strip_extension(f):
    """
    Remove the extension from a file name
    """
    while f.split(".")[-1] in ("fastq","bam","gz"):
        f = ".".join(f.split(".")[:-1])
    return f

def rseqc_mean_insert_size(freq_txt):
    """
    Calculate the mean and SD insert size

    Calculates the mean and standard deviation of the
    inner distance/insert size directly from the
    frequencies output by RSeQC inner_distance.py.
    """
    try:
        frequency_data = {}
        with open(freq_txt) as fp:
            for line in fp:
                values = [int(i) for i in line.rstrip('\n').split('\t')]
                distance = int((values[0] + values[1])/2)
                frequency_data[distance] = values[2]
        # Total number of numbers
        n = sum([frequency_data[i] for i in frequency_data])
        # Mean
        mean = float(sum([i*frequency_data[i] for i in frequency_data]))/float(n)
        # SD
        sd = (sum([((i - mean)**2)*frequency_data[i] for i in frequency_data])/float(n))**0.5
        print "Mean = %f SD = %f" % (mean,sd)
        return (mean,sd)
    except Exception as ex:
        raise Exception("Error trying to get stats on inner distances: "
                        "%s" % ex)

######################################################################
# Command line interface
######################################################################

if __name__ == "__main__":
    
    # Process command line
    p = argparse.ArgumentParser(version=__version__)
    p.add_argument("fastqs",metavar="FASTQ",
                   default=None,
                   nargs="+",
                   help="input Fastq files")
    p.add_argument("-g","--genomedir",
                   dest="genomedir",metavar="GENOMEDIR",
                   default=None,
                   help="path to directory with STAR index "
                   "for genome to use")
    p.add_argument("-r","--ref_gene_model",
                   metavar="REF_GENE_MODEL",
                   default=None,
                   help="path to reference gene model file "
                   "(BED format)")
    p.add_argument("-s","--subset",
                   type=int,
                   default=10000,
                   help="use a random subset of read pairs "
                   "from the input Fastqs; set to zero to "
                   "use all reads (default: 10000)")
    p.add_argument("-n","--nthreads",
                   dest="nthreads",
                   default=1,type=int,
                   help="number of threads to use when "
                   "running STAR")
    p.add_argument("-o","--outdir",
                   metavar="OUTPUT_DIR",
                   default="rseqc_metrics.out",
                   help="path to output directory (default: "
                   "<CWD>/rseqc_metrics.out)")
    p.add_argument("-p","--project_name",
                   metavar="NAME",
                   default="RSEQC",
                   help="basename for gene body coverage "
                   "outputs (default: 'RSEQC')")
                   
    args = p.parse_args()

    # Input Fastqs
    fastqs = [os.path.abspath(fq) for fq in args.fastqs]

    # Reference genes BED file
    ref_gene_model = os.path.abspath(args.ref_gene_model)

    # STAR genome dir
    genomedir = os.path.abspath(args.genomedir)

    # Set up and run the pipeline
    try:
        ppl = RSeQCPipeline()
        retval = ppl.run(fastqs,
                         ref_gene_model,
                         genomedir,
                         subset=args.subset,
                         nthreads=args.nthreads,
                         project_name=args.project_name,
                         outdir=args.outdir)
    except Exception as ex:
        logging.fatal("RSeQC pipeline failed: %s" % ex)
        retval = 1

    # Finish
    sys.exit(retval)
    
