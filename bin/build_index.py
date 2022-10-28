#!/usr/bin/env python
#
#     build_index.py: utility to build indexes for aligners
#     Copyright (C) University of Manchester 2022 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import argparse
from bcftbx.JobRunner import fetch_runner
from auto_process_ngs.indexes import IndexBuilder
from auto_process_ngs.settings import Settings
from auto_process_ngs import get_version

__version__ = get_version()
__settings = Settings()

#######################################################################
# Constants
#######################################################################

_ALIGNERS = {
    'bowtie': '1.0.0',
    'bowtie2': '2.4.1',
    'star': '2.7.7a'
}

#######################################################################
# Main program
#######################################################################

def main():
    # List of aligner names
    aligners = tuple(sorted(_ALIGNERS.keys()))
    # Create parser to handle command line
    p = argparse.ArgumentParser(description="Generate indexes for "
                                "aligners")
    p.add_argument('-v','--version',action='version',
                   version="%(prog)s "+__version__)
    p.add_argument('aligner',metavar="ALIGNER",
                   choices=list(_ALIGNERS),
                   help="aligner to build index for (one of %s)" %
                   ', '.join(["'%s'" % aligner
                              for aligner in _ALIGNERS]))
    p.add_argument('fasta',metavar="FASTA",
                   help="FASTA file with sequence")
    p.add_argument('annotation',metavar="ANNOTATION",nargs='?',
                   help="annotation file (for use with STAR)")
    p.add_argument('-o',dest='out_dir',action='store',
                   help="output directory for indexes")
    # Bowtie-specific options
    bowtie = p.add_argument_group('Bowtie-specific options')
    bowtie.add_argument("--ebwt_base",metavar='NAME',action='store',
                        dest="ebwt_base",default=None,
                        help="specify basename for output .ebwt files "
                        "(defaults to FASTA file basename)")
    # Bowtie2-specific options
    bowtie2 = p.add_argument_group('Bowtie-specific options')
    bowtie2.add_argument("--bt2_base",metavar='NAME',action='store',
                         dest="bt2_base",default=None,
                         help="specify basename for output .bt2 files "
                         "(defaults to FASTA file basename)")
    # STAR-specific options
    star = p.add_argument_group('STAR-specific options')
    star.add_argument("--overhang",metavar='N',action='store',
                      dest="overhang",default=100,
                      help="set value for STAR --sjdbOverhang "
                      "option (default: 100)")
    # Advanced options
    advanced = p.add_argument_group('Advanced options')
    advanced.add_argument('-V','--aligner-version',metavar='VERSION',
                          action='store',
                          dest="aligner_version",default=None,
                          help="specify the version of the aligner to "
                          "target (only works if conda dependency "
                          "resolution is configured)")
    advanced.add_argument('-r','--runner',metavar='RUNNER',action='store',
                          dest="runner",default=None,
                          help="explicitly specify runner definition for "
                          "building the index. RUNNER must be a valid job "
                          "runner specification e.g. "
                          "'GEJobRunner(-pe smp.pe 8)' (default: use "
                          "appropriate runner from configuration)")
    p.add_argument_group()
    args = p.parse_args()

    # Aligner
    aligner = args.aligner.lower()

    # Aligner version
    aligner_version = args.aligner_version
    if aligner_version is None:
        aligner_version = _ALIGNERS[aligner]
    print("Requested %s version: %s" % (aligner,aligner_version))

    # Acquire runner
    if args.runner:
        # Use runner supplied on command line
        runner = fetch_runner(args.runner)
    else:
        if aligner in ("bowtie","bowtie2"):
            runner = __settings.runners.fastq_screen
        elif aligner == "star":
            runner = __settings.runners.star
    print("Job runner: %s" % runner)

    # Set up index builder
    builder = IndexBuilder(runner,
                           use_conda=Settings().conda.enable_conda,
                           conda_env_dir=Settings().conda.env_dir)

    # Build indexes
    if aligner == "bowtie":
        builder.bowtie(args.fasta,args.out_dir,
                       ebwt_basename=args.ebwt_base,
                       bowtie_version=aligner_version)
    elif aligner == "bowtie2":
        builder.bowtie2(args.fasta,args.out_dir,
                        bt2_basename=args.bt2_base,
                        bowtie2_version=aligner_version)
    elif aligner == "star":
        builder.STAR(args.fasta,args.annotation,args.out_dir,
                     overhang=args.overhang,
                     star_version=aligner_version)

if __name__ == '__main__':
    main()
