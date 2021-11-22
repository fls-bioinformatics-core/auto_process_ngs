"""Description

Setup script to install auto_process_ngs

Copyright (C) University of Manchester 2013-2021 Peter Briggs

"""

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = []
for pattern in ('bin/*.py','bin/*.sh',):
    scripts.extend(glob(pattern))

# Installation requirements
install_requires = ['cloudpickle',
                    'configparser',
                    'genomics-bcftbx',
                    'matplotlib',
                    'pandas',
                    'pillow',
                    'psutil']

# If we're on ReadTheDocs then we can reduce this
# to a smaller set (to avoid build timeouts)
import os
if os.environ.get('READTHEDOCS') == 'True':
    install_requires = []

# Setup for installation etc
from setuptools import setup
import auto_process_ngs
setup(name = "auto_process_ngs",
      version = auto_process_ngs.get_version(),
      description = 'Automated processing of NGS data from Illumina sequencers',
      long_description = """Utilities to help with the automated processing,
      QC and management of data from Illumina Next Generation Sequencing (NGS)
      platforms""",
      url = 'https://github.com/fls-bioinformatics-core/auto_process_ngs',
      maintainer = 'Peter Briggs',
      maintainer_email = 'peter.briggs@manchester.ac.uk',
      packages = ['auto_process_ngs',
                  'auto_process_ngs.barcodes',
                  'auto_process_ngs.bcl2fastq',
                  'auto_process_ngs.cli',
                  'auto_process_ngs.commands',
                  'auto_process_ngs.icell8',
                  'auto_process_ngs.qc',],
      license = 'AFL-3',
      # Pull in dependencies
      install_requires = install_requires,
      # Enable 'python setup.py test'
      test_suite='nose.collector',
      tests_require=['nose'],
      # Scripts
      scripts = scripts,
      # Configuration file for QC
      data_files = [('config',['config/auto_process.ini.sample']),
                    ('templates',['templates/README.template.sample']),
                    ('etc/bash_completion.d',['scripts/auto_process-completion.bash']),],
      classifiers=[
          "Development Status :: 4 - Beta",
          "Environment :: Console",
          "Intended Audience :: End Users/Desktop",
          "Intended Audience :: Science/Research",
          "Intended Audience :: Developers",
          "License :: OSI Approved :: Academic Free License (AFL)",
          "Operating System :: POSIX :: Linux",
          "Operating System :: MacOS",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Programming Language :: Python :: 3",
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
      ],
      include_package_data=True,
      zip_safe = False)
