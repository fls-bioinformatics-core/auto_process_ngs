"""Description

Setup script to install auto_process_ngs

Copyright (C) University of Manchester 2013-20 Peter Briggs

"""

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = []
for pattern in ('bin/*.py','bin/*.sh',):
    scripts.extend(glob(pattern))

# Installation requirements
install_requires = ['configparser',
                    'pillow',
                    'pandas',
                    'cloudpickle',
                    'psutil',
                    'future',
                    'genomics-bcftbx']
# Handle matplotlib for Python2/3
import sys
if sys.version_info < (3,):
    install_requires.append('matplotlib<=2.2.3')
else:
    install_requires.append('matplotlib')

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
                  'auto_process_ngs.cli',
                  'auto_process_ngs.commands',
                  'auto_process_ngs.icell8',
                  'auto_process_ngs.qc',],
      license = 'Artistic License',
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
      include_package_data=True,
      zip_safe = False)
