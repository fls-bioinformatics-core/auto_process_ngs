"""Description

Setup script to install auto_process_ngs

Copyright (C) University of Manchester 2013-15 Peter Briggs

"""

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = []
for pattern in ('bin/*.py','bin/*.sh',):
    scripts.extend(glob(pattern))

# Setup for installation etc
from setuptools import setup
import auto_process_ngs
setup(name = "auto_process_ngs",
      version = auto_process_ngs.get_version(),
      description = 'Automatic processing of NGS sequencing data from Illumina sequencers',
      long_description = """Utilities to help with the automated processing and management
      of sequencing data from Illumina’s HiSEQ and MiSEQ platforms""",
      url = 'https://bitbucket.org/pjbriggs/auto_process_ngs',
      maintainer = 'Peter Briggs',
      maintainer_email = 'peter.briggs@manchester.ac.uk',
      packages = ['auto_process_ngs'],
      license = 'Artistic License',
      # Pull in dependencies
      # See http://stackoverflow.com/questions/19738085/why-isnt-setup-py-dependency-links-doing-anything for info on use of 'dependency_links'
      # Note that pip 1.5 needs --process-dependency-links for these
      # to work; for pip 1.6 even this will be removed so then you must
      # do pip install -r requirements.txt first
      install_requires = ['genomics'],
      dependency_links=['git+https://github.com/fls-bioinformatics-core/genomics.git@devel#egg=genomics'],
      # Enable 'python setup.py test'
      test_suite='nose.collector',
      tests_require=['nose'],
      # Scripts
      scripts = scripts,
      # Configuration file for QC
      data_files = [('config',['config/settings.ini.sample']),],
      include_package_data=True,
      zip_safe = False)
