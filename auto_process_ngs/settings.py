# Settings for auto_process module
#
version = "0.0.81"

import os
import sys
import logging
import JobRunner
from utils import AutoProcessConfigParser
from bcf_utils import AttributeDictionary

# Locate settings file
__install_path = os.path.abspath(os.path.normpath(
    os.path.join(os.path.dirname(sys.argv[0]),'..')))
__config_file = os.path.join(__install_path,'settings.ini')
if not os.path.exists(__config_file):
    # No settings.ini file - try to make one
    logging.warning("No local settings file in %s" % __install_path)
    logging.warning("Attempting to make a copy from sample settings file")
    __sample_config_file = os.path.join(__install_path,'settings.ini.sample')
    if not os.path.exists(__sample_config_file):
        raise Exception,"No sample config file: %s" % __sample_config_file
    try:
        open(__config_file,'w').write(open(__sample_config_file,'r').read())
    except Exception,ex:
        raise Exception,"Failed to create settings file: %s" % ex
    logging.warning("Created new file %s" % __config_file)
    logging.warning("Edit configuration settings and rerun")
    sys.exit(1)
#
# Import site-specific settings from local version
config = AutoProcessConfigParser()
config.read(__config_file)
#
# General parameters
general = AttributeDictionary()
general['default_runner'] = config.get('general','default_runner','SimpleJobRunner')
general['max_concurrent_jobs'] = config.getint('general','max_concurrent_jobs',12)
#
# bcl2fastq
bcl2fastq = AttributeDictionary()
bcl2fastq['nprocessors'] = config.getint('bcl2fastq','nprocessors',1)
#
# Define runners for specific jobs
runners = AttributeDictionary()
for name in ('bcl2fastq','qc','stats',):
    runners[name] = config.getrunner('runners',name,general.default_runner)
#
# Information for archiving analyses
# dirn should be a directory in the form [[user@]host:]path]
archive = AttributeDictionary()
archive['dirn'] = config.get('archive','dirn',None)
archive['log'] = config.get('archive','log',None)
archive['group'] = config.get('archive','group',None)
archive['chmod'] = config.get('archive','chmod',None)
#
# Information for uploading QC reports
# dirn should be a directory in the form [[user@]host:]path]
qc_web_server = AttributeDictionary()
qc_web_server['dirn'] = config.get('qc_web_server','dirn',None)
qc_web_server['url'] = config.get('qc_web_server','url',None)
#
# Show settings
def show_dictionary(name,d):
    """Print the contents of a dictionary

    """
    print "[%s]" % name
    for key in d:
        print "\t%s = %s" % (key,d[key])
##show_dictionary('general',general)
##show_dictionary('bcl2fastq',bcl2fastq)
##show_dictionary('runners',runners)
##show_dictionary('archive',archive)
##show_dictionary('qc_web_server',qc_web_server)

