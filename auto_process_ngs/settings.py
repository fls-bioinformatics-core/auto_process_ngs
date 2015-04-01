# Settings for auto_process module
#

import os
import sys
import logging
import bcftbx.JobRunner as JobRunner
from bcftbx.utils import AttributeDictionary
from config import Config

# Locate settings files
# Search in 'config' subdir of installation location, then
# installation location, then in current directory
__install_path = os.path.abspath(os.path.normpath(
    os.path.join(os.path.dirname(sys.argv[0]),'..')))
__config_file_path = (os.path.join(__install_path,'config'),
                      __install_path,
                      os.getcwd(),)
__config_file = None
__sample_config_file = None
for path in __config_file_path:
    __config_file = os.path.join(path,'settings.ini')
    if os.path.exists(__config_file):
        # Located settings file
        break
    # No settings file here, look for a sample version
    if __sample_config_file is None:
        __sample_config_file = __config_file + '.sample'
        if not os.path.exists(__sample_config_file):
            __sample_config_file = None
    # Reset config file to keep looking
    __config_file = None

# No settings.ini file - try to make one
if __config_file is None:
    logging.warning("No local settings file found in %s" % ', '.join(__config_file_path))
    if __sample_config_file is not None:
        logging.warning("Attempting to make a copy from sample settings file")
        __config_file = os.path.splitext(__sample_config_file)[0]
        try:
            open(__config_file,'w').write(open(__sample_config_file,'r').read())
            logging.warning("Created new file %s" % __config_file)
            logging.warning("Edit configuration settings and rerun")
        except Exception,ex:
            raise Exception("Failed to create %s: %s" % (__config_file,ex))
    else:
        raise Exception("No sample config file found")
    sys.exit(1)
#
# Import site-specific settings from local version
config = Config()
config.read(__config_file)
#
# General parameters
general = AttributeDictionary()
general['default_runner'] = config.get('general','default_runner','SimpleJobRunner')
general['max_concurrent_jobs'] = config.getint('general','max_concurrent_jobs',12)
#
# modulefiles
modulefiles = AttributeDictionary()
modulefiles['make_fastqs'] = config.get('modulefiles','make_fastqs')
modulefiles['run_qc'] = config.get('modulefiles','run_qc')
#
# bcl2fastq
bcl2fastq = AttributeDictionary()
bcl2fastq['nprocessors'] = config.getint('bcl2fastq','nprocessors',1)
#
# fastq_stats
fastq_stats = AttributeDictionary()
fastq_stats['nprocessors'] = config.getint('fastq_stats','nprocessors',1)
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
#
# Report configuration settings
def report_settings():
    """
    Report the settings read from the config file
    """
    print "Settings from %s:" % __config_file
    show_dictionary('general',general)
    show_dictionary('modulefiles',modulefiles)
    show_dictionary('bcl2fastq',bcl2fastq)
    show_dictionary('runners',runners)
    show_dictionary('archive',archive)
    show_dictionary('qc_web_server',qc_web_server)

