#!/usr/bin/env python
#
#     envmod.py: wrapper for setting up environment with modulefiles
#     Copyright (C) University of Manchester 2014-2015 Peter Briggs
#
########################################################################
#
# envmod.py
#
#########################################################################

__version__ = "2.1.0"

"""envmod

Wrapper to set up the environment using modulefiles:

http://modules.sourceforge.net/

This is based on code from the Stack Overflow Q&A at:

http://modules.sourceforge.net/

It exposes the 'module' command provided by the environment module
package, but also provides some wrapper commands which could be
used preferentially:

- load(module[,module,...])  : loads specified modulefile(s)
- unload(module[,module...]) : unloads modulefile(s)
- unload_all()               : unloads all currently loaded modulefiles
- loaded()                   : lists loaded modulefiles

Example usage:

>>> import envmod
>>> envmod.load('apps/fastqc/0.11.2')

Attempting to load a non-existent environment module will raise
an exception.

Note that environment modules already loaded in the shell from
where Python was invoked will also be available; these modules can
be unloaded within the Python environment.

"""

#######################################################################
# Setup
#######################################################################

__ENVMODULES__ = False

import os
import logging
modules_python_dirs = ('/usr/share/Modules/init/',
                       '/usr/share/modules/init/',
                       '/opt/clusterware/opt/Modules/init/',
                       '/opt/clusterware/opt/modules/init/')
for d in modules_python_dirs:
    modules_python = os.path.join(d,'python.py')
    if not os.path.isfile(modules_python):
        modules_python = None
    else:
        break
modulecmd_dirs = ('/usr/bin',
                  '/opt/clusterware/opt/Modules/bin',
                  '/opt/clusterware/opt/modules/bin',)
for d in modulecmd_dirs:
    modulecmd = os.path.join(d,'modulecmd')
    if not os.path.isfile(modulecmd):
        modulecmd = None
    else:
        break
if modules_python is None or modulecmd is None:
    # Nothing to load
    logging.debug("No 'python.py' file in any of %s" % str(modules_python_dirs))
else:
    # Load the code from the env modules python file
    try:
        execfile(modules_python)
        __ENVMODULES__ = True
    except Exception as ex:
        logging.debug("Exception executing code from %s: %s: %s" % 
                      (modules_python,ex.__class__.__name__,ex))

# Function to check if env modules are available
def _got_envmodules():
    return __ENVMODULES__

# Add our own version of the 'module' function to deal with
# the python implementation being broken on some systems (e.g.
# Ubuntu 14.04)
def module(*args):
    if not _got_envmodules():
        return
    if type(args[0]) == type([]):
        args = args[0]
    else:
        args = list(args)
    (output,error) = subprocess.Popen([modulecmd,'python']+args, stdout=subprocess.PIPE).communicate()
    exec output

#######################################################################
# Functions
#######################################################################

def loaded():
    """Return a list of loaded environment modules

    Loaded modulefiles are extracted from the
    'LOADEDMODULES' environment variable.

    """
    if not _got_envmodules():
        return []
    try:
        loadedmodules = os.environ['LOADEDMODULES']
        if not loadedmodules:
            return []
        return loadedmodules.split(':')
    except KeyError:
        # No loaded modulefiles?
        return []

def load(*args):
    """Load one or more environment modules

    Raises an exception if the module cannot be loaded.

    """
    for arg in args:
        if arg:
            print("Loading environment module: %s" % arg)
            module('load',arg)
            if arg not in loaded():
                raise Exception("Failed to load environment module '%s'" % arg)

def unload(*args):
    """Unload environment modules

    """
    for arg in args:
        if arg:
            print("Unloading environment module: %s" % arg)
            module('unload',arg)

def unload_all():
    """Unload all environment modules that are currently loaded

    """
    unload(*loaded())

