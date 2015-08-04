#!/bin/env python
#
#     envmod.py: wrapper for setting up environment with modulefiles
#     Copyright (C) University of Manchester 2014-2015 Peter Briggs
#
########################################################################
#
# envmod.py
#
#########################################################################

__version__ = "2.0.0"
__patched__ = False

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

Note for Ubuntu 14.04
---------------------

The Python implementation of the environment modules package
appears to be broken on Ubuntu 14.04. If this detected then the
module will issue a debug message and attempt to patch it with
it's own version of the 'module' function.

To find out if a patched has been applied, check the __patched__
property of the module e.g.

>>> from auto_process_ngs import envmod
>>> envmod.__patched__
True

"""

#######################################################################
# Setup
#######################################################################

import os
import logging
modules_python_dirs = ('/usr/share/Modules/init/',
                       '/usr/share/modules/init/',)
for d in modules_python_dirs:
    modules_python = os.path.join(d,'python.py')
    if not os.path.isfile(modules_python):
        modules_python = None
    else:
        break
if modules_python is None:
    # Nothing to load
    raise ImportError("No 'python.py' file in any of %s" % modules_python_dirs)

# Load the code from the env modules python file
try:
    execfile(modules_python)
except Exception,ex:
    raise ImportError("Exception executing code from %s: %s: %s" % 
                      (modules_python,ex.__class__.__name__,ex))

# Check for broken 'module' function from python.py (e.g. Ubuntu 14.04)
try:
    module('list')
except TypeError:
    # Patch the module function
    logging.debug("Patching broken 'module' function from %s" % modules_python)
    def module(*args):
	if type(args[0]) == type([]):
            args = args[0]
	else:
            args = list(args)
        (output,error) = subprocess.Popen(['/usr/bin/modulecmd','python']+args, stdout=subprocess.PIPE).communicate()
        exec output
    # Indicate that the module is patched
    __patched__ = True

#######################################################################
# Functions
#######################################################################

def loaded():
    """Return a list of loaded environment modules

    Loaded modulefiles are extracted from the
    'LOADEDMODULES' environment variable.

    """
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
            print "Loading environment module: %s" % arg
            module('load',arg)
            if arg not in loaded():
                raise Exception("Failed to load environment module '%s'" % arg)

def unload(*args):
    """Unload environment modules

    """
    for arg in args:
        if arg:
            print "Unloading environment module: %s" % arg
            module('unload',arg)

def unload_all():
    """Unload all environment modules that are currently loaded

    """
    unload(*loaded())

