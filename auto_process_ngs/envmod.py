#!/bin/env python
#
#     envmod.py: wrapper for setting up environment with modulefiles
#     Copyright (C) University of Manchester 2014 Peter Briggs
#
########################################################################
#
# envmod.py
#
#########################################################################

__version__ = "1.0.1"

"""envmod

Wrapper to set up the environment using modulefiles:

http://modules.sourceforge.net/

This is based on code from the Stack Overflow Q&A at:

http://modules.sourceforge.net/

It exposes the 'module' command provided by the environment module
package, but also provides some wrapper commands which could be
used preferentially:

- load(module[,module,...])
- unload(module[,module...])

Example usage:

>>> import envmod
>>> envmod.load('apps/fastqc/0.11.2')

"""

#######################################################################
# Setup
#######################################################################

import os
modules_python = '/usr/share/Modules/init/python.py'
if not os.path.isfile(modules_python):
    # Nothing to load
    raise ImportError("No file %s" % modules_python)

# Load the code from the env modules python file
try:
    execfile(modules_python)
except Exception,ex:
    raise ImportError("Exception executing code from %s: %s: %s" % 
                      (modules_python,ex.__class__.__name__,ex))

#######################################################################
# Functions
#######################################################################

def load(*args):
    """Load one or more environment modules

    """
    for arg in args:
        if arg:
            print "Loading environment module: %s" % arg
            module('load',arg)

def unload(*args):
    """Unload environment modules

    """
    for arg in args:
        if arg:
            print "Unloading environment module: %s" % arg
            module('unload',arg)

