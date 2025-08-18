#!/usr/bin/env python
#
#     decorators.py: decorator functions for AutoProcessor
#     Copyright (C) University of Manchester 2023 Peter Briggs
#

"""
Implements decorators for use with the ``AutoProcessor`` class.
"""

#######################################################################
# Imports
#######################################################################

import time
import logging
from functools import wraps

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Functions
#######################################################################

def add_command(name,f):
    """
    Add a method to an AutoProcessor class

    Implements an '@add_command' decorator which can be
    used to add a function to an AutoProcessor class as
    a new method (aka 'command').

    When the method is invoked additional output is added
    to report the run ID and the command name, and to mark
    the start and end of the command execution.

    Additionally the command execution is wrapped in a
    ``try/except`` block.

    For example:

    >>> def hello(cls):
    ...    print("Hello %s" % cls.person)
    ...
    >>> @command("greeting",hello)
    ... class Example:
    ...   def __init__(self):
    ...      self.person = "World"
    ...   def __str__(self):
    ...      return "Example: '%s'" % self.person
    ...
    >>> Example().greeting()
    [2023-01-11 09:36:04] Example: 'World'
    [2023-01-11 09:36:04] Running 'greeting' command
    Hello World
    [2023-01-11 09:36:04] greeting: finished

    The function must accept a class instance as the
    first argument.

    Arguments:
      name (str): name of the command (which will be the
        method name when added)
      f (object): callable object (e.g. function) that
        will be invoked by the command. The first
        argument of the callable must be an
        AutoProcessor-like class
    """
    @wraps(f)
    def wrapped_func(*args,**kws):
        # Wraps execution of the supplied
        # function to trap exceptions and
        # add additional commentary
        try:
            cls = args[0]
            print("[%s] %s" % (timestamp(),str(cls)))
        except Exception as ex:
            pass
        print("[%s] Running '%s' command" % (timestamp(),name))
        try:
            ret = f(*args,**kws)
        except Exception as ex:
            logger.fatal("%s: %s" % (name,ex))
            ret = 1
        else:
            print("[%s] %s: finished" % (timestamp(),name))
        if ret is None:
            # Assume success
            ret = 0
        return ret
    def timestamp():
        # Return the current time
        return time.strftime("%Y-%m-%d %H:%M:%S")
    def wrapper(cls):
        # Adds the supplied function to
        # to the class
        setattr(cls,name,wrapped_func)
        return cls
    return wrapper
