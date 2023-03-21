#!/usr/bin/env python
#
#     utils: utility classes & funcs for auto_process_ngs module
#     Copyright (C) University of Manchester 2013-2023 Peter Briggs
#
########################################################################
#
# utils.py
#
#########################################################################

"""
Utility classes and functions to support auto_process_ngs module.

Classes:

- OutputFiles:
- BufferedOutputFiles:
- ZipArchive:
- ProgressChecker:
- FileLock:
- FileLockError:
- Location: extracts information from a location specifier

Functions:

- fetch_file:
- bases_mask_is_paired_end:
- get_organism_list:
- normalise_organism_name:
- split_user_host_dir:
- get_numbered_subdir:
- find_executables:
- parse_version:
- parse_samplesheet_spec:
- pretty_print_rows:
- write_script_file:
- edit_file:
- paginate:

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import re
import logging
import zipfile
import gzip
import pydoc
import tempfile
import operator
import time
import fcntl
from urllib.request import urlopen
from urllib.error import URLError
from .command import Command
import bcftbx.utils as bcf_utils
from bcftbx.Md5sum import md5sum

# Module specific logger
logger = logging.getLogger(__name__)

######################################################################
# Magic numbers
######################################################################

MAX_OPEN_FILES = 100
DEFAULT_BUFFER_SIZE = 8192

#######################################################################
# Classes
#######################################################################

class OutputFiles:
    """Class for managing multiple output files

    Usage:

    Create a new OutputFiles instance:
    >>> fp = OutputFiles()

    Set up files against keys:
    >>> fp.open('file1','first_file.txt')
    >>> fp.open('file2','second_file.txt')

    Write content to files:
    >>> fp.write('file1','some content for first file')
    >>> fp.write('file2','content for\nsecond file')

    Append content to an existing file:
    >>> fp.open('file3','third_file.txt',append=True)
    >>> fp.write('file2','appended content')

    Check if key exists and associated file handle is
    available for writing:
    >>> 'file1' in fp
    True
    >>> 'file3' in fp
    False

    Finish and close all open files
    >>> fp.close()

    Reopen and append to a previously opened and closed
    file:
    >>> fp.open('file4','fourth_file.txt')
    >>> fp.write('file4','some content')
    >>> fp.close('file4')
    >>> fp.open('file4',append=True)
    >>> fp.write('file4','more content')

    """
    def __init__(self,base_dir=None):
        """Create a new OutputFiles instance

        Arguments:
          base_dir (str): optional 'base' directory
            which files will be created relative to

        """
        self._fp = dict()
        self._file = dict()
        self._base_dir = base_dir

    def open(self,name,filen=None,append=False):
        """Open a new output file

        'name' is the handle used to reference the
        file when using the 'write' and 'close' methods.

        'filen' is the name of the file, and is unrelated
        to the handle. If not supplied then 'name' must
        be associated with a previously closed file (which
        will be reopened).

        If 'append' is True then append to an existing
        file rather than overwriting (i.e. use mode 'at'
        instead of 'wt').

        """
        if append:
            mode = 'at'
        else:
            mode = 'wt'
        if filen is None:
            filen = self.file_name(name)
        elif self._base_dir is not None:
            filen = os.path.join(self._base_dir,filen)
        else:
            filen = os.path.abspath(filen)
        self._file[name] = filen
        self._fp[name] = open(filen,mode)

    def write(self,name,s):
        """Write content to file (newline-terminated)

        Writes 's' as a newline-terminated string to the
        file that is referenced with the handle 'name'.

        """
        self._fp[name].write("%s\n" % s)

    def file_name(self,name):
        """Get the file name associated with a handle

        NB the file name will be available even if the
        file has been closed.

        Raises KeyError if the key doesn't exist.

        """
        return self._file[name]

    def close(self,name=None):
        """Close one or all open files

        If a 'name' is specified then only the file matching
        that handle will be closed; with no arguments all
        open files will be closed.

        """
        if name is not None:
            self._fp[name].close()
            del(self._fp[name])
        else:
            names = list(self._fp.keys())
            for name in names:
                self.close(name)

    def __contains__(self,name):
        return name in self._fp

    def __len__(self):
        return len(self._fp.keys())

class BufferedOutputFiles(OutputFiles):
    """
    Class for managing multiple output files with buffering

    Version of the 'OutputFiles' class which buffers
    writing of data, to reduce number of underlying write
    operations.

    Usage is similar to OutputFiles, with additional
    'bufsize' argument which can be used to set the
    buffer size to use.
    """
    def __init__(self,base_dir=None,bufsize=DEFAULT_BUFFER_SIZE,
                 max_open_files=MAX_OPEN_FILES):
        """Create a new BufferedOutputFiles instance

        Arguments:
          base_dir (str): optional 'base' directory
            which files will be created relative to
          bufsize (int): optional buffer size; when
            data exceeds this size then it will be
            written to disk
          max_open_files (int): optional limit on the
            number of files that the instance can
            keep open internally at any one time
        """
        OutputFiles.__init__(self,base_dir=base_dir)
        self._bufsize = bufsize
        self._buffer = dict()
        self._mode = dict()
        self._max_open_files = max_open_files

    def open(self,name,filen=None,append=False):
        """Open a new output file

        'name' is the handle used to reference the
        file when using the 'write' and 'close' methods.

        'filen' is the name of the file, and is unrelated
        to the handle. If not supplied then 'name' must
        be associated with a previously closed file (which
        will be reopened).

        If the filename ends with '.gz' then the
        associated file will automatically be written as
        a gzip-compressed file.

        If 'append' is True then append to an existing
        file rather than overwriting (i.e. use mode 'at'
        instead of 'wt').

        """
        if append:
            mode = 'at'
        else:
            mode = 'wt'
        if filen is None:
            filen = self.file_name(name)
        elif self._base_dir is not None:
            filen = os.path.join(self._base_dir,filen)
        else:
            filen = os.path.abspath(filen)
        self._file[name] = filen
        self._mode[name] = mode
        if not name in self._buffer:
            self._buffer[name] = ""

    def fp(self,name):
        try:
            return self._fp[name]
        except KeyError:
            # Close a file we have too many open at once
            # (to avoid IOError [Errno 24])
            if len(self._fp) == self._max_open_files:
                # Arbitrarily close file attached to first
                # handle in the list of keys
                name0 = list(self._fp.keys())[0]
                self.close(name0)
                # Reset the mode to 'append', so the contents
                # aren't clobbered if the file is reopened
                # again later
                self._mode[name0] = 'a'
            if self._file[name].endswith('.gz'):
                open_func = gzip.open
            else:
                open_func = open
            fp = open_func(self._file[name],self._mode[name])
            self._fp[name] = fp
            return fp

    def write(self,name,s):
        """Write content to file (newline-terminated)

        Writes 's' as a newline-terminated string to the
        file that is referenced with the handle 'name'.

        """
        self._buffer[name] += "%s\n" % s
        if len(self._buffer[name]) >= self._bufsize:
            self.dump_buffer(name)

    def dump_buffer(self,name):
        self.fp(name).write(self._buffer[name])
        self._buffer[name] = ""

    def close(self,name=None):
        """Close one or all open files

        If a 'name' is specified then only the file matching
        that handle will be closed; with no arguments all
        open files will be closed.

        """
        if name is not None:
            if self._buffer[name]:
                self.dump_buffer(name)
            try:
                self._fp[name].close()
                del(self._fp[name])
            except KeyError:
                pass
        else:
            names = self._file.keys()
            for name in names:
                self.close(name)

class ZipArchive:
    """
    Utility class for creating .zip archive files

    Example usage:

    >>> z = ZipArchive('test.zip',relpath='/data')
    >>> z.add('/data/file1') # Add a single file
    >>> z.add('/data/dir2/') # Add a directory and all contents
    >>> z.close()  # to write the archive

    """
    def __init__(self,zip_file,contents=None,relpath=None,prefix=None):
        """
        Make an new zip archive instance

        Arguments:
          zip_file (str): path to the zip file to be created
          contents (list): list of file and/or directory paths
            which will be added to the zip file
          relpath (str): optional, if specified then this path
            will be stripped from the leading path for each item
            before being written (see also 'prefix')
          prefix (str): optional, if specified then this path
            will be prepended to the names of the items written
            to the archive. The prepending takes place after the
            relpath argument has been applied

        """
        self._zipfile = zipfile.ZipFile(zip_file,'w',
                                        allowZip64=True)
        self._relpath = relpath
        self._prefix = prefix
        if contents is not None:
            for item in contents:
                self.add(item)

    def add(self,item):
        """
        Add an item (file or directory) to the zip archive
        """
        item = os.path.abspath(item)
        if os.path.isfile(item):
            # Add file
            self.add_file(item)
        elif os.path.isdir(item):
            # Add directory and contents
            self.add_dir(item)
        else:
            raise Exception("ZipArchive: unknown item type for '%s'"
                            % item)

    def add_file(self,filen):
        """
        Add a file to the zip archive
        """
        if self._relpath:
            zip_pth = os.path.relpath(filen,self._relpath)
            if zip_pth.startswith('..'):
                raise Exception("Error adding a file outside %s (%s)" %
                                (self._relpath,filen))
        else:
            zip_pth = filen
        if self._prefix:
            zip_pth = os.path.join(self._prefix,zip_pth)
        self._zipfile.write(filen,zip_pth)

    def add_dir(self,dirn):
        """
        Recursively add a directory and its contents
        """
        for item in os.listdir(dirn):
            f = os.path.join(dirn,item)
            if os.path.isdir(f):
                self.add_dir(f)
            else:
                self.add_file(f)

    def close(self):
        self._zipfile.close()

    def __del__(self):
        self.close()

class ProgressChecker:
    """
    Check if an index is a multiple of a value or percentage

    Utility class to help with reporting progress of iterations
    over large numbers of items.

    Typically progress would only be reported after a certain
    number or percentage of items have been consumed; the
    ProgressChecker can be used to check if this number or
    percentage has been reached.

    Example usage: to report after every 100th item:

    >>> progress = ProgressChecker(every=100)
    >>> for i in range(10000):
    >>>    if progress.check(i):
    >>>       print("Item %d" % i)

    To report every 5% of items:

    >>> nitems = 10000
    >>> progress = ProgressChecker(percent=5,total=nitems)
    >>> for i in range(nitems):
    >>>    if progress.check(i):
    >>>       print("Item %d (%.2f%%)" % (i,progress.percent(i)))
    """
    def __init__(self,every=None,percent=None,total=None):
        """
        Create a new ProgressChecker instance

        Arguments:
          every (int): specify interval number of items
            for reporting
          percent (float): specify a percentage interval
          total (int): total number of items (must be
            provided if using `percent`)
        """
        if every is None:
            every = max(int(float(total)*float(percent)/100.0),1)
        self._every = int(every)
        self._total = total

    def check(self,i):
        """
        Check index to see if it matches the interval

        Arguments:
          i (int): index to check

        Returns:
          Boolean: True if index matches the interval,
            False if not.
        """
        return (i%self._every == 0)

    def percent(self,i):
        """
        Convert index to a percentage

        Arguments:
          i (int): index to convert

        Returns:
          Float: index expressed as a percentage of the
            total number of items.
        """
        return float(i)/float(self._total)*100.0

class FileLock:
    """
    Class for locking filesystem objects across processes

    Usage:

    >>> # Make new FileLock instance for cwd
    >>> lock = FileLock('.')
    >>> # Check that instance doesn't hold the lock
    >>> lock.has_lock
    False
    >>> # Acquire the lock
    >>> lock.acquire()
    >>> lock.has_lock
    True
    >>> # Release the lock
    >>> lock.release()

    If another FileLock instance holds the lock (within this
    process, or within another) then ``FileLock.acquire()``
    will raise an immediate ``FileLockError``; however,
    specifying a timeout period means that the instance will
    keep retrying to acquire the lock until either it is
    successful, or the timeout period is exceeded (again
    returning a ``FileLockError`` exception).

    The FileLock class can also be used as a context
    manager, e.g.

    >>> with FileLock('.'):
    ...   # Lock pwd while you do stuff

    Uses the fcntl module (see
    https://docs.python.org/3.6/library/fcntl.html,
    https://stackoverflow.com/a/32650956 and
    https://stackoverflow.com/a/55011593)

    Arguments:
      f (str): path to filesystem object (file or directory)
        to lock
      timeout (float): optional, specifies a timeout period
        after which failure to acquire the lock raises an
        exception
    """
    def __init__(self,f,timeout=None):
        # Filesystem object
        self._f = os.path.abspath(f)
        # Corresponding file descriptor
        self._lockfd = None
        # Timeout
        self._timeout = timeout

    def __enter__(self):
        """
        Called by 'with' statement to enter runtime context

        Attempts to acquire lock on specified file system
        object, and returns itself.
        """
        self.acquire()
        return self

    def __exit__(self,*args):
        """
        Called when execution leaves 'with' code block

        Releases lock on file system object that was
        acquired by ``__enter__``.
        """
        self.release()

    def acquire(self,timeout=None):
        """
        Acquire the lock on the filesystem object

        Arguments:
          timeout (float): optional, specifies a
            timeout period after which failure to
            acquire the lock raises an exception
            (NB overrides timeout set on
            instantiation).
        """
        start_time = time.time()
        if timeout is None:
            timeout = self._timeout
        if timeout is not None:
            interval = max(0.01,timeout/100.0)
        else:
            interval = None
        lockfd = os.open(self._f,os.O_RDONLY)
        while not self.has_lock:
            try:
                # Try to obtain an exclusive lock on the object
                fcntl.flock(lockfd,fcntl.LOCK_EX | fcntl.LOCK_NB)
                self._lockfd = lockfd
            except BlockingIOError:
                # The object is already locked by something else
                if timeout is not None and \
                   (time.time() - start_time) <= timeout:
                    # Timeout specified but not yet exceeded
                    # Wait a bit before trying again
                    logger.debug("FileLock: waiting for lock on %s" %
                                 self._f)
                    time.sleep(interval)
                else:
                    # No timeout, or timeout period exceeded
                    # so give up
                    break
        if not self.has_lock:
            # Failed to acquire the lock so clean up
            # and raise exception
            os.close(lockfd)
            msg = "Failed to acquire lock for %s" % self._f
            if timeout:
                msg += ": timed out after %s seconds" % timeout
            raise FileLockError(msg)

    def release(self):
        """
        Release the lock on the filesystem object
        """
        if self.has_lock:
            os.close(self._lockfd)
            self._lockfd = None

    @property
    def has_lock(self):
        """
        Check if the FileLock instance holds the lock
        """
        return (self._lockfd is not None)

class FileLockError(Exception):
    """
    Exceptions associated with the FileLock class
    """

class Location:
    """
    Class for examining a file-system location specifier

    A location specifier can be a local or a remote file or
    directory. The general form is:

    ``[[user@]server:]path``

    For a local location, only the 'path' component needs to
    be supplied.

    For a remote location, 'server' and 'path' must be
    supplied, while 'user' is optional.

    Alternatively the location can be a URL identifier of
    the form:

    ``protocol://server/path``

    The following properties are available:

    - user: the user name (or None if not specified)
    - server: the server name (or None if not specified)
    - path: the path component
    - is_remote: True if the location is on a remote host,
      False if it is local (or if it is a URL)
    - is_url: True if the location points to a URL
    - url: the URL identifier, if the location points to
      a URL (or None if not a URL)
    - protocol: the URL protocol (or None if not a URL)

    Arguments:
       location (str): location specifer of the form
          '[[user@]server:]path'
    """
    def __init__(self,location):
        self._location = location
        self._user = None
        self._server = None
        self._path = None
        self._url = None
        self._protocol = None
        for protocol in ('http','https','ftp','ftps','file',):
            protocol_ = "%s://" % protocol
            if self._location.startswith(protocol_):
                self._url = self._location
                self._protocol = protocol
                self._server = self._location[len(protocol_):].split('/')[0]
                if not self._server:
                    self._server = None
                    self._path = self._location[len(protocol_):]
                else:
                    self._path = '/'.join(self._location[len(protocol_):].split('/')[1:])
                break
        if not self._url:
            self._user,self._server,self._path = split_user_host_dir(
                self._location)
    @property
    def user(self):
        """
        Return 'user' part of '[[user@]server:]path'
        """
        return self._user
    @property
    def server(self):
        """
        Return 'server' part of '[[user@]server:]path'
        """
        return self._server
    @property
    def path(self):
        """
        Return 'path' part of '[[user@]server:]path'
        """
        return self._path
    @property
    def is_remote(self):
        """
        Check if location is on a remote server
        """
        return (self._server is not None and not self.is_url)
    @property
    def is_url(self):
        """
        Check if location is a URL
        """
        return (self._url is not None)
    @property
    def url(self):
        """
        Return path as a URL (or None if not a URL)
        """
        return self._url
    @property
    def protocol(self):
        """
        Return URL protocol (or None if not a URL)
        """
        return self._protocol
    def __repr__(self):
        return self._location


#######################################################################
# Functions
#######################################################################

def fetch_file(src,dest=None):
    """
    Fetch a copy of a file from an arbitrary location

    Gets a copy of a file which can be specified as either
    a local or a remote path (i.e. using the syntax
    ``[[USER@]HOST:]PATH``) or as URL.

    If a destination file name is not supplied then the
    destination file name will be that of the source file
    and will be copied to the current working directory;
    if the supplied destination is a directory then the
    destination file name will be that of the source file
    and will be copied to that directory.

    The destination must be a local file or directory.

    Arguments:
      src (str): path or URL of the file to be copied
      dest (str): optional, local destination file name
        or directory to copy ``src`` to

    Returns:
      String: path to the copy.
    """
    logger.debug("Fetching file '%s'" % src)
    src = Location(src)
    # Sort out destination
    if not dest:
        # Use the source file name
        dest = os.path.basename(src.path)
    elif os.path.isdir(dest):
        # Specified destination is a directory
        dest = os.path.join(dest,
                             os.path.basename(src.path))
    dest = os.path.abspath(dest)
    # Copy the file
    if src.is_url:
        # File is a URL
        with tempfile.TemporaryFile('w+b') as fp:
            # Get URL contents
            try:
                urlfp = urlopen(src.url)
                fp.write(urlfp.read())
            except URLError as ex:
                # Failed to download from URL
                raise Exception("Error fetching file from '%s': %s" %
                                (src.url,ex))
            # Copy content to local file
            fp.seek(0)
            with open(dest,'wb') as fpp:
                fpp.write(fp.read())
    else:
        # File is on a local or remote server
        if src.is_remote:
            src = str(src)
            cp = 'scp'
        else:
            src = os.path.abspath(src.path)
            cp = 'cp'
        copy_command = Command(cp,src,dest)
        logger.debug("Running '%s'" % copy_command)
        #status = copy_command.run_subprocess()
        status,output = copy_command.subprocess_check_output()
        if status != 0:
            raise Exception("Error fetching file from '%s': %s" %
                            (src,output))
    return dest

def bases_mask_is_paired_end(bases_mask):
    # Determine if run is paired end based on bases mask string
    non_index_reads = []
    for read in bases_mask.split(','):
        try:
            read.index('I')
        except ValueError:
            non_index_reads.append(read)
    if len(non_index_reads) == 2:
        # Paired end
        return True
    elif len(non_index_reads) < 2:
        # Single end
        return False
    else:
        # An error?
        raise Exception("Bad bases mask '%s'?" % bases_mask)

def get_organism_list(organisms):
    """
    Return a list of normalised organism names

    Normalisation consists of converting names to
    lower case and spaces to underscores.

    E.g.

    "Human,Mouse" -> ['human','mouse']
    "Xenopus tropicalis" -> ['xenopus_tropicalis']

    Arguments:
      organisms (str): string with organism names
        separated by commas

    Returns:
      List: list of normalised organism names
    """
    if not organisms:
        return []
    return [normalise_organism_name(x)
            for x in str(organisms).split(',')]

def normalise_organism_name(name):
    """
    Return normalised organism name

    Normalisation consists of converting names to
    lower case and spaces to underscores.

    E.g.

    "Human" -> 'human'
    "Xenopus tropicalis" -> 'xenopus_tropicalis'

    Arguments:
      name (str): organism name

    Returns:
      String: normalised organism name
    """
    # Make lowercase, split on commas and replace whitespace
    # with single underscore, then strip leading/trailing
    # underscores
    return re.sub(r'\s+','_',str(name).lower()).strip('_')

def split_user_host_dir(location):
    # Split a location of the form [[user@]host:]dir into its
    # user, hostname and directory components
    try:
        location = location.strip()
    except AttributeError:
        # Not a string?
        logger.error("Bad input to split_user_host_dir: '%s'" % location)
        return (None,None,None)
    if not location:
        return (None,None,None)
    try:
        location.index(':')
        location,dirn = location.split(':')
        try:
            location.index('@')
            user,host = location.split('@')
        except ValueError:
            user = None
            host = location
    except ValueError:
        user = None
        host = None
        dirn = location
    return (user,host,dirn)

def get_numbered_subdir(name,parent_dir=None,full_path=False):
    """
    Return a name for a new numbered log subdirectory

    Generates the name for a numbered subdirectory.

    Subdirectories are named as NNN_<name>  e.g.
    001_setup, 002_make_fastqs etc.

    'Gaps' are ignored, so the number associated with
    the new name will be one plus the highest index
    that already exists.

    **Note that a directory is not created** - this
    must be done by the calling subprogram. As a
    result there is the possibility of a race
    condition.

    Arguments:
      name (str): name for the subdirectory
        (typically the name of the processing
        stage that will produce logs to be
        written to the subdirs
      parent_dir (str): path to the parent
        directory where the indexed directory
        would be created; defaults to CWD if
        not set
      full_path (bool): if True then return the
        full path for the new subdirectory;
        default is to return the name relative
        to the parent directory

    Returns:
      String: name for the new log subdirectory
        (will be the full path if 'full_path'
        was specified).
    """
    # Sort out parent directory
    if parent_dir is None:
        parent_dir = os.getcwd()
    parent_dir = os.path.abspath(parent_dir)
    # Get the highest number from the names of
    # any other existing numbered subdirs
    i = 0
    for d in bcf_utils.list_dirs(parent_dir):
        try:
            i = max(i,int(d.split('_')[0]))
        except ValueError:
            pass
    # Generate and return name/path
    subdir = "%03d_%s" % (i+1,str(name))
    if full_path:
        subdir = os.path.join(parent_dir,subdir)
    return subdir

def find_executables(names,info_func,reqs=None,paths=None):
    """
    List available executables matching list of names

    By default searches the PATH for the executables listed
    in 'names', using the supplied 'info_func' to acquire
    package names and versions of each, returns a list of
    executables with the full path, package and version.

    'info_func' is a function that must be supplied by the
    calling subprogram. Its signature should look like:

    >>> def info_func(p):
    ...   # Determine full_path, package_name and
    ...   # version
    ...   # Then return these as a tuple
    ...   return (full_path,package_name,version)

    The 'reqs' argument allows a specific version or range
    of versions to be requested; in this case the returned
    list will only contain those packages which satisfy
    the requested versions.

    A range of version specifications can be requested by
    separating multiple specifiers with a comma - for
    example '>1.8.3,<2.16'.

    The full set of operators is:

    - ==, >, >=, <=, <

    If no versions are requested then the packages will
    be returned in PATH order; otherwise they will be
    returned in version order (highest to lowest).

    Arguments:
      names (list): list of executable names to look for.
        These can be full paths or executables with no
        leading paths
      info_func (function): function to use to get tuples
        of (full_path,package_name,version) for an
        executable
      reqs (str): optional version requirement expression
        (for example '>=1.8.4'). If supplied then only
        executables fulfilling the requirement will be
        returned. If no operator is supplied then '=='
        is implied.
      paths (list): optional set of directory paths to
        search when looking for executables. If not
        supplied then the set of paths specified in
        the PATH environment variable will be searched.

    Returns:
      List: full paths to executables matching specified
        criteria.

    """
    # Search paths
    if paths is None:
        paths = os.environ['PATH'].split(os.pathsep)
    # Search for executables
    available_exes = []
    for path in paths:
        if os.path.isfile(path):
            path = os.path.dirname(path)
        for name in names:
            prog_path = os.path.abspath(os.path.join(path,name))
            if bcf_utils.PathInfo(prog_path).is_executable:
                available_exes.append(prog_path)
    # Filter on requirement
    if reqs:
        # Loop over ranges
        for req in reqs.split(','):
            logger.debug("Filtering on expression: %s" % req)
            # Determine operator and version
            req_op = None
            req_version = None
            for op in ('==','>=','<=','>','<'):
                if req.startswith(op):
                    req_op = op
                    req_version = req[len(op):].strip()
                    break
            if req_version is None:
                req_op = '=='
                req_version = req.strip()
            logger.debug("Required version: %s %s" % (req_op,req_version))
            if req_op == '==':
                op = operator.eq
            elif req_op == '>=':
                op = operator.ge
            elif req_op == '>':
                op = operator.gt
            elif req_op == '<':
                op = operator.lt
            elif req_op == '<=':
                op = operator.le
            # Filter the available executables on version
            logger.debug("Pre filter: %s" % available_exes)
            logger.debug("Versions  : %s" % [info_func(x)[2]
                                              for x in available_exes])
            available_exes = list(
                filter(lambda x: op(
                    parse_version(info_func(x)[2]),
                    parse_version(req_version)),
                       available_exes))
            logger.debug("Post filter: %s" % available_exes)
        # Sort into version order, highest to lowest
        available_exes.sort(
            key=lambda x: parse_version(info_func(x)[2]),
            reverse=True)
        logger.debug("Post sort: %s" % available_exes)
        logger.debug("Versions : %s" % [info_func(x)[2]
                                 for x in available_exes])
    if available_exes:
        print("%s: found available packages:" % ','.join(names))
        for i,package in enumerate(available_exes):
            package_info = info_func(package)
            print("%s %s\t%s\t%s" % (('*' if i == 0 else ' '),
                                     package_info[0],
                                     package_info[1],
                                     package_info[2]))
    else:
        print("%s: No packages found" % ','.join(names))
    return available_exes

def parse_version(s):
    """
    Split a version string into a tuple for comparison

    Given a version string of the form e.g. "X.Y.Z",
    return a tuple of the components e.g. (X,Y,Z)

    Where possible components will be coverted to
    integers.

    If the version string is empty then the version
    number will be set to an arbitrary negative
    integer.

    Typically the result from this function would not
    be used directly, instead it is used to compare
    two versions, for example:

    >>> parse_version("2.17") < parse_version("1.8")
    False

    Arguments:
      s (str): version string

    Returns:
      Tuple: tuple of the version string
    """
    if s == "":
        # Essentially versionless; set to an
        # arbitrarily small integer
        s = "-99999"
    items = []
    for i in s.split('.'):
        try:
            i = int(i)
        except ValueError:
            pass
        items.append(i)
    return tuple(items)

def parse_samplesheet_spec(s):
    """
    Split sample sheet line specification into components

    Given a sample sheet line specification of the form
    "[LANES:][COL=PATTERN:]VALUE", split into the
    component parts and return as a tuple:

    (VALUE,LANES,COL,PATTERN)

    where 'LANES' is a list of lanes (or 'None', if not
    present), and the others will be strings (or 'None',
    if not present).

    Arguments:
      s (str): specification string

    Returns:
      Tuple: the extracted components as
        (VALUE,LANES,COL,PATTERN).
    """
    # Initialise
    value = None
    lanes = None
    column = None
    pattern = None
    # Split input string
    items = s.split(':')
    # Value is at the end
    value = items[-1]
    items = items[:-1]
    # Optional column matching
    if items:
        if '=' in items[-1]:
            # Pattern to match
            column,pattern = items[-1].split('=')
            items = items[:-1]
    # Optional lanes
    if items:
        # Lane specification
        lanes = bcf_utils.parse_lanes(items[0])
    # Return as tuple
    return (value,lanes,column,pattern)

def pretty_print_rows(data,prepend=False):
    """Format row-wise data into 'pretty' lines

    Given 'row-wise' data (in the form of a list of lists),
    for example:

    [['hello','A salutation'],[goodbye','The End']]

    formats into a string of text where lines are
    newline-separated and the 'fields' are padded with
    spaces so that they line up left-justified in
    columns, for example:

    hello   A salutation
    goodbye The End

    Arguments:
      data: row-wise data as a list of lists
      prepend: (optional), if True then columns
        are right-justified (i.e. padding is
        added before each value).

    """
    # Get maximum field widths for each column
    widths = []
    for row in data:
        for i in range(len(row)):
            width = len(str(row[i]))
            try:
                widths[i] = max(width,widths[i])
            except IndexError:
                widths.append(width)
    # Build output
    output = []
    for row in data:
        line = []
        for item,width in zip([str(x) for x in row],widths):
            padding = ' '*(width-len(item))
            if prepend:
                line.append(padding + item)
            else:
                line.append(item + padding)
        output.append(' '.join(line))
    return '\n'.join(output)

def write_script_file(script_file,contents,append=False,shell=None):
    """Write command to file

    Arguments:
      script_file (str): path of file to write command to
      contents (str): content to write to the file
      append (bool): optional, if True and script_file exists
        then append content (default is to overwrite existing
        contents) 
      shell: optional, if set then defines the shell to
        specify after '!#'

    """
    if append:
        mode = 'at'
    else:
        mode = 'wt'
    with open(script_file,mode=mode) as fp:
        if (not append) and (shell is not None):
            fp.write("#!%s\n" % shell)
        fp.write("%s\n" % contents)
    os.chmod(script_file,0o775)

def edit_file(filen,editor="vi",append=None):
    """
    Send a file to an editor

    Creates a temporary copy of a file and opens an
    editor to allow the user to make changes. Any
    edits are saved back to the original file.

    Arguments:
      filen (str): path to the file to be edited
      editor (str): optional, editor command to be used
        (will be overriden by user's EDITOR environment
        variable even if set). Defaults to 'vi'.
      append (str): optional, if set then append the
        supplied text to the end of the file before
        editing. NB the text will only be kept if the
        user saves a change to the file in the editor.

    """
    # Acquire an editor command
    try:
        editor = os.environ["EDITOR"]
    except KeyError:
        pass
    if editor is None:
        logger.critical("No editor specified!")
        return
    # Make a temporary copy for editing
    f,tmpfile = tempfile.mkstemp()
    os.fdopen(f).close()
    with open(tmpfile,'wt') as fp:
        if os.path.exists(filen):
            fp.write(open(filen,'r').read())
        else:
            fp.write()
        if append:
            fp.write("%s\n" % str(append))
    checksum = md5sum(tmpfile)
    # Build command line to run the editor
    editor = str(editor).split(' ')
    edit_cmd = Command(editor[0],*editor[1:])
    edit_cmd.add_args(tmpfile)
    edit_cmd.run_subprocess()
    # Finished
    if md5sum(tmpfile) != checksum:
        with open(filen,'wt') as fp:
            fp.write(open(tmpfile,'r').read())
            os.remove(tmpfile)
    else:
        logger.warning("no changes to write")

def paginate(text):
    """
    Send text to stdout with pagination

    If the function detects that the stdout is an interactive
    terminal then the supplied text will be piped via a
    paginator command.

    The pager command will be the default for ``pydoc``, but
    can be over-ridden by the ``PAGER`` environment variable.

    If stdout is not a terminal (for example if it's being
    set to a file, or piped to another command) then the
    pagination is skipped.

    Arguments:
      text (str): text to be printed using pagination

    """
    # If stdout is a terminal
    if os.isatty(sys.stdout.fileno()):
        # Acquire a pager command
        try:
            pager = os.environ["PAGER"]
        except KeyError:
            pager = None
        # Output the prediction with paging
        if pager is not None:
            pydoc.pipepager(text,cmd=pager)
        else:
            pydoc.pager(text)
    else:
        # Stdout not a terminal
        print(text)
