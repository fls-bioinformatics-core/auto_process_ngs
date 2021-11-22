#!/usr/bin/env python
#
#     download_fastqs.py: utility to download fastq files from web server
#     Copyright (C) University of Manchester 2015,2021 Peter Briggs
#
#########################################################################
#
# download_fastqs.py
#
#########################################################################

"""
Utility to download Fastq files from web server.

"""

#######################################################################
# Imports
#######################################################################

import argparse
import tempfile
import re
import sys
import os
import io
try:
    import hashlib
except ImportError:
    # hashlib not available, use deprecated md5 module
    import md5
from urllib.request import urlopen

#######################################################################
# Constants
#######################################################################

BLOCKSIZE = 1024*1024

#######################################################################
# Functions
#######################################################################

def print_(s,end='\n'):
    # Implement a local version of 'print' function which
    # mimics the 'end' argument from Python3
    sys.stdout.write(u"%s%s" % (s,end))

def download_file(url,dest):
    # See http://stackoverflow.com/a/22776/579925
    u = urlopen(url)
    f = io.open(dest,'wb')
    meta = u.info()
    file_size = int(meta.get("Content-Length"))
    print_("Downloading: %s Bytes: %s " % (dest, file_size))
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"% 10d  [%3.2f%%]" % (file_size_dl,
                                        file_size_dl*100.0/file_size)
        status = status + chr(8)*(len(status)+1)
        print_(status,end='')
    f.close()
    print_("")

def check_md5sum(filen,md5):
    try:
        c = hashlib.md5()
    except NameError:
        c = md5.new()
    f = io.open(filen,'rb')
    for block in iter(lambda: f.read(BLOCKSIZE),b''):
        c.update(block)
    c = c.digest()
    try:
        # Py3 method
        chksum = ("%02x"*len(c)) % tuple(map(int,c))
    except ValueError:
        # Fallback to Py2 method
        chksum = ("%02x"*len(c)) % tuple(map(ord,c))
    return (chksum == md5)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Download checksum file and fastqs from "
        "URL into current directory (or directory DIR, if "
        "specified), and verify the downloaded files against "
        "the checksum file.")
    p.add_argument('url',metavar="URL",
                   help="URL with checksum file and fastqs")
    p.add_argument('dir',metavar="DIR",
                   nargs='?',help="directory to put downloaded "
                   "fastqs into (defaults to current directory)")
    args = p.parse_args()
    if args.dir:
        dest_dir = os.path.abspath(args.dir)
    else:
        dest_dir = os.getcwd()
    print_("Moving to %s" % dest_dir)
    try:
        os.chdir(dest_dir)
    except Exception as ex:
        if not os.path.isdir(dest_dir):
            sys.stderr.write("ERROR directory doesn't exist?\n")
        else:
            sys.stderr.write("ERROR %s\n" % ex)
        sys.exit(1)
    print_("Fetching index from %s" % args.url)
    index_page = str(urlopen(args.url).read())
    print_("Locating chksum file",end='')
    for line in index_page.split('\n'):
        chksum_file = re.search('[A-Za-z0-9_\-]*\.chksums',line)
        if chksum_file:
            break
    if not chksum_file:
        sys.stderr.write("ERROR failed to find a chksum file\n")
        sys.exit(1)
    print_(chksum_file.group(0))
    chksum_file = chksum_file.group(0)
    # Download chksum file
    print_("Fetching %s" % chksum_file)
    download_file(os.path.join(args.url,chksum_file),chksum_file)
    # Get file list and checksums
    chksums = dict()
    with io.open(chksum_file,'rt') as fp:
        for line in fp:
            chksum,filen = line.strip('\n').split()
            chksums[filen] = chksum
    file_list = sorted(list(chksums.keys()))
    # Flag for checking MD5 sums
    checksum_errors = False
    # Download the files
    print_("Downloading %d files" % len(file_list))
    for f in file_list:
        if os.path.exists(f):
            print_("%s: file exists, skipping download" % f)
        else:
            download_file(os.path.join(args.url,f),f)
        if check_md5sum(f,chksums[f]):
            print_("%s: checksum OK" % f)
        else:
            print_("%s: checksum FAILED" % f)
            checksum_errors = True
    if checksum_errors:
        sys.stderr.write("ERROR one or more checksums failed, see above\n")
        sys.exit(1)
    else:
        print_("All checksums verified ok")

