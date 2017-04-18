#!/usr/bin/env python
#
#     download_fastqs.py: utility to download fastq files from web server
#     Copyright (C) University of Manchester 2015 Peter Briggs
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

import optparse
import tempfile
from urllib2 import urlopen
import re
import sys
import os
try:
    import hashlib
except ImportError:
    # hashlib not available, use deprecated md5 module
    import md5

#######################################################################
# Constants
#######################################################################

BLOCKSIZE = 1024*1024

#######################################################################
# Functions
#######################################################################

def download_file(url,dest):
    # See http://stackoverflow.com/a/22776/579925
    u = urlopen(url)
    f = open(dest,'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print "Downloading: %s Bytes: %s" % (dest, file_size)
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8)*(len(status)+1)
        print status,
    f.close()

def check_md5sum(filen,md5):
    try:
        c = hashlib.md5()
    except NameError:
        c = md5.new()
    f = open(filen,'rb')
    for block in iter(lambda: f.read(BLOCKSIZE),''):
        c.update(block)
    c = c.digest()
    chksum = ("%02x"*len(c)) % tuple(map(ord,c))
    return (chksum == md5)

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = optparse.OptionParser(usage="%prog [OPTIONS] URL [DIR]",
                              description="Download checksum file and fastqs from "
                              "URL into current directory (or directory DIR, if "
                              "specified), and verify the downloaded files against "
                              "the checksum file.")
    options,args = p.parse_args()
    if not len(args):
        p.error("Supply a URL to download Fastqs from, and an option target d")
    elif len(args) > 2:
        p.error("Too many arguments, use -h for usage information")
    url = args[0]
    if len(args) == 2:
        dest_dir = os.path.abspath(args[1])
        print "Moving to %s" % dest_dir
        try:
            os.chdir(dest_dir)
        except Exception,ex:
            if not os.path.isdir(dest_dir):
                sys.stderr.write("ERROR directory doesn't exist?\n")
            else:
                sys.stderr.write("ERROR %s\n" % ex)
            sys.exit(1)
    print "Fetching index from %s" % url
    index_page = urlopen(url).read()
    print "Locating chksum file",
    for line in index_page.split('\n'):
        chksum_file = re.search('[A-Za-z0-9_]*\.chksums',line)
        if chksum_file:
            break
    if not chksum_file:
        sys.stderr.write("ERROR failed to find a chksum file\n")
        sys.exit(1)
    print chksum_file.group(0)
    chksum_file = chksum_file.group(0)
    # Download chksum file
    print "Fetching %s" % chksum_file
    download_file(os.path.join(url,chksum_file),chksum_file)
    # Get file list and checksums
    chksums = dict()
    with open(chksum_file,'r') as fp:
        for line in fp:
            chksum,filen = line.strip('\n').split()
            chksums[filen] = chksum
    file_list = chksums.keys()
    file_list.sort()
    # Flag for checking MD5 sums
    checksum_errors = False
    # Download the files
    print "Downloading %d files" % len(file_list)
    for f in file_list:
        if os.path.exists(f):
            print "%s: file exists, skipping download" % f
        else:
            download_file(os.path.join(url,f),f)
        if check_md5sum(f,chksums[f]):
            print "%s: checksum OK" % f
        else:
            print "%s: checksum FAILED" % f
            checksum_errors = True
    if checksum_errors:
        sys.stderr.write("ERROR one or more checksums failed, see above\n")
        sys.exit(1)
    else:
        print "All checksums verified ok"

