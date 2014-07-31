#!/bin/env python
#
# Prepare an entry for the logging file
#
# Iterates over subdirectories of the supplied top-level directory
# looking for README files containing information about each project,
# extracts the data and creates a description suitable for the logging
# file
import sys
import os
if __name__ == "__main__":
    try:
        top_dir = sys.argv[1]
    except IndexError:
        sys.stderr.write("Need to supply a directory\n")
        sys.exit(1)
    logging_entry = []
    print "Examining %s" % top_dir
    for f in os.listdir(top_dir):
        print "Looking at %s" % f
        readme = os.path.join(top_dir,f,'README')
        if os.path.isfile(readme):
            data = dict()
            fp = open(readme,'rU')
            for line in fp:
                items = line.rstrip().split('\t')
                ##print "%s" % items
                try:
                    data[items[0]] = items[1]
                except IndexError:
                    pass
            try:
                logging_entry.append("%s: %s %s %s (PI: %s) [%s] %s" % (f,
                                                                        data['User'],
                                                                        data['Organism'],
                                                                        data['Library type'],
                                                                        data['PI'],
                                                                        data['Library prep'],
                                                                        data['Samples']))
            except KeyError:
                pass
    print '; '.join(logging_entry)
