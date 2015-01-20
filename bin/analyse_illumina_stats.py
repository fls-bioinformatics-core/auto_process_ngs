#!/bin/env python
import sys
import bcftbx.TabFile as tf
from bcftbx.IlluminaData import IlluminaFastq

if __name__ == '__main__':
    # Read in data
    stats = tf.TabFile(sys.argv[1],first_line_is_header=True)
    data = dict()
    for line in stats:
        # Collect sample name, lane etc
        fq = IlluminaFastq(line['Fastq'])
        if fq.read_number != 1:
            # Only interested in R1 reads
            continue
        lane = fq.lane_number
        sample = "%s/%s" % (line['Project'],line['Sample'])
        nreads = line['Nreads']
        # Update information in dict
        if lane not in data:
            data[lane] = dict()
        try:
            data[lane][sample] += nreads
        except KeyError:
            data[lane][sample] = nreads
    # Get list of lanes
    lanes = [int(x) for x in data.keys()]
    lanes.sort()
    # Report
    for lane in lanes:
        print "\nLane %d" % lane
        samples = data[lane].keys()
        samples.sort()
        total_reads = sum([data[lane][x] for x in samples])
        print "Total reads = %d" % total_reads
        for sample in samples:
            nreads = data[lane][sample]
            print "- %s\t%d\t%.1f%%" % (sample,nreads,float(nreads)/total_reads*100.0)
    
            
