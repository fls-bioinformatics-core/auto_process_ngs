#!/usr/bin/env python
#
# QC plot generation
import os
import tempfile
from math import ceil
from matplotlib import pyplot as plt
from PIL import Image
from bcftbx.htmlpagewriter import PNGBase64Encoder
from .fastqc import FastqcData
from .fastqc import FastqcSummary
from .fastq_screen import Fastqscreen
from .fastq_stats import FastqQualityStats

# Colours taken from http://www.rapidtables.com/web/color/RGB_Color.htm
RGB_COLORS = {
    'black': (0,0,0),
    'blue': (0,0,255),
    'cornflowerblue': (100,149,237),
    'cyan': (0,255,255),
    'darkyellow1': (204,204,0),
    'green': (0,128,0),
    'grey': (145,145,145),
    'lightgrey': (211,211,211),
    'maroon': (128,0,0),
    'navyblue': (0,0,153),
    'orange': (255,165,0),
    'red': (255,0,0),
    'yellow': (255,255,0),
    'white': (255,255,255),
}

HEX_COLORS = {
    'green': '#008000',
    'orange': '#FFA500',
    'red': '#FF0000',
}

def encode_png(png_file):
    """
    Return Base64 encoded string for a PNG
    """
    return "data:image/png;base64," + \
        PNGBase64Encoder().encodePNG(png_file)

def screenplot(screen_files,outfile,threshold=None):
    """
    Generate plot of FastqScreen outputs

    Arguments:
      screen_files (list): list of paths to one or more
        ...screen.txt files from FastqScreen
      outfile (str): path to output file,outfile
      threshold (float): minimum percentage of mapped
        reads (below which library is excluded)

    """
    # Read in the screen data
    screens = []
    for screen_file in screen_files:
        screens.append(Fastqscreen(screen_file))
    nscreens = len(screens)
    # Plot the data
    plt.figure(1)
    for i,screen_data in enumerate(screens):
        # Create a sub-plot
        plt.subplot(nscreens,1,i+1)
        # Filter on threshold
        if threshold:
            screen_data = filter(lambda x: float(x['%Unmapped'])
                                 <= (100.0-threshold),
                                 screen_data)
        # Make a stacked bar chart
        plt.grid(True)
        x = xrange(len(screen_data))
        for mapping,color in (('%Multiple_hits_multiple_libraries','#800000'),
                              ('%One_hit_multiple_libraries','#000099'),
                              ('%Multiple_hits_one_library','r'),
                              ('%One_hit_one_library','b'),):
            data = [r[mapping] for r in screen_data]
            plt.barh(x,data,color=color,label=mapping)
        # Add the library names
        plt.yticks([x+0.5 for x in xrange(len(screen_data))],
                   [r['Library'] for r in screen_data])
        # Set axis limits on current plot
        plt.gca().set_xlim([0,100])
        plt.gca().set_ylim([0,len(screen_data)])
        # Only set legend for last plot
        if i == 0:
            plt.legend(loc=4)
    # Write out plot
    plt.savefig(outfile)
    
def uscreenplot(screen_files,outfile=None,inline=None):
    """
    Generate 'micro-plot' of FastqScreen outputs

    Arguments:
      screen_files (list): list of paths to one or more
        ...screen.txt files from FastqScreen
      outfile (str): path to output file

    """
    # Mappings
    mappings = ('%One_hit_one_library',
                '%Multiple_hits_one_library',
                '%One_hit_multiple_libraries',
                '%Multiple_hits_multiple_libraries',)
    # Colours
    colors = (RGB_COLORS['blue'],
              RGB_COLORS['navyblue'],
              RGB_COLORS['red'],
              RGB_COLORS['maroon'])
    # Read in the screen data
    screens = []
    for screen_file in screen_files:
        screens.append(Fastqscreen(screen_file))
    nscreens = len(screens)
    # Make a small stacked bar chart
    bbox_color = (145,145,145)
    barwidth = 4
    width = nscreens*50
    n_libraries_max = max([len(s) for s in screens])
    height = (n_libraries_max + 1)*(barwidth + 1)
    img = Image.new('RGB',(width,height),"white")
    pixels = img.load()
    # Process each screen in turn
    for nscreen,screen in enumerate(screens):
        xorigin = nscreen*50
        xend = xorigin+50-1
        yend = height-1
        # Draw a box around the plot
        for i in xrange(xorigin,xorigin+50):
            pixels[i,0] = bbox_color
            pixels[i,yend] = bbox_color
        for j in xrange(height):
            pixels[xorigin,j] = bbox_color
            pixels[xend,j] = bbox_color
        # Draw the stacked bars for each library
        for n,library in enumerate(screen.libraries):
            data = filter(lambda x: x['Library'] == library,screen)[0]
            x = xorigin
            y = n*(barwidth+1) + 1
            # Get the total percentage for the stack
            total_percent = sum([data[m] for m in mappings])
            if total_percent > 2.0:
                # Plot the stack as-is
                for mapping,rgb in zip(mappings,colors):
                    # Round up to nearest pixel (so that non-zero
                    # percentages are always represented)
                    npx = int(ceil(data[mapping]/2.0))
                    for i in xrange(x,x+npx):
                        for j in xrange(y,y+barwidth):
                            pixels[i,j] = rgb
                    x += npx
            elif total_percent > 0.25:
                # Small non-zero values can't be represented
                # accurately so just plot a placeholder
                max_mapped = 0.0
                for mapping,rgb in zip(mappings,colors):
                    if data[mapping] > max_mapped:
                        max_rgb = rgb
                for j in xrange(y,y+barwidth):
                    pixels[xorigin,j] = max_rgb
        # Add 'no hits'
        x = xorigin
        y = n_libraries_max*(barwidth+1) + 1
        npx = int(screen.no_hits/2.0)
        for i in xrange(x,x+npx):
            for j in xrange(y,y+barwidth):
                pixels[i,j] = bbox_color
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ufastqscreen.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        os.rename(tmp_plot,outfile)
    os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def uboxplot(fastqc_data=None,fastq=None,
             outfile=None,inline=None):
    """
    Generate FASTQ per-base quality 'micro-boxplot'

    'Micro-boxplot' is a thumbnail version of the per-base
    quality boxplots for a FASTQ file.

    Arguments:
       fastqc_data (str): path to a ``fastqc_data.txt``
        file
       outfile (str): path to output file

    Returns:
       String: path to output PNG file

    """
    # Boxplots need: mean, median, 25/75th and 10/90th quantiles
    # for each base
    max_qual = 41
    height = max_qual + 1
    fastq_stats = FastqQualityStats()
    if fastqc_data is not None:
        fastq_stats.from_fastqc_data(fastqc_data)
    elif fastq is not None:
        fastq_stats.from_fastq(fastq)
    else:
        raise Exception("supply path to fastqc_data.txt or fastq file")
    # To generate a bitmap in Python see:
    # http://stackoverflow.com/questions/20304438/how-can-i-use-the-python-imaging-library-to-create-a-bitmap
    #
    # Initialise output image instance
    img = Image.new('RGB',(fastq_stats.nbases,height),"white")
    pixels = img.load()
    # Create colour bands for different quality ranges
    for i in xrange(0,fastq_stats.nbases,2):
        for j in xrange(0,20):
            pixels[i,max_qual-j-1] = (230,175,175)
        for j in xrange(20,30):
            pixels[i,max_qual-j-1] = (230,215,175)
        for j in xrange(30,max_qual):
            pixels[i,max_qual-j-1] = (175,230,175)
    # Draw a box around the outside
    box_color = RGB_COLORS['grey']
    for i in xrange(fastq_stats.nbases):
        pixels[i,0] = box_color
        pixels[i,height-1] = box_color
    for j in xrange(height):
        pixels[0,j] = box_color
        pixels[fastq_stats.nbases-1,j] = box_color
    # For each base position determine stats
    for i in xrange(fastq_stats.nbases):
        #print "Position: %d" % i
        for j in xrange(fastq_stats.p10[i],fastq_stats.p90[i]):
            # 10th-90th percentile coloured cyan
            pixels[i,max_qual-j] = RGB_COLORS['grey']
        for j in xrange(fastq_stats.q25[i],fastq_stats.q75[i]):
            # Interquartile range coloured yellow
            pixels[i,max_qual-j] = RGB_COLORS['darkyellow1']
        # Median coloured red
        pixels[i,max_qual-int(fastq_stats.median[i])] = RGB_COLORS['red']
        # Mean coloured black
        pixels[i,max_qual-int(fastq_stats.mean[i])] = RGB_COLORS['blue']
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".uboxplot.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        os.rename(tmp_plot,outfile)
    os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def ufastqcplot(summary_file,outfile=None,inline=False):
    """
    Make a 'micro' summary plot of FastQC output

    The micro plot is a small PNG which represents the
    summary results from each FastQC module in a
    matrix, with rows representing the modules and
    three columns representing the status ('PASS', 'WARN'
    and 'FAIL', from left to right).

    For example (in text form):

          ==
    ==
    ==
       ==
    ==

    indictaes that the status of the first module is
    'FAIL', the 2nd, 3rd and 5th are 'PASS', and the
    4th is 'WARN'.

    Arguments:
      summary_file (str): path to a FastQC
        'summary.txt' output file
      outfile (str): path for the output PNG

    """
    status_codes = {
        'PASS' : { 'index': 0,
                   'color': 'green',
                   'rgb': RGB_COLORS['green'],
                   'hex': HEX_COLORS['green'] },
        'WARN' : { 'index': 1,
                   'color': 'orange',
                   'rgb': RGB_COLORS['orange'],
                   'hex': HEX_COLORS['orange'] },
        'FAIL' : { 'index': 2,
                   'color': 'red',
                   'rgb': RGB_COLORS['red'],
                   'hex': HEX_COLORS['red'] },
        }
    fastqc_summary = FastqcSummary(summary_file)
    # Initialise output image instance
    nmodules = len(fastqc_summary.modules)
    img = Image.new('RGB',(30,4*nmodules),"white")
    pixels = img.load()
    # For each test: put a mark depending on the status
    for im,m in enumerate(fastqc_summary.modules):
        code = status_codes[fastqc_summary.status(m)]
        # Make the mark
        x = code['index']*10 + 1
        #y = 4*nmodules - im*4 - 3
        y = im*4 + 1
        for i in xrange(x,x+8):
            for j in xrange(y,y+3):
                #print "%d %d" % (i,j)
                pixels[i,j] = code['rgb']
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ufastqc.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        os.rename(tmp_plot,outfile)
    os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def ustackedbar(data,outfile=None,inline=False,bbox=True,
                height=20,length=100,colors=None):
    """
    Make a 'micro' stacked bar chart

    A 'stacked' bar consists of a bar divided into
    sections, with each section of proportional length
    to the corresponding value.

    Arguments:
      data (List): list or tuple of data values
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
      bbox (boolean): if True then draw a bounding box
        around the plot
      height (int): height of the bar in pixels
      length (int): length of the bar in pixels
      colors (List): list or tuple of color values
    """
    # Initialise
    bgcolor = "black"
    if colors is None:
        colors = sorted(list(RGB_COLORS.keys()))
    # Create the image
    img = Image.new('RGB',(length,height),bgcolor)
    pixels = img.load()
    # Normalise the data
    total = float(sum(data))
    ndata = [int(float(d)/total*float(length)) for d in data]
    # Reset the last value
    ndata[-1] = length - sum(ndata[:-1])
    # Create the plot
    p = 0
    for ii,d in enumerate(ndata):
        color = colors[ii%len(colors)]
        try:
            color = RGB_COLORS[color]
        except KeyError:
            pass
        for i in xrange(p,p+d):
             for j in xrange(height):
                pixels[i,j] = color
        p += d
    # Overlay a bounding box
    if bbox:
        for i in xrange(0,length):
            pixels[i,0] = RGB_COLORS[bgcolor]
            pixels[i,height-1] = RGB_COLORS[bgcolor]
        for j in xrange(0,height):
            pixels[0,j] = RGB_COLORS[bgcolor]
            pixels[length-1,j] = RGB_COLORS[bgcolor]
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ubar.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        os.rename(tmp_plot,outfile)
    os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile
