#!/usr/bin/env python
#
# QC plot generation
import os
import shutil
import tempfile
import logging
from math import ceil
from PIL import Image
from builtins import range
from bcftbx.htmlpagewriter import PNGBase64Encoder
from .fastqc import FastqcData
from .fastqc import FastqcSummary
from .fastq_screen import Fastqscreen
from .fastq_stats import FastqQualityStats
from .fastq_strand import Fastqstrand

# Module specific logger
logger = logging.getLogger(__name__)

# Colours taken from http://www.rapidtables.com/web/color/RGB_Color.htm
RGB_COLORS = {
    'black': (0,0,0),
    'blue': (0,0,255),
    'cornflowerblue': (100,149,237),
    'cyan': (0,255,255),
    'darkyellow1': (204,204,0),
    'gainsboro': (220,220,220), # Very light grey!
    'green': (0,128,0),
    'grey': (145,145,145),
    'lightgrey': (211,211,211),
    'lightblue': (173,216,230),
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
    return "data:image/png;base64,%s" % \
        PNGBase64Encoder().encodePNG(png_file)
    
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
        for i in range(xorigin,xorigin+50):
            pixels[i,0] = bbox_color
            pixels[i,yend] = bbox_color
        for j in range(height):
            pixels[xorigin,j] = bbox_color
            pixels[xend,j] = bbox_color
        # Draw the stacked bars for each library
        for n,library in enumerate(screen.libraries):
            data = list(filter(lambda x:
                               x['Library'] == library,screen))[0]
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
                    # Don't exceed plot limit
                    if x+npx > width:
                        npx = width-x
                    for i in range(x,x+npx):
                        for j in range(y,y+barwidth):
                            pixels[i,j] = rgb
                    x += npx
            elif total_percent > 0.25:
                # Small non-zero values can't be represented
                # accurately so just plot a placeholder
                max_mapped = 0.0
                for mapping,rgb in zip(mappings,colors):
                    if data[mapping] > max_mapped:
                        max_rgb = rgb
                for j in range(y,y+barwidth):
                    pixels[xorigin,j] = max_rgb
        # Add 'no hits'
        x = xorigin
        y = n_libraries_max*(barwidth+1) + 1
        npx = int(screen.no_hits/2.0)
        for i in range(x,x+npx):
            for j in range(y,y+barwidth):
                pixels[i,j] = bbox_color
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ufastqscreen.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
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
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file

    Returns:
       String: path to output PNG file
    """
    # Boxplots need: mean, median, 25/75th and 10/90th quantiles
    # for each base
    max_qual = 41
    fastq_stats = FastqQualityStats()
    if fastqc_data is not None:
        try:
            fastq_stats.from_fastqc_data(fastqc_data)
        except Exception as ex:
            logger.warning("uboxplot: failed to load Fastqc data: %s" % ex)
            raise ex
    elif fastq is not None:
        fastq_stats.from_fastq(fastq)
    else:
        raise Exception("supply path to fastqc_data.txt or fastq file")
    # Sweep the data to check for the maximum quality score
    for stats in (fastq_stats.p10,
                  fastq_stats.q25,
                  fastq_stats.q75,
                  fastq_stats.p90,
                  fastq_stats.median,
                  fastq_stats.mean):
        for qual in stats:
            if qual > max_qual:
                max_qual = int(ceil(qual))
                logger.warning("uboxplot: setting max quality to %d" %
                               max_qual)
    # To generate a bitmap in Python see:
    # http://stackoverflow.com/questions/20304438/how-can-i-use-the-python-imaging-library-to-create-a-bitmap
    #
    # Initialise output image instance
    height = max_qual + 1
    img = Image.new('RGB',(fastq_stats.nbases,height),"white")
    pixels = img.load()
    # Create colour bands for different quality ranges
    for i in range(0,fastq_stats.nbases,2):
        for j in range(0,20):
            pixels[i,max_qual-j-1] = (230,175,175)
        for j in range(20,30):
            pixels[i,max_qual-j-1] = (230,215,175)
        for j in range(30,max_qual):
            pixels[i,max_qual-j-1] = (175,230,175)
    # Draw a box around the outside
    box_color = RGB_COLORS['grey']
    for i in range(fastq_stats.nbases):
        pixels[i,0] = box_color
        pixels[i,height-1] = box_color
    for j in range(height):
        pixels[0,j] = box_color
        pixels[fastq_stats.nbases-1,j] = box_color
    # For each base position determine stats
    for i in range(fastq_stats.nbases):
        #print("Position: %d" % i)
        try:
            for j in range(fastq_stats.p10[i],fastq_stats.p90[i]):
                # 10th-90th percentile coloured cyan
                pixels[i,max_qual-j] = RGB_COLORS['grey']
        except TypeError:
            pass
        try:
            for j in range(fastq_stats.q25[i],fastq_stats.q75[i]):
                # Interquartile range coloured yellow
                pixels[i,max_qual-j] = RGB_COLORS['darkyellow1']
        except TypeError:
            pass
        # Median coloured red
        try:
            median = int(fastq_stats.median[i])
            pixels[i,max_qual-median] = RGB_COLORS['red']
        except TypeError:
            pass
        # Mean coloured black
        pixels[i,max_qual-int(fastq_stats.mean[i])] = RGB_COLORS['blue']
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".uboxplot.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
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

    indicates that the status of the first module is
    'FAIL', the 2nd, 3rd and 5th are 'PASS', and the
    4th is 'WARN'.

    Arguments:
      summary_file (str): path to a FastQC
        'summary.txt' output file
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
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
        for i in range(x,x+8):
            for j in range(y,y+3):
                #print("%d %d" % (i,j))
                pixels[i,j] = code['rgb']
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ufastqc.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
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
    try:
        ndata = [int(float(d)/total*float(length)) for d in data]
    except ZeroDivisionError:
        # Total was zero
        ndata = [0 for d in data]
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
        for i in range(p,p+d):
             for j in range(height):
                pixels[i,j] = color
        p += d
    # Overlay a bounding box
    if bbox:
        for i in range(0,length):
            pixels[i,0] = RGB_COLORS[bgcolor]
            pixels[i,height-1] = RGB_COLORS[bgcolor]
        for j in range(0,height):
            pixels[0,j] = RGB_COLORS[bgcolor]
            pixels[length-1,j] = RGB_COLORS[bgcolor]
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ubar.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
        os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def ustrandplot(fastq_strand_out,outfile=None,inline=False,
                height=25,width=50,fg_color=None,dynamic=False):
    """
    Make a 'micro' chart for strandedness


    This micro plot is a small PNG which summarises the
    results from fastq_strand.py as two horizontal bars
    (one for forward, one for reverse) with the lengths
    representing the psuedo-percentages of each.

    For example (in text form):

    =
    ==========

    If the fastq_strand results included multiple genomes
    then there will be one pair of bars for each genome.

    Arguments:
      fastq_strand_out (str): path to a fastq_strand
        output file
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
      height (int): height of the plot in pixels
      width (int): width of the plot in pixels
      fg_color (tuple): tuple of RGB values to use for
        the foreground colour of the bars
      dynamic (boolean): if True then the height of the
        plot will be increased for each additional
        genome in the output fastq_strand file
    """
    # Get the raw data
    data = Fastqstrand(fastq_strand_out)
    # Adjust the plot height for "dynamic" mode
    ngenomes = len(data.genomes)
    if dynamic and ngenomes:
        height *= ngenomes
    # Set the largest percentage from the data
    max_percent = 100
    for genome in data.genomes:
        for p in ('forward','reverse'):
            max_percent = max(max_percent,data.stats[genome][p])
    # Set the spacing between bars in the plot
    spacing = 3
    # Set the width of the bars
    if ngenomes:
        bar_width = int(float(float(height)/ngenomes - 3*spacing)/2.0)
    # Colour
    if fg_color is None:
        fg_color = RGB_COLORS['black']
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS['white'])
    pixels = img.load()
    # Plot bars for the forward and reverse percentages
    # for each genome
    for ii,genome in enumerate(data.genomes):
        # Forward strand
        bar_length = int(data.stats[genome].forward/
                         max_percent*(width-4))
        for i in range(2,bar_length+2):
            start = int(ii*float(height)/ngenomes) + spacing
            end = start + bar_width
            for j in range(start,end):
                pixels[i,j] = fg_color
        # Pad the remainder of the bar
        for i in range(bar_length+2,width-2):
            start = int(ii*float(height)/ngenomes) + spacing
            end = start + bar_width
            for j in range(start,end):
                pixels[i,j] = RGB_COLORS['lightgrey']
        # Reverse strand
        bar_length = max(int(data.stats[genome].reverse/
                         max_percent*(width-4)),1)
        for i in range(2,bar_length+2):
            start = int((float(ii)+0.5)*float(height)/ngenomes) + spacing
            end = start + bar_width
            for j in range(start,end):
                pixels[i,j] = fg_color
        # Pad the remainder of the bar
        for i in range(bar_length+2,width-2):
            start = int((float(ii)+0.5)*float(height)/ngenomes) + spacing
            end = start + bar_width
            for j in range(start,end):
                pixels[i,j] = RGB_COLORS['lightgrey']
    # Output the plot to file
    fp,tmp_plot = tempfile.mkstemp(".ustrand.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
        os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def useqlenplot(dist,masked_dist=None,min_len=None,max_len=None,
                outfile=None,inline=False,height=None,
                bg_color="gainsboro",bbox_color="white",
                seq_color="black",masked_color="red"):
    """
    Make a 'micro' plot of sequence length

    Given a sequence length distribution, create a histogram-style
    plot where the numbers of sequences with different lengths
    are shown (similar to the 'Sequence Length Distribution' plot
    from FastQC).

    Optionally if a distribution of masked reads is also supplied
    then these data will be overlayed on top.

    The distributions should be supplied as dictionaries or
    mappings where the keys are sequence lengths and the
    corresponding values are the number of sequences.

    Arguments:
      dist (mapping): mapping of sequence lengths to numbers
        of sequences, giving the distribution of sequence
        lengths
      masked_dist (mapping): optional, mapping of sequence
        lengths to numbers of masked reads
      min_len (int): optional, set the lower limit of the
        plot (otherwise defaults to the lowest length
        present in the distribution)
      max_len (int): optional, set the upper limit of the
        plot (otherwise defaults to the highest length
        present in the distribution)
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
      height (int): height of the plot in pixels
      bg_color (str): name of colour to use for the
        background
      bbox_color (str): name of colour to use for the
        bounding box
      seq_color (str): name of colour to use for the
        sequence distribution
      masked_color (str): name of color to use for masked
        sequence distribution
    """
    # Set height of plot
    if height is None:
        height = 20
    # Get ordered list of sequence lengths
    seq_lens = sorted(dist.keys())
    # Get max number of reads for a single length
    max_reads = max([dist[x] for x in dist])
    # Get min and max read lengths & size of plot
    if max_len is None:
        max_len = max(dist)
    if min_len is None:
        min_len = min(dist)
    width = max_len - min_len + 3
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS[bg_color])
    pixels = img.load()
    # Add bounding box
    for i in range(0,width):
        pixels[i,0] = RGB_COLORS[bbox_color]
        pixels[i,height-1] = RGB_COLORS[bbox_color]
    for j in range(0,height):
        pixels[0,j] = RGB_COLORS[bbox_color]
        pixels[width-1,j] = RGB_COLORS[bbox_color]
    # Create the plot
    for seq_len in seq_lens:
        # Plot the sequencing lengths
        i = seq_len - min_len
        if dist[seq_len]:
            ndata = max(1,int(float(dist[seq_len])/max_reads*(height-2)))
            for j in range(1,ndata+1):
                pixels[i+1,height-j-1] = RGB_COLORS[seq_color]
        # Overlay the masked read data
        if masked_dist:
            try:
                ndata = max(1,int(float(masked_dist[seq_len])/
                                  max_reads*(height-2)))
                for j in range(1,ndata+1):
                    pixels[i+1,height-j-1] = RGB_COLORS[masked_color]
            except KeyError:
                pass
    # Output the plot
    return make_plot(img,
                     outfile=outfile,
                     inline=inline,
                     ext=".useqlen.png")

def ureadcountplot(nreads,nmasked=None,npadded=None,max_reads=None,
                   outfile=None,inline=False,width=50,height=6,
                   bg_color="white",fg_color="green",masked_color="red",
                   padded_color="orange",fill_color="lightgrey"):
    """
    Make a 'micro' plot summarising read counts and masking

    Given a total number of reads/sequences in a Fastq file
    plus the number of those sequences which are masked (i.e.
    completely composed of Ns) and padded (i.e. have one or
    more trailing Ns), plots a horizontal bar indicating the
    sequence composition.

    By default the numbers are normalised so that the total
    number of reads fills the bar; however if a maximum read
    count is also supplied then the normalisation is relative
    to that maximum (so the plot also indicates the relative
    size of the Fastq compared to the maximum read count).

    Arguments:
      nreads (int): number of reads in the Fastq
      nmasked (int): number of masked reads
      npadded (int): number of padded reads
      max_reads (int): maximum number of reads (e.g. in
        all Fastqs) for normalisation
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
      width (int): width of the plot in pixels
      height (int): height of the plot in pixels
      bg_color (str): name of colour to use for the
        background
      fg_color (str): name of the colour for plotting
        unmasked, unpadded reads
      masked_color (str): name of the colour for plotting
        masked read fraction
      padded_color (str): name of the colour for plotting
        padded read fraction
      fill_color (str): name of the colour for filling
        the unoccupied remainder of the bar
    """
    # Initialise
    if not max_reads:
        max_reads = nreads
    if not nmasked:
        nmasked = 0
    if not npadded:
        npadded = 0
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS[bg_color])
    pixels = img.load()
    # Plot masked content
    start = 0
    if nmasked and nmasked > max_reads/100.0:
        end = max(1,int(float(nmasked)/float(max_reads)*width))
        for i in range(start,end):
            for j in range(1,height-1):
                pixels[i,j] = RGB_COLORS[masked_color]
        start = end
    else:
        end = start
    # Plot padded content
    if npadded and npadded > max_reads/100.0:
        end = max(end+1,end+int(float(npadded)/float(max_reads)*width))
        for i in range(start,end):
            for j in range(1,height-1):
                pixels[i,j] = RGB_COLORS[padded_color]
        start = end
    # Plot remaining read content
    if nreads == max_reads:
        end = width
    else:
        end = int(float(nreads)/float(max_reads)*width)
    for i in range(start,end):
        for j in range(1,height-1):
            pixels[i,j] = RGB_COLORS[fg_color]
    # Fill remainder of the plot
    for i in range(end,width):
        for j in range(1,height-1):
            pixels[i,j] = RGB_COLORS[fill_color]
    # Output the plot
    return make_plot(img,
                     outfile=outfile,
                     inline=inline,
                     ext=".ureadcount.png")

def uadapterplot(adapter_content,adapter_names=None,outfile=None,
                 inline=False,height=25,bar_width=6,spacing=2,
                 multi_bar=True):
    """
    Make a 'micro' plot summarising adapter content

    The plot consists of vertical bar(s) which indicate the
    relative presence of each adapter  class in the sequence
    data.

    In 'multi-bar' mode the plot has one bar for each
    adapter class; in 'single-bar' mode the plot combines
    all data into a single bar.

    The adapter content should be supplied as a dictionary
    where the keys are adapter names and the corresponding
    adapter content is expressed as a decimal fraction.

    Arguments:
      adapter_content (mapping): dictionary mapping
        adapter names to adapter content
      adapter_names (list): optional, list of adapter
        classes; if provided then defines the order in
        which the adapters appear in the plot
        (otherwise taken from the keys in the mapping)
      outfile (str): path for the output PNG
      inline (bool): if True then returns the PNG
        as base64 encoded string rather than as a file
      height (int): height of the plot in pixels
      bar_width (int): width of each bar representing
        content for an adapter class, in pixels
      spacing (int): spacing between each bar, in pixels
      multi_bar (bool): if True then make a multi-bar
        plot (one bar per adapter class); otherwise
        make a single bar plot (all adapter data in a
        single bar)
    """
    # Width of plot required for each bar
    width = bar_width + spacing*2
    # Colours for each adapter
    fg_colors = ('red','blue','green','black')
    # Get adapter names
    if not adapter_names:
        adapter_names = sorted(adapter_content.keys())
    # Number of bars and total width of plot
    if multi_bar:
        nbars = len(adapter_names)
    else:
        nbars = 1
    width *= nbars
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS['white'])
    pixels = img.load()
    # Add a baseline on the plot
    for j in range(spacing,width-spacing):
        pixels[j,height-2] = RGB_COLORS['lightgrey']
    # Generate the plot
    if multi_bar:
        # Plot a bar for each adapter
        for ii,adapter in enumerate(adapter_names):
            # Set colour based on adapter content
            fg_color = fg_colors[ii%4]
            # Length of bar represents adapter content
            bar_length = int(adapter_content[adapter]*(height-4))
            # Draw the coloured part of the bar
            for i in range(2,bar_length+2):
                start = int(ii*float(width)/nadapters) + spacing
                end = start + bar_width
                for j in range(start,end):
                    pixels[j,height-i] = RGB_COLORS[fg_color]
            # Pad the remainder of the bar
            for i in range(bar_length+2,height-2):
                start = int(ii*float(width)/nadapters) + spacing
                end = start + bar_width
                for j in range(start,end):
                    pixels[j,height-i] = RGB_COLORS['lightgrey']
    else:
        # Single bar with all data
        start = 1
        for ii,adapter in enumerate(adapter_names):
            # Set colour based on adapter content
            fg_color = fg_colors[ii%4]
            # Length of bar represents adapter content
            bar_length = int(adapter_content[adapter]*(height-2))
            # Draw the coloured part of the bar
            end = start + bar_length
            for i in range(start,end):
                for j in range(1,width-2):
                    pixels[j,height-i] = RGB_COLORS[fg_color]
            start = end
        # Pad the remainder of the bar
        for i in range(start,height-2):
            for j in range(1,width-2):
                pixels[j,height-i] = RGB_COLORS['lightgrey']
    # Output the plot
    return make_plot(img,
                     outfile=outfile,
                     inline=inline,
                     ext=".uadapters.png")

def make_plot(img,outfile=None,inline=False,ext=".plot.png"):
    """
    Internal: output PNG plots from Image objects

    Arguments:
      img (Image): image to output plot for
      outfile (str): path to output file to write
        PNG to (if None then no file will be
        created)
      inline (bool): if True then return base64
        encoded string for the plot
      ext (str): optional, extension to use for
        temporary file
    """
    fp,tmp_plot = tempfile.mkstemp(ext)
    img.save(tmp_plot)
    os.fdopen(fp).close()
    if inline:
        encoded_plot = encode_png(tmp_plot)
    if outfile is not None:
        shutil.move(tmp_plot,outfile)
    else:
        os.remove(tmp_plot)
    if inline:
        return encoded_plot
    else:
        return outfile

def _tiny_png(outfile,width=4,height=4,
             bg_color=RGB_COLORS['white'],
             fg_color=RGB_COLORS['blue']):
    """
    Create a small checkered PNG for testing
    """
    img = Image.new('RGB',(width,height),bg_color)
    pixels = img.load()
    for i in range(int(width/2)):
        for j in range(int(height/2)):
            pixels[i,j] = fg_color
    for i in range(int(width/2),width):
        for j in range(int(height/2),height):
            pixels[i,j] = fg_color
    fp,tmp_plot = tempfile.mkstemp(".tiny.png")
    img.save(tmp_plot)
    os.fdopen(fp).close()
    shutil.move(tmp_plot,outfile)
    return outfile
