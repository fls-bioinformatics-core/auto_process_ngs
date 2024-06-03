#!/usr/bin/env python
#
# QC plot generation
import os
import shutil
import tempfile
import logging
import collections
from math import ceil
from math import floor
from PIL import Image
from builtins import range
from bcftbx.htmlpagewriter import PNGBase64Encoder
from bcftbx.utils import AttributeDictionary
from .fastqc import FastqcData
from .fastqc import FastqcSummary
from .fastq_screen import Fastqscreen
from .fastq_stats import FastqQualityStats
from .fastq_strand import Fastqstrand

# Data
from .constants import SEQUENCE_DEDUP_CUTOFFS

# Module specific logger
logger = logging.getLogger(__name__)

# Colours taken from http://www.rapidtables.com/web/color/RGB_Color.htm
RGB_COLORS = AttributeDictionary(
    black=(0,0,0),
    blue=(0,0,255),
    cornflowerblue=(100,149,237),
    cyan=(0,255,255),
    darkyellow1=(204,204,0),
    gainsboro=(220,220,220), # Very light grey!
    green=(0,128,0),
    grey=(145,145,145),
    lightgrey=(211,211,211),
    lightblue=(173,216,230),
    maroon=(128,0,0),
    navyblue=(0,0,153),
    orange=(255,165,0),
    red=(255,0,0),
    yellow=(255,255,0),
    white=(255,255,255),
)

HEX_COLORS = {
    'green': '#008000',
    'orange': '#FFA500',
    'red': '#FF0000',
}

class Plot:
    """
    Utility class for creating pixel-based microplots

    Provides methods for building up plots out of smaller
    elements (lines, blocks, single points and sets of
    points). It also provides methods for adding
    bounding boxes and background striping.

    Usage:

    >>> p = Plot(25,25)
    >>> p.bbox(RGB.grey)
    >>> p.plot({ x:x for x in range(25)},RGB.red)
    >>> p.save("example.png")
    >>> p.encoded_png()
    'data:shjdhbcchb...'

    Arguments:
      width (int): width of the plot canvas (pixels)
      height (int): height of the plot canvas (pixels)
      bgcolor (str): background color (default: 'white')
    """
    def __init__(self,width,height,bgcolor='white'):
        self._width = int(width)
        self._height = int(height)
        self._bgcolor = bgcolor
        self._pixels = dict()
        for x in range(self._width):
            self._pixels[x] = dict()

    def set_pixel(self,x,y,color):
        """
        Set the color of a single pixel

        Arguments:
          x (int): x-position of pixel
          y (int): y-position of pixel
          color (tuple): tuple specifying an RGB color
        """
        x = self._int(x)
        y = self._int(y)
        try:
            if (0 <= y < self._height):
                self._pixels[x][y] = color
            else:
                raise KeyError
        except KeyError:
            print("%s,%s: out-of-bounds" % (x,y))

    def hline(self,y,color):
        """
        Draw a horizontal line

        Argument:
          y (int): y-axis position of line
          color (tuple): tuple specifying an RGB color
        """
        for x in range(self._width):
            self.set_pixel(x,y,color)

    def vline(self,x,color,llen=None):
        """
        Draw a vertical line

        Argument:
          y (int): y-axis position of line
          color (tuple): tuple specifying an RGB color
          llen (int): length of the line (defaults to
            plot canvas height)
        """
        if llen is None:
            llen = self._height
        llen = self._int(llen)
        for y in range(llen):
            self.set_pixel(x,y,color)

    def block(self,xy1,xy2,color,color2=None):
        """
        Draw a rectangular block

        The block is defined by corners xy1 and xy2.

        Arguments:
          xy1 (tuple): pair of (x,y) positions defining
            first corner of the block
          xy2 (tuple): pair of (x,y) positions defining
            second corner of the block
          color (tuple): tuple specifying an RGB color
            to fill the block with
        """
        x1 = self._int(min((xy1[0],xy2[0])))
        x2 = self._int(max((xy1[0],xy2[0])))
        y1 = self._int(min((xy1[1],xy2[1])))
        y2 = self._int(max((xy1[1],xy2[1])))
        for x in range(x1,x2):
            if color2 is not None and x%2 != 0:
                color_ = color2
            else:
                color_ = color
            for y in range(y1,y2):
                self.set_pixel(x,y,color_)

    def stripe(self,color1,color2):
        """
        Fill the plot with vertical stripes

        Arguments:
          color1 (tuple): tuple specifying first RGB color
          color2 (tuple): tuple specifying second RGB color
        """
        for x in range(self._width):
            if x%2 == 0:
                self.vline(x,color1)
            else:
                self.vline(x,color2)

    def bbox(self,color):
        """
        Draw a single-pixel bounding box

        Arguments:
          color (tuple): tuple specifying the border RGB
            color
        """
        self.hline(0,color)
        self.hline(self._height-1,color)
        self.vline(0,color)
        self.vline(self._width-1,color)

    def plot(self,data,color,fill=False,interpolation="mean"):
        """
        Plot arbitrary data along the x-axis

        If data limits exceed the canvas size then the data
        will be scaled in both x and y directions to fit into
        the plot canvas, with the 'interpolation' argument
        specifying the method to use for handling data from
        multiple points which are combined into a single bin:

        - 'mean' (the default) assigns the mean value of
          all the points that are put into the same bin
          (resulting in a single point plotted for each bin)
        - 'minmax' plots a vertical line for each bin based
          on the minimum and maximum values assigned to the
          bin

        (NB the 'minmax' interpolation can smooth out
        discontinuities in the plotted data that can appear
        when using the 'mean' method with datasets which are
        much wider than the plot width.)

        Setting the 'fill' argument to True when using the
        'mean' interpolation method fills the area under
        each plotted point; it is ignored when using 'minmax'.

        Arguments:
          data (mapping): set of x-axis positions mapping
            to corresponding y-axis positions
          color (tuple): tuple specifying an RGB color
          fill (bool): if True then also fill the
            points under each point (default: no fill)
          interpolation (str): interpolation mode to use
            when plotting rebinned data (can be one of
            'mean' or 'minmax')
        """
        # Check there is data
        if not data:
            logger.warning("Plot.plot(): empty dataset?")
            return
        # Rebin the data along the x-axis
        data_ = self.rebin_data(data)
        if interpolation == "mean":
            # Normalise the data along the y-axis
            data_ = self.normalise_data(data_.mean)
            # Plot data converted to integer values
            for x in data_:
                y = self._int(data_[x])
                if not fill:
                    # Plot single point
                    self.set_pixel(x,y,color)
                else:
                    # Fill area under plotted point
                    for yy in range(0,y):
                        self.set_pixel(x,yy,color)
        elif interpolation == "minmax":
            # Normalise the data along the y-axis
            maxval = max([data_.min[x] for x in data_.min]+
                         [data_.max[x] for x in data_.max])
            data1 = self.normalise_data(data_.min,maxval=maxval)
            data2 = self.normalise_data(data_.max,maxval=maxval)
            # Plot data ranges
            for x in data1:
                y1 = self._int(min(data1[x],data2[x]))
                y2 = self._int(max(data1[x],data2[x]))
                for y in range(y1,y2+1):
                    self.set_pixel(x,y,color)
        else:
             raise Exception("%s: unrecognised interpolation mode"
                             % interpolation)

    def plot_range(self,data1,data2,color):
        """
        Plot two sets of data and fill the region inbetween

        If data limits exceed the canvas size then the data
        will be scaled in both x and y directions to fit into
        the plot canvas.

        Arguments:
          data1 (mapping): set of x-axis positions mapping
            to corresponding y-axis positions for first
            dataset
          data2 (mapping): set of x-axis positions mapping
            to corresponding y-axis positions for second
            dataset
          color (tuple): tuple specifying an RGB color
        """
        # Plot two sets of data and fill the region
        # inbetween
        if not data1 or not data2:
            logger.warning("Plot.plot_range(): empty dataset(s)?")
            return
        # Rebin the data along the x-axis
        data1_ = self.rebin_data(data1).mean
        data2_ = self.rebin_data(data2).mean
        # Normalise the data along the y-axis
        data1_ = self.normalise_data(data1_)
        data2_ = self.normalise_data(data2_)
        # Fill region contained within two datasets
        for x in data1_:
            if x not in data2_:
                continue
            y1 = self._int(min(data1_[x],data2_[x]))
            y2 = self._int(max(data1_[x],data2_[x]))
            for y in range(y1,y2+1):
                self.set_pixel(x,y,color)

    def bar(self,data,xy1,xy2,colors):
        """
        Draw a horizontal stacked barplot

        The bar is defined by corners xy1 and xy2 (cf the
        'block' method).

        Arguments:
          data (sequence): list of values corresponding
            to the size of each section in the 'stack'
          xy1 (tuple): pair of (x,y) positions defining
            first corner of the bar
          xy2 (tuple): pair of (x,y) positions defining
            second corner of the bar
          colors (sequence) list of RGB colour tuples,
            to assign to each section in the 'stack';
            if the number of sections exceeds the number
            of colours then colours are reused in order
            from the beginning
        """
        # Extract limits of the bar
        x1 = xy1[0]
        y1 = xy1[1]
        x2 = xy2[0]
        y2 = xy2[1]
        # Bar length (pixels)
        bar_length = x2 - x1
        # Convert data to set of cumulative sums
        cumul = 0
        data_ = list()
        for d in data:
            cumul += d
            data_.append(cumul)
        # Convert sums to set of integer block lengths
        # and reverse order (so longest blocks are drawn
        # first, and subsequent blocks are overlaid)
        data_ = list(reversed([ int(float(d)/cumul*bar_length)
                                for d in data_ ]))
        # Number of data points
        ndata = len(data_)
        # Number of colours
        ncolors = len(colors)
        colors_ = list(reversed(colors))
        # Draw blocks on top of each other (longest first)
        for i,d in enumerate(data_):
            ii = (ndata - i - 1) % ncolors
            self.block((x1,y1),
                       (x1+d,y2),
                       colors[ii])

    def rebin_data(self,data):
        """
        Scales an arbitrary dataset along x-axis

        Scales the supplied dataset so that the x-axis
        data fits into the plot size, by treating the
        plot x-values as 'bins' and combining data points
        from the dataset that fall into the same bin.

        Arguments:
          data (mapping): set of x-axis positions mapping
            to corresponding y-axis positions

        Returns:
          NamedTuple: named tuple with three datasets
            referenced by the keys 'mean' (mean values),
            'min' (minimum values) and 'max' (maximum
            values) for each position.
        """
        # Get max x-value in the data
        xvals = list(data.keys())
        maxval = max(xvals)
        # Rebin data along x-axis
        data_ = dict()
        mean_data = dict()
        min_data = dict()
        max_data = dict()
        if maxval > self._width:
            for x in xvals:
                xx = self._int(float(x)/float(maxval)*(self._width-1))
                if xx not in data_:
                    data_[xx] = list()
                data_[xx].append(data[x])
            for x in data_:
                # Average data values in each bin
                mean_data[x] = sum(data_[x])/len(data_[x])
                # Min/max data values in each bin
                min_data[x] = min(data_[x])
                max_data[x] = max(data_[x])
        else:
            mean_data = { x:data[x] for x in data }
            min_data = { x:data[x] for x in data }
            max_data = { x:data[x] for x in data }
        # Return using a named tuple
        rebinneddata = collections.namedtuple("rebinneddata",
                                              "mean min max")
        return rebinneddata(mean=mean_data,
                            min=min_data,
                            max=max_data)

    def normalise_data(self,data,maxval=None):
        """
        Scales an arbitrary dataset along y-axis

        Arguments:
          data (mapping): set of x-axis positions mapping
            to corresponding y-axis positions

        Returns:
          Mapping: copy of the dataset with normalised
            y-values
        """
        # Get max y-value in the data
        if maxval is None:
            maxval = max([data[x] for x in data])
        # Normalise data along y-axis
        if maxval > self._height:
            data_ = { x:float(data[x])/float(maxval)*(self._height-1)
                      for x in data }
        else:
            data_ = { x:data[x] for x in data }
        return data_

    def _int(self,i):
        """
        Internal: convert value to the nearest integer
        """
        # Convert to nearest integer
        ii = floor(float(i))
        if (float(i) - float(ii)) < 0.5:
            return int(ii)
        else:
            return int(ii) + 1

    def _img(self):
        """
        Internal: return a Pillow 'Image' instance for the plot
        """
        img = Image.new('RGB',(self._width,self._height),
                        self._bgcolor)
        pixels = img.load()
        for x in range(self._width):
            for y in range(self._height):
                try:
                    pixels[x,self._height-y-1] = self._pixels[x][y]
                except KeyError:
                    pass
        return img

    def save(self,f,ext=".png"):
        """
        Save a copy of the plot to file

        Arguments:
          f (str): path to output file
          ext (str): optional, specify the plot extension
            (defaults to '.png')
        """
        return make_plot(self._img(),
                         outfile=f,
                         ext=ext)

    def encoded_png(self):
        """
        Return plot PNG as base64 encoded string
        """
        return make_plot(self._img(),inline=True)

def encode_png(png_file):
    """
    Return Base64 encoded string for a PNG
    """
    return "data:image/png;base64,%s" % \
        PNGBase64Encoder().encodePNG(png_file)
    
def uscreenplot(screen_files,outfile=None,screen_width=None,
                inline=None,use_legacy_colours=False):
    """
    Generate 'micro-plot' of FastqScreen outputs

    Arguments:
      screen_files (list): list of paths to one or more
        ...screen.txt files from FastqScreen
      outfile (str): path to output file
      screen_width (int): optional, set the width for
        each screen plot
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
      use_legacy_colours (boolean): if True then use the
        original colour palette from FastqScreen
        (default: False, use mimick the current colour
        palette)
    """
    # Mappings
    mappings = ('%One_hit_one_library',
                '%Multiple_hits_one_library',
                '%One_hit_multiple_libraries',
                '%Multiple_hits_multiple_libraries',)
    # Colours
    if use_legacy_colours:
        colors = (RGB_COLORS['blue'],
                  RGB_COLORS['navyblue'],
                  RGB_COLORS['red'],
                  RGB_COLORS['maroon'])
    else:
        colors = ((146,197,222),
                  (5,113,176),
                  (244,165,130),
                  (202,0,32))
    # Width for each screen plot
    if screen_width is None:
        screen_width = 50
    # Read in the screen data
    screens = []
    for screen_file in screen_files:
        screens.append(Fastqscreen(screen_file))
    nscreens = len(screens)
    # Make a small stacked bar chart
    bbox_color = (145,145,145)
    barwidth = 4
    width = nscreens*screen_width
    n_libraries_max = max([len(s) for s in screens])
    height = (n_libraries_max + 1)*(barwidth + 1)
    img = Image.new('RGB',(width,height),"white")
    pixels = img.load()
    # Process each screen in turn
    for nscreen,screen in enumerate(screens):
        xorigin = nscreen*screen_width
        xend = xorigin+screen_width-1
        yend = height-1
        # Draw a box around the plot
        for i in range(xorigin,xorigin+screen_width):
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
                    npx = int(ceil(data[mapping]/100.0*screen_width))
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
        npx = int(screen.no_hits/100.0*screen_width)
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

def uboxplot(fastqc_data=None,fastq=None,max_width=None,
             outfile=None,inline=None):
    """
    Generate FASTQ per-base quality 'micro-boxplot'

    'Micro-boxplot' is a thumbnail version of the per-base
    quality boxplots for a FASTQ file.

    Arguments:
      fastqc_data (str): path to a ``fastqc_data.txt``
        file
      fastq (str): path to a FASTQ file (quality stats
        will be extracted directly if ``fastqc_data``
        is not supplied)
      max_width (int): maximum width of plot in pixels;
        if ``None`` then by default width will be the
        number of bases, otherwise plots that would exceed
        this width will be scaled to fit
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
    # Set up data dictionaries for plotting
    p10 = { i:fastq_stats.p10[i] for i in range(fastq_stats.nbases) }
    p90 = { i:fastq_stats.p90[i] for i in range(fastq_stats.nbases) }
    q25 = { i:fastq_stats.q25[i] for i in range(fastq_stats.nbases) }
    q75 = { i:fastq_stats.q75[i] for i in range(fastq_stats.nbases) }
    median = { i:fastq_stats.median[i] for i in range(fastq_stats.nbases) }
    mean = { i:fastq_stats.mean[i] for i in range(fastq_stats.nbases) }
    # Initialise plot
    height = max_qual + 1
    width = fastq_stats.nbases
    if max_width:
        width = min(width,max_width)
    boxplot = Plot(width,height)
    # Create colour bands for different quality ranges
    boxplot.block((0,0),
                  (fastq_stats.nbases,20),
                  (230,175,175),RGB_COLORS.white)
    boxplot.block((0,20),
                  (fastq_stats.nbases,30),
                  (230,215,175),RGB_COLORS.white)
    boxplot.block((0,30),
                  (fastq_stats.nbases,max_qual),
                  (175,230,175),RGB_COLORS.white)
    # Add bounding box
    boxplot.bbox(RGB_COLORS.grey)
    # 10th-90th percentile (grey)
    boxplot.plot_range(p10,p90,RGB_COLORS.grey)
    # Interquartile range (yellow)
    boxplot.plot_range(q25,q75,RGB_COLORS.darkyellow1)
    # Plot median (red)
    boxplot.plot(median,RGB_COLORS.red)
    # Plot mean (blue)
    boxplot.plot(mean,RGB_COLORS.blue)
    # Output the plot image
    if inline:
        return boxplot.encoded_png()
    else:
        return boxplot.save(outfile)

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
        height = 50
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
    # Stripe the background
    for i in range(0,width,2):
        for j in range(0,height):
            pixels[i,j] = RGB_COLORS['white']
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
                   outfile=None,inline=False,width=50,height=12,
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

def uduplicationplot(total_deduplicated_percentage,height=None,
                     width=None,mode='dup',style='fancy',
                     warn_cutoff=None,fail_cutoff=None,
                     outfile=None,inline=False):
    """
    Make a 'micro' plot summarising sequence duplication

    Given the percentage of reads after deduplication (as
    calculated by FastQC's "Sequence Duplication Levels"
    module), plots a horizontal bar indicating the level
    of sequence (de)duplication.

    Two modes are available:

    - 'dup' shows the fraction of sequences removed after
      deduplication
    - 'dedup' shows the fraction of sequences remaining after
      deduplication

    Two styles are supported:

    - 'fancy' produces a multi-colour plot where the fraction
      of unique sequences is coloured according to pass,
      warn or fail cut-off levels, with the background to the
      bar striped with colours according to the cut-offs
    - 'simple' produces a two-colour plot where the fraction
      of unique sequences is shown in blue and the remainder
      in red

    Arguments:
      total_deduplicated_percentage (float): percentage of
        sequences remaining after deduplication
      height (int): height of the plot in pixels
      width (int): width of the plot in pixels
      style (str): either 'fancy' (default) or 'simple'
      mode (str): either 'dup' (default) or 'dedup'
      warn_cutoff (float): fraction of unique sequences
        below which the plot should indicate a warning
      fail_cutoff (float): fraction of unique sequences
        below which the plot should indicate a failure
      outfile (str): path for the output PNG
      inline (boolean): if True then returns the PNG
        as base64 encoded string rather than as a file
    """
    # Set mode for plot
    if mode == 'dedup':
        show_dedup = True
    elif mode == 'dup':
        show_dedup = False
    else:
        raise Exception("uduplicationplot: unknown mode '%s'"
                        % mode)
    # Set cutoffs for plot
    if warn_cutoff is None:
        # Set to 30%
        warn_cutoff = SEQUENCE_DEDUP_CUTOFFS['warn']
    if fail_cutoff is None:
        # Set to 20%
        fail_cutoff = SEQUENCE_DEDUP_CUTOFFS['fail']
    # Set plot limits
    if width is None:
        width = 50
    if height is None:
        height = 12
    # Set parameters according to style
    if style == 'simple':
        stripe_bg = False
        use_cutoffs = False
    elif style == 'fancy':
        stripe_bg = True
        use_cutoffs = True
    else:
        raise Exception("uduplicationplot: unknown style '%s'"
                        % style)
    # Set colours
    bg_color = "white"
    fg_color = "blue"
    fg_color_warn = "orange"
    fg_color_fail = "red"
    fg_color_fill = "red"
    if use_cutoffs:
        # Set colour based on cutoffs
        if total_deduplicated_percentage < fail_cutoff*100.0:
            dedup_color = fg_color_fail
        elif total_deduplicated_percentage < warn_cutoff*100.0:
            dedup_color = fg_color_warn
        else:
            dedup_color = fg_color
    else:
        # Ignore cutoffs
        dedup_color = fg_color
    # Scale cutoffs to plot width
    warn_cutoff = int(warn_cutoff*float(width))
    fail_cutoff = int(fail_cutoff*float(width))
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS[bg_color])
    pixels = img.load()
    # Create background striping
    if stripe_bg:
        for i in range(0,width,2):
            if i < fail_cutoff:
                stripe_color = (230,175,175) # "red"
            elif i < warn_cutoff:
                stripe_color = (230,215,175) # "orange"
            else:
                stripe_color = (175,225,255) # "blue"
            for j in range(1,height-1):
                if show_dedup:
                    ii = i
                else:
                    ii = width - i - 1
                pixels[ii,j] = stripe_color
    # Plot fraction of deduplicated sequences
    frac_dedup_seqs = float(total_deduplicated_percentage)/100.0
    frac_dup_seqs = 1.0 - frac_dedup_seqs
    start = 0
    if frac_dedup_seqs > 0:
        if show_dedup:
            end = max(1,int(frac_dedup_seqs*width))
        else:
            end = max(1,int(frac_dup_seqs*width))
        for i in range(start,end):
            for j in range(start,end):
                for j in range(1,height-1):
                    pixels[i,j] = RGB_COLORS[dedup_color]
        start = end
    else:
        end = start
    # Fill remainder of the plot with solid colour
    if not stripe_bg:
        for i in range(end,width):
            for j in range(1,height-1):
                pixels[i,j] = RGB_COLORS[fg_color_fail]
    # Output the plot
    return make_plot(img,
                     outfile=outfile,
                     inline=inline,
                     ext=".uduplication.png")

def uadapterplot(adapter_content,adapter_names=None,outfile=None,
                 inline=False,height=40,bar_width=10,
                 use_legacy_colours=False):
    """
    Make a 'micro' plot summarising adapter content

    The plot consists of a stacked vertical bar, with each
    stacked element indicating the relative portion of each
    adapter class in the sequence data.

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
      use_legacy_colours (boolean): if True then use the
        original colour palette from FastQC adapter plot
        (default: False, use mimick the current colour
        palette)
    """
    # Width of plot required for each bar
    width = bar_width + 4
    # Colours for each adapter
    if use_legacy_colours:
        fg_colors = (RGB_COLORS['red'],
                     RGB_COLORS['blue'],
                     RGB_COLORS['green'],
                     RGB_COLORS['black'])
    else:
        # Based on Tol colour scheme from
        # https://davidmathlogic.com/colorblind/
        fg_colors = ((136,34,85),
                     (51,34,136),
                     (17,119,51),
                     (221,204,119),
                     (68,170,153),
                     (170,68,153),
                     (204,102,119),
                     (136,204,238))
    # Number of colours
    ncolors = len(fg_colors)
    # Get adapter names
    if not adapter_names:
        adapter_names = sorted(adapter_content.keys())
    # Create the image
    img = Image.new('RGB',(width,height),RGB_COLORS['white'])
    pixels = img.load()
    # Stripe the background
    for i in range(1,height,2):
        for j in range(1,width-2):
            pixels[j,i] = RGB_COLORS['lightgrey']
    # Generate the plot
    start = 1
    for ii,adapter in enumerate(adapter_names):
        # Set colour based on adapter content
        fg_color = fg_colors[ii%ncolors]
        # Length of bar represents adapter content
        bar_length = int(adapter_content[adapter]*(height-2))
        # Draw the coloured part of the bar
        end = start + bar_length
        for i in range(start,end):
            for j in range(1,width-2):
                pixels[j,height-i] = fg_color
        start = end
    # Output the plot
    return make_plot(img,
                     outfile=outfile,
                     inline=inline,
                     ext=".uadapters.png")

def uinsertsizeplot(data,outfile=None,inline=False):
    """
    Return a mini-plot with the Picard insert size histogram

    Arguments:
      data (dict): dictionary mapping insert sizes to
        associated number of alignments (from Picard
        CollectInsertSizeMetrics)
      outfile (str): path for the output PNG
      inline (bool): if True then returns the PNG
        as base64 encoded string rather than as a file
    """
    p = Plot(50,40)
    p.stripe(RGB_COLORS.lightgrey,RGB_COLORS.white)
    p.plot(data,RGB_COLORS.red,fill=True)
    # Output the plot
    if inline:
        return p.encoded_png()
    else:
        return p.save(outfile,
                      ext=".insertsize.png")

def ucoverageprofileplot(data,outfile=None,inline=False):
    """
    Return a mini-plot of the Qualimap gene body coverage profile

    Arguments:
      data (dict): dictionary mapping transcript positions
        (percentile) to associated mean coverage depth (from
        Qualimap's RNA-seq analysis)
      outfile (str): path for the output PNG
      inline (bool): if True then returns the PNG
        as base64 encoded string rather than as a file
    """
    p = Plot(50,40)
    p.stripe(RGB_COLORS.gainsboro,RGB_COLORS.white)
    p.plot(data,RGB_COLORS.red,interpolation="minmax")
    # Output the plot
    if inline:
        return p.encoded_png()
    else:
        return p.save(outfile,
                      ext=".coverageprofile.png")

def ugenomicoriginplot(data,width=100,height=40,outfile=None,
                       inline=False):
    """
    Return a mini barplot of the Qualimap genomic origin of reads

    Arguments:
      data (dict): dictionary mapping genomic origin names
        to the associated percentage of reads (from Qualimap
        rnaseq 'Genomic Origin of Reads')
      height (int): plot height in pixels
      width (int): plot width in pixels
      outfile (str): path for the output PNG
      inline (bool): if True then returns the PNG
        as base64 encoded string rather than as a file
    """
    # Internal constants
    names = ['exonic',
             'intronic',
             'intergenic',
             'overlapping exon',
             'rRNA']
    # Colours taken from 'IBM' palette
    # See e.g. https://davidmathlogic.com/colorblind/
    colors = {
        'exonic': (100,143,255),
        'intronic': (120,94,240),
        'intergenic': (254,97,0),
        'overlapping exon': (220,38,127),
        'rRNA': (255,176,0)
    }
    # Filter list of names to just those in the data
    names = [n for n in names if n in data]
    # Set up barplot
    bar_length = width - 4
    bar_height = height - 10
    # Determine limits of the bar
    x1 = (width - bar_length)/2
    y1 = (height - bar_height)/2
    x2 = (width - bar_length)/2 + bar_length
    y2 = (height - bar_height)/2 + bar_height
    # Create plot
    p = Plot(width,height)
    p.bar([data[name][1] for name in names],
          (x1,y1),
          (x2,y2),
          [colors[name] for name in names])
    # Output the plot
    if inline:
        return p.encoded_png()
    else:
        return p.save(outfile,
                      ext=".genomicorigin.png")

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
