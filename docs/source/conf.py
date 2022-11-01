# -*- coding: utf-8 -*-
#
# auto_process_ngs documentation build configuration file, created by
# sphinx-quickstart on Fri Feb 27 10:23:25 2015.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os
sys.path.append('../')

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.insert(0, os.path.abspath('.'))

# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'auto_process_ngs'
copyright = u'2015-2022 University of Manchester'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), '..')))
from auto_process_ngs import get_version
#
# The short X.Y version.
version = get_version()
# The full version, including alpha/beta/rc tags.
release = get_version()

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# See https://stackoverflow.com/questions/18969093/how-to-include-the-toctree-in-the-sidebar-of-each-page
html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'auto_process_ngsdoc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'auto_process_ngs.tex', u'auto\\_process\\_ngs Documentation',
   u'Peter Briggs', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'auto_process_ngs', u'auto_process_ngs Documentation',
     [u'Peter Briggs'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'auto_process_ngs', u'auto_process_ngs Documentation',
   u'Peter Briggs', 'auto_process_ngs', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# -- Make example plot images for QC report -----------------------------------

from auto_process_ngs.qc import plots
imagesdir = os.path.join("images","auto","qc")
if not os.path.exists(imagesdir):
    os.makedirs(imagesdir)

# Read count plot examples
img = os.path.join(imagesdir,"read_count_uplot.png")
plots.ureadcountplot(3045784,nmasked=0,npadded=0,
                     max_reads=3665120,
                     outfile=img)

img = os.path.join(imagesdir,"read_count_uplot_small.png")
plots.ureadcountplot(1360307,nmasked=0,npadded=0,
                     max_reads=3665120,
                     outfile=img)

img = os.path.join(imagesdir,"read_count_uplot_masking_and_padding1.png")
plots.ureadcountplot(58192145,nmasked=1087126,npadded=1455460,
                     max_reads=58192145,
                     outfile=img)

img = os.path.join(imagesdir,"read_count_uplot_masking_and_padding2.png")
plots.ureadcountplot(50181731,nmasked=21227806,npadded=28246553,
                     max_reads=58192145,
                     outfile=img)

# Sequence length distribution plot examples
img = os.path.join(imagesdir,"seq_dist_uplot.png")
plots.useqlenplot({ 65: 3.0,
                    66: 6.0,
                    67: 2.0,
                    68: 1.0,
                    69: 16.0,
                    70: 218.0,
                    71: 307.0,
                    72: 1289.0,
                    73: 1559577.0,
                    74: 5165990.0,
                    75: 10549970,
                    76: 20008064 },
                  min_len=26,max_len=76,outfile=img)

img = os.path.join(imagesdir,"seq_dist_uplot_slewed.png")
plots.useqlenplot({ 10: 2017916,
                    11: 167219,
                    12: 223764,
                    13: 387221,
                    14: 244739,
                    15: 247888,
                    16: 272688,
                    17: 478475,
                    18: 543203,
                    19: 805492,
                    20: 2681762,
                    21: 7473567,
                    22: 12133952,
                    23: 5344314,
                    24: 2067501,
                    25: 342139,
                    26: 224964,
                    27: 349178,
                    28: 319042,
                    29: 332929,
                    30: 290437,
                    31: 401247,
                    32: 497299,
                    33: 327686,
                    34: 388090,
                    35: 115688,
                    36: 75495,
                    37: 46890,
                    38: 14520,
                    39: 11042,
                    40: 3499,
                    41: 2577,
                    42: 3838,
                    43: 9538,
                    44: 14095,
                    45: 14476,
                    46: 852,
                    47: 445,
                    48: 307,
                    49: 378,
                    50: 325,
                    51: 398,
                    52: 483,
                    53: 431,
                    54: 252,
                    55: 405,
                    56: 285,
                    57: 336,
                    58: 709,
                    59: 773,
                    60: 1697,
                    61: 12008,
                    62: 8188,
                    63: 9502,
                    64: 2204,
                    65: 4104,
                    66: 7917,
                    67: 9872,
                    68: 13565,
                    69: 8609,
                    70: 1932,
                    71: 29161,
                    72: 7255,
                    73: 4784,
                    74: 6592,
                    75: 127139,
                    76: 75135 },
                  min_len=10,max_len=76,outfile=img)

img = os.path.join(imagesdir,"seq_dist_uplot_untrimmed.png")
plots.useqlenplot({ 98: 3045784 },min_len=48,max_len=98,outfile=img)

# Sequence duplication plot examples
img = os.path.join(imagesdir,"duplication_uplot_pass.png")
plots.uduplicationplot(45.08,outfile=img)

img = os.path.join(imagesdir,"duplication_uplot_warn.png")
plots.uduplicationplot(25.73,outfile=img)

img = os.path.join(imagesdir,"duplication_uplot_fail.png")
plots.uduplicationplot(9.63,outfile=img)

img = os.path.join(imagesdir,"duplication_uplot_bg.png")
plots.uduplicationplot(0.0,outfile=img)

# Adapter content plot examples
img = os.path.join(imagesdir,"adapter_uplot_no_adptrs.png")
plots.uadapterplot({'Illumina Universal Adapter': 0.0,
                    'Illumina Small RNA Adapter': 0.0,
                    'Nextera Transposase Sequence': 0.0,
                    'SOLID Small RNA Adapter': 0.0 },
                   ('Illumina Universal Adapter',
                    'Illumina Small RNA Adapter',
                    'Nextera Transposase Sequence',
                    'SOLID Small RNA Adapter'),outfile=img)

img = os.path.join(imagesdir,"adapter_uplot_adptrs_sml.png")
plots.uadapterplot({'Illumina Universal Adapter': 0.0,
                    'Illumina Small RNA Adapter': 0.0,
                    'Nextera Transposase Sequence': 0.2,
                    'SOLID Small RNA Adapter': 0.0 },
                   ('Illumina Universal Adapter',
                    'Illumina Small RNA Adapter',
                    'Nextera Transposase Sequence',
                    'SOLID Small RNA Adapter'),outfile=img)

img = os.path.join(imagesdir,"adapter_uplot_adptrs_lrg.png")
plots.uadapterplot({'Illumina Universal Adapter': 0.57,
                    'Illumina Small RNA Adapter': 0.0,
                    'Nextera Transposase Sequence': 0.0,
                    'SOLID Small RNA Adapter': 0.0 },
                   ('Illumina Universal Adapter',
                    'Illumina Small RNA Adapter',
                    'Nextera Transposase Sequence',
                    'SOLID Small RNA Adapter'),outfile=img)

# Picard insert sizes
img = os.path.join(imagesdir,"picard_insert_size_uplot.png")
plots.uinsertsizeplot({ 25: 1, 27: 1, 28: 3, 29: 1, 30: 3, 31: 2, 32: 4,
                        33: 3, 34: 2, 35: 4, 36: 8, 37: 5, 38: 3, 39: 7,
                        40: 5, 41: 4, 42: 5, 43: 9, 44: 2, 45: 14, 46: 4,
                        47: 10, 48: 10, 49: 13, 50: 12, 51: 13, 52: 13,
                        53: 15, 54: 14, 55: 20, 56: 24, 57: 24, 58: 18,
                        59: 20, 60: 29, 61: 31, 62: 39, 63: 46, 64: 38,
                        65: 52, 66: 66, 67: 53, 68: 67, 69: 76, 70: 71,
                        71: 89, 72: 93, 73: 79, 74: 104, 75: 119, 76: 127,
                        77: 137, 78: 170, 79: 212, 80: 161, 81: 226,
                        82: 221, 83: 243, 84: 263, 85: 285, 86: 284,
                        87: 297, 88: 334, 89: 348, 90: 362, 91: 438,
                        92: 432, 93: 475, 94: 522, 95: 555, 96: 555,
                        97: 659, 98: 599, 99: 665, 100: 644, 101: 710,
                        102: 726, 103: 685, 104: 786, 105: 776, 106: 787,
                        107: 921, 108: 971, 109: 978, 110: 994, 111: 920,
                        112: 1052, 113: 959, 114: 963, 115: 990, 116: 989,
                        117: 1072, 118: 1205, 119: 1207, 120: 1205,
                        121: 1256, 122: 1086, 123: 1081, 124: 948,
                        125: 1102, 126: 1061, 127: 1212, 128: 1238,
                        129: 1112, 130: 1116, 131: 911, 132: 904, 133: 773,
                        134: 896, 135: 956, 136: 813, 137: 847, 138: 768,
                        139: 794, 140: 798, 141: 814, 142: 696, 143: 817,
                        144: 741, 145: 732, 146: 720, 147: 767, 148: 667,
                        149: 739, 150: 695, 151: 749, 152: 737, 153: 623,
                        154: 647, 155: 660, 156: 598, 157: 609, 158: 528,
                        159: 568, 160: 533, 161: 588, 162: 518, 163: 550,
                        164: 513, 165: 502, 166: 478, 167: 484, 168: 502,
                        169: 506, 170: 562, 171: 629, 172: 544, 173: 443,
                        174: 500, 175: 435, 176: 390, 177: 399, 178: 430,
                        179: 375, 180: 333, 181: 348, 182: 343, 183: 315,
                        184: 318, 185: 298, 186: 339, 187: 300, 188: 297,
                        189: 326, 190: 334, 191: 317, 192: 296, 193: 285,
                        194: 240, 195: 260, 196: 258, 197: 216, 198: 294,
                        199: 284, 200: 206, 201: 236, 202: 214, 203: 183,
                        204: 209, 205: 180, 206: 224, 207: 221, 208: 198,
                        209: 184, 210: 178, 211: 195, 212: 177, 213: 191,
                        214: 172, 215: 165, 216: 159, 217: 156, 218: 148,
                        219: 141, 220: 128, 221: 167, 222: 133, 223: 165,
                        224: 149, 225: 142, 226: 158, 227: 129, 228: 123,
                        229: 132, 230: 122, 231: 110, 232: 112, 233: 137,
                        234: 97, 235: 104, 236: 103, 237: 108, 238: 93,
                        239: 103, 240: 92, 241: 82, 242: 110, 243: 81,
                        244: 101, 245: 79, 246: 105, 247: 83, 248: 92,
                        249: 84, 250: 75, 251: 59, 252: 76, 253: 61,
                        254: 73, 255: 57, 256: 84, 257: 86, 258: 71,
                        259: 69, 260: 74, 261: 70, 262: 67, 263: 64,
                        264: 66, 265: 73, 266: 74, 267: 59, 268: 70,
                        269: 48, 270: 67, 271: 56, 272: 51, 273: 59,
                        274: 73, 275: 60, 276: 62, 277: 56, 278: 46,
                        279: 48, 280: 54, 281: 43, 282: 50, 283: 55,
                        284: 40, 285: 45, 286: 28, 287: 49, 288: 39,
                        289: 40, 290: 45, 291: 46, 292: 43, 293: 44,
                        294: 49, 295: 32, 296: 46, 297: 40, 298: 29,
                        299: 30, 300: 43, 301: 31, 302: 27, 303: 43,
                        304: 32, 305: 31, 306: 31, 307: 28, 308: 35,
                        309: 33, 310: 34, 311: 33, 312: 29, 313: 27,
                        314: 35, 315: 29, 316: 31, 317: 23, 318: 30,
                        319: 22, 320: 27, 321: 33, 322: 22, 323: 24,
                        324: 20, 325: 27, 326: 17, 327: 22, 328: 22,
                        329: 28, 330: 21, 331: 31, 332: 20, 333: 34,
                        334: 24, 335: 20, 336: 12, 337: 17, 338: 19,
                        339: 27, 340: 13, 341: 25, 342: 12, 343: 11,
                        344: 16, 345: 11, 346: 15, 347: 18, 348: 15,
                        349: 14, 350: 17, 351: 18, 352: 13, 353: 22,
                        354: 16, 355: 12, 356: 11, 357: 16, 358: 17,
                        359: 13, 360: 15, 361: 11, 362: 15, 363: 20,
                        364: 17, 365: 15, 366: 13, 367: 8, 368: 6,
                        369: 11, 370: 10, 371: 13, 372: 9, 373: 6,
                        374: 10, 375: 20, 376: 10, 377: 17, 378: 10,
                        379: 7, 380: 13, 381: 15, 382: 10, 383: 12,
                        384: 10, 385: 13, 386: 8, 387: 8, 388: 9, 389: 9,
                        390: 10, 391: 12, 392: 6, 393: 9, 394: 8, 395: 4,
                        396: 13, 397: 8, 398: 7, 399: 10, 400: 7, 401: 5,
                        402: 14, 403: 8, 404: 10, 405: 2, 406: 12, 407: 5,
                        408: 7, 409: 8, 410: 15, 411: 8, 412: 5, 413: 8,
                        414: 9, 415: 13, 416: 10, 417: 5, 418: 8, 419: 7,
                        420: 9, 421: 8, 422: 9, 423: 9, 424: 7, 425: 6,
                        426: 6, 427: 3, 428: 6, 429: 3, 430: 8, 431: 6,
                        432: 6, 433: 2, 434: 3, 435: 4, 436: 3, 437: 9,
                        438: 9, 439: 8, 440: 7, },outfile=img)

# Qualimap coverage
img = os.path.join(imagesdir,"qualimap_gene_body_coverage_uplot.png")
plots.ucoverageprofileplot({0: 69.0, 1: 82.6, 2: 98.7, 3: 113.9,
                            4: 129.5, 5: 137.5, 6: 147.9, 7: 157.6,
                            8: 163.3, 9: 163.3, 10: 174.2, 11: 184.6,
                            12: 186.9, 13: 188.3, 14: 185.1, 15: 171.5,
                            16: 168.9, 17: 166.3, 18: 164.7, 19: 169.9,
                            20: 172.2, 21: 170.6, 22: 172.7, 23: 182.1,
                            24: 188.9, 25: 191.5, 26: 199.4, 27: 203.1,
                            28: 200.6, 29: 197.7, 30: 205.4, 31: 197.0,
                            32: 195.3, 33: 200.7, 34: 189.9, 35: 186.5,
                            36: 187.3, 37: 179.8, 38: 174.1, 39: 178.7,
                            40: 184.3, 41: 188.1, 42: 187.1, 43: 189.5,
                            44: 181.2, 45: 171.5, 46: 167.8, 47: 165.0,
                            48: 162.8, 49: 162.5, 50: 161.3, 51: 159.4,
                            52: 159.1, 53: 164.8, 54: 173.8, 55: 171.2,
                            56: 167.4, 57: 169.4, 58: 168.2, 59: 166.7,
                            60: 164.9, 61: 165.2, 62: 163.3, 63: 164.8,
                            64: 165.8, 65: 160.9, 66: 156.1, 67: 154.3,
                            68: 153.4, 69: 153.7, 70: 151.0, 71: 142.1,
                            72: 131.0, 73: 129.4, 74: 128.0, 75: 128.5,
                            76: 125.8, 77: 124.7, 78: 122.7, 79: 119.4,
                            80: 118.7, 81: 117.4, 82: 116.7, 83: 116.5,
                            84: 113.9, 85: 113.3, 86: 117.4, 87: 119.8,
                            88: 120.1, 89: 113.7, 90: 114.2, 91: 122.7,
                            92: 124.3, 93: 115.0, 94: 102.9, 95: 88.8,
                            96: 55.9, 97: 37.7, 98: 23.3, 99: 7.5},
                           outfile=img)

# Qualimap genomic origin of reads
img = os.path.join(imagesdir,"qualimap_genomic_origin_reads.png")
plots.ugenomicoriginplot({'exonic': (83408, 79.74),
                          'intronic': (18374, 17.57),
                          'intergenic': (2818, 2.69),
                          'overlapping exon': (1767, 1.69),
                          'rRNA': (0, 0.0)},
                         outfile=img)

# -- Make command reference documents ------------------------------------------

if not os.path.exists("reference"):
    os.mkdir("reference")

commandref = os.path.join(os.getcwd(),"reference","commands.rst")
with open(commandref,'w') as commands:
    commands.write("""
``auto_process`` commands
=========================

.. note::

   This documentation has been auto-generated from the
   command help

``auto_process.py`` implements the following commands:

.. contents:: :local:

""")

import subprocess
for subcmd in ("setup",
               "make_fastqs",
               "analyse_barcodes",
               "setup_analysis_dirs",
               "run_qc",
               "publish_qc",
               "archive",
               "report",
               "samplesheet",
               "merge_fastq_dirs",
               "update_fastq_stats",
               "import_project",
               "config",
               "params",
               "metadata",
               "readme",
               "clone"):
    # Capture the output
    help_text_file = "%s.help" % subcmd
    with open(help_text_file,'w') as fp:
        subprocess.call(['auto_process.py',subcmd,'--help'],stdout=fp)
    # Write into the document
    with open(commandref,'a') as fp:
        help_text = open(help_text_file,'r').read()
        title = "%s" % subcmd
        ref = ".. _commands_%s:" % subcmd
        fp.write("%s\n\n%s\n%s\n\n::\n\n" % (ref,
                                             title,
                                             "*"*len(title)))
        for line in help_text.split('\n'):
            fp.write("    %s\n" % line)
        os.remove(help_text_file)

# -- Make utility reference documents ------------------------------------------

utilityref = os.path.join(os.getcwd(),"reference","utilities.rst")
with open(utilityref,'w') as utilities:
    utilities.write("""
Utilities
=========

.. note::

   This documentation has been auto-generated from the
   command help

In addition to the main ``auto_process.py`` command, a number of utilities
are available:

.. contents:: :local:

""")
for utility in ("analyse_barcodes.py",
                "assign_barcodes.py",
                "audit_projects.py",
                "build_index.py",
                "concat_fastqs.py",
                "barcode_splitter.py",
                "download_fastqs.py",
                "demultiplex_icell8_atac.py",
                "fastq_statistics.py",
                "icell8_contamination_filter.py",
                "icell8_report.py",
                "icell8_stats.py",
                "manage_fastqs.py",
                "process_icell8.py",
                "run_qc.py",
                "split_icell8_fastqs.py",
                "transfer_data.py",
                "update_project_metadata.py"):
    # Capture the output
    help_text_file = "%s.help" % utility
    with open(help_text_file,'w') as fp:
        subprocess.call([utility,'--help'],
                        stdout=fp,
                        stderr=subprocess.STDOUT)
    # Write into the document
    with open(utilityref,'a') as fp:
        help_text = open(help_text_file,'r').read()
        title = "%s" % utility
        ref = ".. _utilities_%s:" % os.path.splitext(utility)[0]
        fp.write("%s\n\n%s\n%s\n\n::\n\n" % (ref,
                                             title,
                                             "*"*len(title)))
        if not help_text:
            help_text = "No output from --help command?"
        for line in help_text.split('\n'):
            fp.write("    %s\n" % line)
        os.remove(help_text_file)

# -- Make developers reference documents ---------------------------------------

# Fetch a list of modules
# See https://stackoverflow.com/a/1708706/579925
from pkgutil import walk_packages
def get_modules(pkg):
    modlist = []
    for importer,modname,ispkg in walk_packages(path=pkg.__path__,
                                                prefix=pkg.__name__+'.',
                                                onerror=lambda x: None):
        modname = '.'.join(modname.split('.')[1:])
        modlist.append(modname)
    return modlist

import auto_process_ngs
modlist = get_modules(auto_process_ngs)

# Ensure the 'developers' subdir exists
devdocdir = os.path.join(os.getcwd(),"developers")
if not os.path.exists(devdocdir):
    print("Making %s" % devdocdir)
    os.mkdir(devdocdir)

# Ensure the 'api_docs' subdir exists
api_doc_dir = os.path.join(os.getcwd(),
                           "developers",
                           "api_docs")
if not os.path.exists(api_doc_dir):
    print("Making %s" % api_doc_dir)
    os.mkdir(api_doc_dir)

# Generate documents for each module
for modname in modlist:
    docname = modname.replace('.','_')
    docfile = os.path.join(api_doc_dir,
                           "%s.rst" % docname)
    print("Generating %s" % docfile)
    with open(docfile,'w') as doc:
         title = "``auto_process_ngs.%s``" % modname
         doc.write("""%s
%s

.. automodule:: auto_process_ngs.%s
   :members:
""" % (title,'='*len(title),modname))

# Generate an index
api_index = os.path.join(api_doc_dir,"index.rst")
print("Writing %s" % api_index)
with open(api_index,'w') as doc:
    doc.write("""=============================
Developers' API documentation
=============================

.. toctree::
   :maxdepth: 2

""")
    for modname in modlist:
        doc.write("   %s\n" % modname.replace('.','_'))


