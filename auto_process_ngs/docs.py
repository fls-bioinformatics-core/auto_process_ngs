#!/usr/bin/env python
#
#     docs: helpers for auto-process-ngs Sphinx documentation generation
#     Copyright (C) University of Manchester 2025-2026 Peter Briggs
#

"""
Provides utility classes and functions to help with generating Sphinx
documentation for the ``auto_process_ngs`` package.

The helpers are intended to be used within the Sphinx ``conf.py`` file,
to generate content such as tables of applications and protocols which
are then included in the documentation source files, and API and reference
documentation.

Classes:

* ``RstSimpleTable``: class for making reStructuredText 'simple' tables
* ``RstGridTable``: class for making reStructuredText 'grid' tables
* ``ExampleQCPlots``: class for making example QC micro-plots

Functions:

* ``generate_qc_protocols_tbl``: generate table of QC protocols
* ``generate_fq_protocols_tbl``: generate table of Fastq generation protocols
* ``generate_fq_prototocols_for_apps``: table of Fastq generation protocols
  for specified applications
* ``generate_applications_tbl``: generate table of applications
* ``generate_api_docs``: generate developer's API documentation
* ``generate_utility_docs``: generate documentation for utilities
* ``generate_command_docs``: generate documentation for 'auto_process' commands
* ``get_modules``: list Python modules in a package
* ``get_help_text``: fetch text from running program with '--help'

"""

import os
import subprocess
from .applications import fetch_application_data
from .bcl2fastq.protocols import PROTOCOLS as FQ_PROTOCOLS
from .qc.protocols import fetch_protocol_definition
from .qc import plots
from .qc.protocols import QC_PROTOCOLS
from pkgutil import walk_packages


class RstSimpleTable:
    """
    Class for making simple reStructuredText tables.

    See documentation for reStructureText simple tables at:
    https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#simple-tables

    Usage:

    >>> table_data = [
    ...     ['col1_row1','col2_row1','col3_row1'],
    ...     ['col1_row2','col2_row2','col3_row2'],
    ...     ['col1_row3','col2_row3','col3_row3'],
    ... ]
    >>> header = ['Column 1','Column 2','Column 3']
    >>> table = RstSimpleTable(table_data)
    >>> for line in table.construct_table(header=header):
    ...     print(line)
    =========  =========  =========
    Column 1   Column 2   Column 3
    =========  =========  =========
    col1_row1  col2_row1  col3_row1
    col1_row2  col2_row2  col3_row2
    col1_row3  col2_row3  col3_row3
    =========  =========  =========

    Arguments:
        table_data (list): list of table rows, where each
            row is a list of column values.
    """
    def __init__(self, table_data):
        self._table_data = table_data

    def get_field_widths(self, header=None):
        """
        Returns a list of field widths.
        """
        field_widths = []
        for row in self._table_data:
            for i, col in enumerate(row):
                try:
                    field_widths[i] = max(field_widths[i], len(col))
                except IndexError:
                    field_widths.append(len(col))
        if header:
            for i, title in enumerate(header):
                field_widths[i] = max(field_widths[i], len(title))
        return field_widths

    def make_divider(self, field_widths):
        """
        Make a row divider
        """
        divider = []
        for width in field_widths:
            divider.append("="*width)
        return "  ".join(divider)

    def construct_table(self, header=None, indent=""):
        """
        Returns reStructuredText simple table

        Arguments:
            header (list): list of column titles to use in
                the table header (otherwise no header is made)
            indent (str): string to use for indenting each
                line of the table (default: no indentation)
        """
        # Collect the field widths
        field_widths = self.get_field_widths(header=header)
        # Start constructing the table
        table = []
        # Add top divider
        table.append(indent + self.make_divider(field_widths))
        # Add table header
        if header:
            line = []
            for title, width in zip(header, field_widths):
                line.append(title + " "*(width - len(title)))
            table.append(indent + "  ".join(line))
            # Add header divider
            table.append(indent + self.make_divider(field_widths))
        # Add the table contents
        for row in self._table_data:
            line = []
            for col, width in zip(row, field_widths):
                line.append(col + " "*(width - len(col)))
            table.append(indent + "  ".join(line))
        # Closing divider
        table.append(indent + self.make_divider(field_widths))
        return table


class RstGridTable:
    """
    Class for making reStructuredText 'grid' tables.

    See documentation for reStructureText grid tables at:
    https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#grid-tables

    Usage:

    >>> table_data = [
    ...     ['col1_row1','col2_row1','col3_row1'],
    ...     ['col1_row2','col2_row2','col3_row2'],
    ...     ['col1_row3','col2_row3','col3_row3'],
    ... ]
    >>> header = ['Column 1','Column 2','Column 3']
    >>> table = RstGridTable(table_data)
    >>> for line in table.construct_table(header=header):
    ...     print(line)
    +-----------+-----------+-----------+
    | Column 1  | Column 2  | Column 3  |
    +===========+===========+===========+
    | col1_row1 | col2_row1 | col3_row1 |
    +-----------+-----------+-----------+
    | col1_row2 | col2_row2 | col3_row2 |
    +-----------+-----------+-----------+
    | col1_row3 | col2_row3 | col3_row3 |
    +-----------+-----------+-----------+

    Where columns with identical values in adjacent rows are merged:

    >>> table_data = [
    ...     ['block1','col2_row1','col3_row1'],
    ...     ['block1','col2_row2','col3_row2'],
    ...     ['block2','col2_row3','col3_row3'],
    ... ]
    >>> header = ['Block','Column 2','Column 3']
    >>> table = RstGridTable(table_data)
    >>> for line in table.construct_table(header=header):
    ...     print(line)
    +--------+-----------+-----------+
    | Block  | Column 2  | Column 3  |
    +========+===========+===========+
    | block1 | col2_row1 | col3_row1 |
    |        +-----------+-----------+
    |        | col2_row2 | col3_row2 |
    +--------+-----------+-----------+
    | block2 | col2_row3 | col3_row3 |
    +--------+-----------+-----------+

    Arguments:
        table_data (list): list of table rows, where each
            row is a list of column values.
    """
    def __init__(self, table_data):
        self._table_data = table_data

    def get_field_widths(self, header=None):
        """
        Returns a list of field widths.
        """
        field_widths = []
        for row in self._table_data:
            for i, col in enumerate(row):
                try:
                    field_widths[i] = max(field_widths[i], len(col))
                except IndexError:
                    field_widths.append(len(col))
        if header:
            for i, title in enumerate(header):
                field_widths[i] = max(field_widths[i], len(title))
        return field_widths

    def make_header(self, header, field_widths):
        """
        Make the table header
        """
        table_header = [self.make_divider(field_widths)]
        line = []
        for title, width in zip(header, field_widths):
            line.append(title + " "*(width - len(title)))
        table_header.append("| " + " | ".join(line) + " |")
        return table_header

    def make_divider(self, field_widths, divider_char="-"):
        """
        Make a row divider
        """
        divider = []
        for width in field_widths:
            divider.extend(["+", divider_char*(width + 2)])
        return "".join(divider) + "+"

    def construct_table(self, header=None, indent=""):
        """
        Returns reStructuredText table

        Arguments:
            header (list): list of column titles to use in
                the table header (otherwise no header is made)
            indent (str): string to use for indenting each
                line of the table (default: no indentation)
        """
        # Collect the field widths
        field_widths = self.get_field_widths(header=header)
        # Start constructing the table
        table = []
        # Add table header
        if header:
            for line in self.make_header(header, field_widths):
                table.append(indent + line)
        # Add the table contents
        previous_row = None
        for row in self._table_data:
            if previous_row is None:
                # First line
                previous_row = row
                # Add a divider
                if header:
                    # Separate from header
                    divider = self.make_divider(field_widths, divider_char="=")
                else:
                    # No header
                    divider = self.make_divider(field_widths)
            else:
                # Inside table body
                if previous_row[0] != row[0]:
                    # First column data differs from previous row,
                    # so treat as a new "block"
                    # Store this row for next round
                    previous_row = row
                    divider = self.make_divider(field_widths)
                else:
                    # Inside a "block": look at merging rows where values match
                    # within this column
                    new_row = []
                    for col, prev_col, width in zip(row, previous_row, field_widths):
                        if col != prev_col:
                            # Column contents differ between rows so
                            # insert explicit value
                            new_row.append(col)
                        else:
                            # Column contents are the same so insert
                            # 'None' to indicate merged cell
                            new_row.append(None)
                    previous_row = row
                    row = new_row
                    # Create divider
                    divider = ""
                    for i, col in enumerate(row):
                        if col or (i>0 and row[i-1]):
                            divider += "+"
                        else:
                            divider += "|"
                        if col:
                            divider += "-"*(field_widths[i]+2)
                        else:
                            divider += " "*(field_widths[i]+2)
                    if divider.endswith("-"):
                        divider += "+"
                    else:
                        divider += "|"
            # Append the row
            table.append(indent + divider)
            line = ""
            for i, col in enumerate(row):
                if col:
                    # Pad value
                    line += "| " + col + " "*(field_widths[i]-len(col)+1)
                else:
                    # Empty cell
                    line += "| " + " "*(field_widths[i]+1)
            line += "|"
            table.append(indent + line)
        # Closing divider
        table.append(indent + self.make_divider(field_widths))
        return table


class ExampleQCPlots:
    """
    Class for generating example QC plots.

    Provides a set of static methods for generating example
    versions of various QC "microplots" (read counts,
    sequence length distribution, sequence duplication,
    adapter content, Picard insert sizes, Qualimap gene body
    coverage and genomic origin),
    """

    @staticmethod
    def read_counts(images_dir):
        """
        Generate example read count plots
        """
        # Read count plot examples
        img = os.path.join(images_dir, "read_count_uplot.png")
        plots.ureadcountplot(3045784, nmasked=0, npadded=0,
                             max_reads=3665120,
                             outfile=img)
        img = os.path.join(images_dir, "read_count_uplot_small.png")
        plots.ureadcountplot(1360307, nmasked=0, npadded=0,
                             max_reads=3665120,
                             outfile=img)
        img = os.path.join(images_dir, "read_count_uplot_masking_and_padding1.png")
        plots.ureadcountplot(58192145, nmasked=1087126, npadded=1455460,
                             max_reads=58192145,
                             outfile=img)
        img = os.path.join(images_dir, "read_count_uplot_masking_and_padding2.png")
        plots.ureadcountplot(50181731, nmasked=21227806, npadded=28246553,
                             max_reads=58192145,
                             outfile=img)

    @staticmethod
    def sequence_lengths(images_dir):
        """
        Generate example sequence length distribution plots
        """
        img = os.path.join(images_dir, "seq_dist_uplot.png")
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
        img = os.path.join(images_dir, "seq_dist_uplot_slewed.png")
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
        img = os.path.join(images_dir,"seq_dist_uplot_untrimmed.png")
        plots.useqlenplot({ 98: 3045784 }, min_len=48,max_len=98,outfile=img)

    @staticmethod
    def sequence_duplication(images_dir):
        """
        Generate example sequence duplication plots
        """
        # Sequence duplication plot examples
        img = os.path.join(images_dir, "duplication_uplot_pass.png")
        plots.uduplicationplot(45.08, outfile=img)

        img = os.path.join(images_dir, "duplication_uplot_warn.png")
        plots.uduplicationplot(25.73, outfile=img)

        img = os.path.join(images_dir, "duplication_uplot_fail.png")
        plots.uduplicationplot(9.63, outfile=img)

        img = os.path.join(images_dir, "duplication_uplot_bg.png")
        plots.uduplicationplot(0.0, outfile=img)

    @staticmethod
    def adapter_content(images_dir):
        """
        Generate example adapter content plots
        """
        # No adapters
        img = os.path.join(images_dir, "adapter_uplot_no_adptrs.png")
        plots.uadapterplot({"Illumina Universal Adapter": 0.0,
                            "Illumina Small RNA 3' Adapter": 0.0,
                            "Illumina Small RNA 5' Adapter": 0.0,
                            "Nextera Transposase Sequence": 0.0,
                            "PolyA": 0.0,
                            "PolyG": 0.0},
                           ["Illumina Universal Adapter",
                            "Illumina Small RNA 3' Adapter",
                            "Illumina Small RNA 5' Adapter",
                            "Nextera Transposase Sequence",
                            "PolyA",
                            "PolyG"], outfile=img)
        # Small adapter content
        img = os.path.join(images_dir, "adapter_uplot_adptrs_sml.png")
        plots.uadapterplot({"Illumina Universal Adapter": 0.0,
                            "Illumina Small RNA 3' Adapter": 0.0,
                            "Illumina Small RNA 5' Adapter": 0.0,
                            "Nextera Transposase Sequence": 0.2,
                            "PolyA": 0.0,
                            "PolyG": 0.0},
                           ["Illumina Universal Adapter",
                            "Illumina Small RNA 3' Adapter",
                            "Illumina Small RNA 5' Adapter",
                            "Nextera Transposase Sequence",
                            "PolyA",
                            "PolyG"], outfile=img)
        # Large adapter content
        img = os.path.join(images_dir, "adapter_uplot_adptrs_lrg.png")
        plots.uadapterplot({"Illumina Universal Adapter": 0.57,
                            "Illumina Small RNA 3' Adapter": 0.0,
                            "Illumina Small RNA 5' Adapter": 0.0,
                            "Nextera Transposase Sequence": 0.0,
                            "PolyA": 0.0,
                            "PolyG": 0.0},
                           ["Illumina Universal Adapter",
                            "Illumina Small RNA 3' Adapter",
                            "Illumina Small RNA 5' Adapter",
                            "Nextera Transposase Sequence",
                            "PolyA",
                            "PolyG"], outfile=img)

    @staticmethod
    def picard_insert_sizes(images_dir):
        """
        Generate example Picard insert sizes plots
        """
        img = os.path.join(images_dir, "picard_insert_size_uplot.png")
        plots.uinsertsizeplot({25: 1, 27: 1, 28: 3, 29: 1, 30: 3, 31: 2, 32: 4,
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
                               438: 9, 439: 8, 440: 7, }, outfile=img)

    @staticmethod
    def qualimap(images_dir):
        """
        Generate example Qualimp plots
        """
        # Qualimap coverage
        img = os.path.join(images_dir, "qualimap_gene_body_coverage_uplot.png")
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
        img = os.path.join(images_dir, "qualimap_genomic_origin_reads.png")
        plots.ugenomicoriginplot({'exonic': (83408, 79.74),
                                  'intronic': (18374, 17.57),
                                  'intergenic': (2818, 2.69),
                                  'overlapping exon': (1767, 1.69),
                                  'rRNA': (0, 0.0)},
                                 outfile=img)


def generate_fq_protocols_tbl(docfile):
    """
    Generate table with the Fastq generation protocols

    Table lists the Fastq protocols and has columns
    with protocol name, description and read lengths.

    Arguments:
        docfile (str): path to output RST file
    """
    fq_protocols_data = []
    for p in FQ_PROTOCOLS:
        description = FQ_PROTOCOLS[p]["description"]
        reads = []
        for r in ("r1_length", "i1_length", "i2_length", "r2_length", "r3_length"):
            try:
                reads.append(str(FQ_PROTOCOLS[p][r]))
            except KeyError:
                pass
        fq_protocols_data.append([f"``{p}``", description, " | ".join(reads)])
    fq_protocols_data = sorted(fq_protocols_data, key=lambda p: p[0])
    tbl = RstSimpleTable(fq_protocols_data)
    with open(docfile, "wt") as fp:
        fp.write("\n".join(tbl.construct_table(
            header=["Protocol", "Description", "Read lengths"], )))


def generate_fq_protocols_for_apps(docfile, header=["Application", "Protocol"],
                                   tags=None):
    """
    Fastq generation protocols table for specific applications

    Table lists the Fastq protocols and has columns
    with application description and protocol name.

    Arguments:
        docfile (str): path to output RST file
        header (list): list of column titles to use in the table
        tags (list): list of tags to include in the table
    """
    if tags:
        # Only include protocols associated with applications
        # which match the supplied tags
        fq_protocols = set()
        for app in fetch_application_data(tags=tags, expand=True):
            fq_protocols.add(app["fastq_generation"])
    else:
        # All protocols
        fq_protocols = FQ_PROTOCOLS
    fq_protocols_data = []
    for fq_protocol in sorted(fq_protocols):
        fq_protocols_data.append([FQ_PROTOCOLS[fq_protocol]["description"],
                                  f"``{fq_protocol}``"])
    tbl = RstSimpleTable(fq_protocols_data)
    with open(docfile, "wt") as fp:
        fp.write("\n".join(tbl.construct_table(header=header)))


def generate_applications_table(docfile, tags=None):
    """
    Generates table of applications

    The table will only include applications which match
    the specified tags (if supplied), and will have the
    columns 'Platform', 'Library', and optionally a third
    column 'Extensions' (if at least one application has
    extensions defined).

    Arguments:
        docfile (str): path to output RST file
        tags (list): list of tags to specify which
          applications should be included in the table
    """
    header = ["Platform", "Library type"]
    application_data = []
    has_extensions = False
    has_assays = False
    for application in fetch_application_data(tags=tags, expand=True):
        if application["libraries"][0] == "*":
            # Skip any wildcard libraries
            continue
        platform = f"``{application['platforms'][0]}``"
        library = f"``{application['libraries'][0]}``"
        try:
            extensions = application["extensions"]
            if extensions:
                extensions = ", ".join([f"``{ext}``" for ext in extensions])
                has_extensions = True
            else:
                extensions = ""
        except KeyError:
            extensions = ""
        try:
            assays = application["assays"]
            if assays:
                assays = ", ".join([f"``{assay}``" for assay in assays])
                has_assays = True
            else:
                assays = ""
        except KeyError:
            assays = ""
        application_data.append([assays, platform, library, extensions])
    if has_extensions:
        # Update the header
        header.append("Extensions")
    else:
        # Remove the empty "extensions" column
        application_data = [app[:-1] for app in application_data]
    if has_assays:
        # Update the header
        header = ["Assays"] + header
    else:
        # Remove the empty "assays" column
        application_data = [app[1:] for app in application_data]
    tbl = RstGridTable(application_data)
    with open(docfile, "wt") as fp:
        fp.write("\n".join(tbl.construct_table(header=header)))


def generate_qc_protocols_tbl(docfile):
    """
    Generate table with the QC protocols

    Table lists the QC protocols and has columns
    with protocol name and description.

    Arguments:
        docfile (str): path to output RST file
    """
    protocols_data = []
    for qc_protocol in QC_PROTOCOLS:
        p = fetch_protocol_definition(qc_protocol)
        protocols_data.append([f"``{p.name}``", p.description])
    tbl = RstSimpleTable(protocols_data)
    with open(docfile, "wt") as fp:
        fp.write("\n".join(tbl.construct_table(
            header=["QC protocol", "Description"])))


def generate_api_docs(pkg, docdir):
    """
    Generate developer's API documentation

    Arguments:
        pkg (str): package to generate API docs for
        docdir (str): path to directory where docs
          will be written
    """
    # Ensure the output directory exists
    if not os.path.exists(docdir):
        raise Exception("'%s': does not exist" % docdir)

    # Fetch a list of modules
    modlist = get_modules(pkg)

    # Generate documents for each module
    for modname in modlist:
        docname = modname.replace('.', '_')
        docfile = os.path.join(docdir, f"{docname}.rst")
        print("Generating %s" % docfile)
        with open(docfile, "wt") as doc:
            title = f"``{pkg.__name__}.{modname}``"
            doc.write("""%s
%s

.. automodule:: %s.%s
   :members:
""" % (title, '=' * len(title), pkg.__name__, modname))

    # Generate an index
    api_index = os.path.join(docdir, "index.rst")
    print("Writing %s" % api_index)
    with open(api_index, 'w') as doc:
        doc.write("""=============================
Developers' API documentation
=============================

.. toctree::
   :maxdepth: 2

""")
        for modname in modlist:
            doc.write("   %s\n" % modname.replace('.', '_'))


def generate_utility_docs(utilities, docfile):
    """
    Generate documentation for utility scripts

    Arguments:
      utilities (list): list of utility scripts
      docfile (str): path to output document
    """
    with open(docfile, "wt") as doc:
        doc.write("""
Utilities
=========

.. note::

   This documentation has been auto-generated from the
   command help

In addition to the main ``auto_process.py`` command, a number of utilities
are available:

.. contents:: :local:

""")
    for utility in utilities:
        # Get help text
        help_text = get_help_text(utility)
        if not help_text:
            print("No help text available for %s" % utility)
            continue
        print("Writing documentation for utility '%s'" % utility)
        # Write into the document
        with open(docfile, "a") as doc:
            title = "%s" % utility
            ref = ".. _utilities_%s:" % os.path.splitext(utility)[0]
            doc.write("%s\n\n%s\n%s\n\n::\n\n" % (ref,
                                                  title,
                                                  "*"*len(title)))
            for line in help_text.split('\n'):
                doc.write("    %s\n" % line)


def generate_command_docs(cmds, docfile):
    """
    Generate documentation for 'auto_process' commands
    """
    with open(docfile, "wt") as doc:
        doc.write("""
``auto_process`` commands
=========================

.. note::

   This documentation has been auto-generated from the
   command help

``auto_process.py`` implements the following commands:

.. contents:: :local:

""")
    for subcmd in cmds:
        # Build command
        cmd = ["auto_process.py", subcmd]
        cmdline = " ".join(cmd)
        # Get help text
        help_text = get_help_text(*cmd)
        if not help_text:
            print("No help text available for '%s'" % cmdline)
            continue
        print("Writing documentation for '%s'" % cmdline)
        # Write into the document
        with open(docfile, "a") as doc:
            title = "%s" % subcmd
            ref = ".. _commands_%s:" % subcmd
            doc.write("%s\n\n%s\n%s\n\n::\n\n" % (ref,
                                                  title,
                                                  "*"*len(title)))
            for line in help_text.split('\n'):
                doc.write("    %s\n" % line)


def get_modules(pkg, exclude_tests=True):
    """
    Get a list of modules in a package

    See https://stackoverflow.com/a/1708706/579925

    Arguments:
        pkg (str): package name
        exclude_tests (bool): if True (the default)
          then exclude any modules with names stating
          with 'test'

    Returns:
        list: list of module names within the package
    """
    modlist = []
    for importer, modname, ispkg in walk_packages(path=pkg.__path__,
                                                  prefix=pkg.__name__ + '.',
                                                  onerror=lambda x: None):
        modname = '.'.join(modname.split('.')[1:])
        if modname.startswith("test"):
            continue
        modlist.append(modname)
    return modlist


def get_help_text(cmd, *args):
    """
    Get help text for a command

    Given a command (and optionally any additional arguments),
    return a string containing the help text.

    This is acquired by running the command with the '--help'
    argument appended.

    Arguments:
        cmd (str): command name
        args (list): optional additional arguments
    """
    cmd = [cmd] + list(args) + ["--help"]
    help_text_file = "%s.help" % cmd[0]
    try:
        with open(help_text_file, "wt") as fp:
            subprocess.call(cmd,
                            stdout=fp,
                            stderr=subprocess.STDOUT)
        help_text = open(help_text_file, "rt").read()
        return help_text
    finally:
        os.remove(help_text_file)