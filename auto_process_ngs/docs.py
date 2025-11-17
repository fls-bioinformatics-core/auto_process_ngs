#!/usr/bin/env python
#
#     docs: helpers for auto-process-ngs Sphix documentation generation
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Provides utility classes and functions to help with generating Sphinx
documentation for the ``auto_process_ngs`` package.

The helpers are intended to be called from the Sphinx ``conf.py`` file,
to generate content such as tables of applications and protocols which
are then included in the documentation source files.

Classes:

* ``RstTable``: class for making reStructuredText tables

"""


class RstTable:
    """
    Class for making reStructuredText tables.

    Usage:

    >>> table_data = [
    ...     ['col1_row1','col2_row1','col3_row1'],
    ...     ['col1_row2','col2_row2','col3_row2'],
    ...     ['col1_row3','col2_row3','col3_row3'],
    ... ]
    >>> header = ['Column 1','Column 2','Column 3']
    >>> table = RstTable(table_data)
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
    >>> table = RstTable(table_data)
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

    def construct_table(self, header=None):
        """
        Returns reStructuredText table
        """
        # Collect the field widths
        field_widths = self.get_field_widths(header=header)
        # Start constructing the table
        table = []
        # Add table header
        if header:
            for line in self.make_header(header, field_widths):
                table.append(line)
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
            table.append(divider)
            line = ""
            for i, col in enumerate(row):
                if col:
                    # Pad value
                    line += "| " + col + " "*(field_widths[i]-len(col)+1)
                else:
                    # Empty cell
                    line += "| " + " "*(field_widths[i]+1)
            line += "|"
            table.append(line)
        # Closing divider
        table.append(self.make_divider(field_widths))
        return table