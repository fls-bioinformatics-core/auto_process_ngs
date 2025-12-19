#!/usr/bin/env python
#
#     docs: helpers for auto-process-ngs Sphinx documentation generation
#     Copyright (C) University of Manchester 2025 Peter Briggs
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

Functions:

* ``generate_api_docs``: generate developer's API documentation
* ``generate_utility_docs``: generate documentation for utilities
* ``generate_command_docs``: generate documentation for 'auto_process' commands
* ``get_modules``: list Python modules in a package
* ``get_help_text``: fetch text from running program with '--help'

"""

import os
import subprocess
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


def generate_api_docs(pkg, docdir):
    """
    Generate developer's API documentation
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