#!/usr/bin/env python
#
#     docwriter.py: HTML document generation library
#     Copyright (C) University of Manchester 2017-2021 Peter Briggs
#
#########################################################################

"""
``docwriter`` provides a set of classes that can be used to create
HTML documents programmatically.

Example usage: create a new document:

>>> d = Document("Report")
>>> d.add_css_rule("h1 { color: blue; }")
>>> d.add_css_rule("h2 { color: green; }")
>>> d.add_css_rule("h3 { color: red; }")

Add a section:

>>> s = d.add_section("Introduction")
>>> s.add("This is the introduction")

Add a subsection:

>>> subs = s.add_section("Background")
>>> subs.add("This is the background")

Generate the HTML directly:

>>> d.html()

Write the HTML document to a file:

>>> d.write("report.html")

There are also classes to create lists, tables, images, links and
anchors (aka "targets") within documents.

The full list of available classes are:

- Document: construct HTML documents
- Section: generic document section
- Table: construct HTML tables
- List: create ordered and unordered HTML lists
- Img: embed references to images
- Link: embed references to other documents/items
- Target: create anchor for referencing via a link
- Para: wrap heterogeneous items into a single block
- WarningIcon: create an inline warning icon image
"""

#######################################################################
# Imports
#######################################################################

from bcftbx.htmlpagewriter import HTMLPageWriter

#######################################################################
# Module data
#######################################################################

VALID_CSS_ID_CHARS = "-_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

WARNING_ICON_BASE64 = r'iVBORw0KGgoAAAANSUhEUgAAADIAAAAyCAMAAAAp4XiDAAAAkFBMVEX/oAD/////pQD/ogD/rC//587/nwD/3Kr/7Mz/pwD//PX/nQD/9+n/szX/+O3///3/4rr/1Zr/4rf/v2P/xnX/89//5r//tUH/79j/+fL/rSf/sC7/w23/0I//2qb/vFf/yXz/zYT/t07/uVr/yYb/rAf/s0P/5cT/1pj/3q//t0j/zJH/vlz/yHj/157/qhjbULfEAAAAAXRSTlMAQObYZgAAAAFiS0dEAIgFHUgAAAAJcEhZcwAACxMAAAsTAQCanBgAAAAHdElNRQfjBgsMGw801p5CAAAAcklEQVRIx+3VMQ6AMAxD0fIvwv1vicSABIpI7BYGqOe8oVIatzZzDUfU+aIiiAzuEagGVAOqOQ0sexJDRladkDwkInQTKoT3CTXCJF8kyScu7dgzJFxlRv9K51w4xDh9gek8yuPKwqkkp/iserVK/FfZAOKECZLQoDSVAAAAAElFTkSuQmCC'

#######################################################################
# Classes
#######################################################################

class Document:
    """
    Utility class for constructing documents

    Basic usage:

    >>> d = Document("Report")
    >>> d.add_css_rule("h1 { color: red; }")
    >>> s = d.add_section("Introduction")
    >>> s.write("report.html")

    The following properties are available:

    - title: returns the title string
    - css_rules: returns the CSS rules as a
      list

    The 'html' method returns the document
    body rendered as HTML; the 'write'
    method writes the full HTML document to
    file, including the document header, CSS
    rules etc.

    Arguments:
      title (str): title for the document
    """
    def __init__(self,title=None):
        """
        Create a new Document instance
        """
        self._title = title
        self._sections = []
        self._css_rules = []

    @property
    def title(self):
        """
        Get the document's title
        """
        return self._title

    @property
    def css_rules(self):
        """
        Get the CSS rules as a list
        """
        return [rule for rule in self._css_rules]

    def add_section(self,title=None,section=None,name=None,
                    css_classes=None):
        """
        Add a section to the document

        If an existing section is supplied then
        this is appended to the document; otherwise
        a new section is created with the title
        and name supplied.

        Arguments:
          title (str): title for the new section
          section (Section): an existing Section
            instance which will be appended to the
            document
          name (str): internal name for the section;
            if not supplied then the name will be
            generated automatically from the title
          css_classes (list): list or iterable with
            names of CSS classes to associate with
            the new section.

        Returns:
          Section: the new or supplied section as
            appropriate.
        """
        if section is None:
            section = Section(title=title,name=name,
                              css_classes=css_classes)
        self._sections.append(section)
        return section

    def add_css_rule(self,css_rule):
        """
        Add a CSS rule block to the document

        Arguments:
          css_rule (str): block of text to append
            to the CSS rules to write to the final
            document
        """
        self._css_rules.append(css_rule)

    def html(self):
        """
        Generate HTML version of the document contents

        Note that this only generates the "body"
        HTML code; the header and footers will not
        be returned (so this will not be the full
        HTML that is written by the 'write' method).

        Returns:
          String: HTML representation of the document
            content.
        """
        html = []
        if self._title is not None:
            html.append("<h1>%s</h1>" % self._title)
        for section in self._sections:
            try:
                html.append(section.html())
            except AttributeError as ex:
                html.append("<p>Failed to render section: %s</p>" % ex)
        return '\n'.join(html)

    def write(self,outfile):
        """
        Write document contents to a file

        Arguments
          outfile (str): path to file to write
            HTML document to
        """
        html = HTMLPageWriter(self._title)
        for css_rule in self._css_rules:
            html.addCSSRule(css_rule)
        html.add(self.html())
        html.write("%s" % outfile)

class Section:
    """
    Class representing a generic document section

    A 'section' is a container for arbitrary content
    which can include text, images, links, tables,
    or other sections (which are considered to be
    subsections of the containing section).

    A section typically has a title, however
    'anonymous' (untitled) sections are also allowed.

    Normally sections are created via the 'add_section'
    methods of a 'Document' or 'Section' instance,
    however they can also be created independently
    and then added (again via the 'add_section' method).

    Content is added to the section using the 'add'
    method, which can invoked multiple times to
    append additional content.

    When the content is rendered, each item is
    rendered to html via its 'html' method; or if the
    item doesn't have a 'html' method then it will be
    rendered as a string and wrapped in '<p>' tags.

    Example usage:

    >>> d = Document("Report")
    >>> abstract = d.section("Abstract")
    >>> abstract.add("This is the abstract","It's quite short")
    >>> paper = d.section()
    >>> methods = paper.add_subsection("Methods")
    >>> methods.add("This describes the methods used")
    >>> results = paper.add_subsection("Results)
    >>> result.add("Some nice figures")
    >>> result.add(Img("graph.png"))

    Arguments:
      title (str): title text
      name (str): name used for the 'id' of
        the section
      level (int): section heading level
        (defaults to 2)
      css_classes (list): list or iterable with
        names of CSS classes to associate with
        the section
    """
    def __init__(self,title=None,name=None,level=2,
                 css_classes=None):
        """
        Create a new Section instance
        """
        self._title = title
        self._name = name
        self._content = []
        self._css_classes = []
        self._level = level
        if css_classes:
            self.add_css_classes(*css_classes)

    @property
    def name(self):
        """
        Return the name (id) for the section

        """
        if self._name is not None:
            return self._name
        if self._title is not None:
            return sanitize_css_string(self._title)
        else:
            return self._title

    @property
    def title(self):
        """
        Return the title for the section
        """
        return self._title

    @property
    def level(self):
        """
        Return the level of the section
        """
        return self._level

    @property
    def css_classes(self):
        """
        Return the CSS classes as a list
        """
        return [cls for cls in self._css_classes]

    def add_css_classes(self,*classes):
        """
        Associate CSS classes with the section

        Arguments:
          classes (str): the names of one or
            more CSS classes to associate with
            the section when it is rendered
        """
        for css_class in classes:
            self._css_classes.append(css_class)

    def add(self,*args):
        """
        Add content to the section

        Arguments:
          args (object): one or more arbitrary
            items of content to append to the
            section; items must implement either
            a 'html' or '__str__' method
        """
        for content in args:
            self._content.append(content)

    def add_subsection(self,title=None,section=None,name=None,
                       css_classes=None):
        """
        Add subsection within the section

        If an existing section is supplied then this
        is appended to the content; otherwise a new
        Section instance is created and appended.

        Arguments:
          title (str): title text
          section (Section): if supplied then must be
            a Section instance which is inserted into
            the current Section as-is
          name (str): name used for the 'id' of
            the section
          css_classes (list): list or iterable with
            names of CSS classes to associate with
            the new section.

        Returns:
          Section: the appended section.
        """
        if section is None:
            subsection = Section(title=title,name=name,
                                 level=self._level+1,
                                 css_classes=css_classes)
        else:
            subsection = section
        self.add(subsection)
        return subsection

    def html(self):
        """
        Generate HTML version of the section

        Returns:
          String: HTML representation of the section
            content.
        """
        if self._title is None and \
           not self._content and \
           not self._css_classes:
            return ""
        div = "<div"
        if self.name:
            div += " id='%s'" % self.name
        if self._css_classes:
            div += " class='%s'" % ' '.join(self._css_classes)
        div += ">"
        html = [div,]
        if self._title is not None:
            html.append("<h%d>%s</h%d>" % (self._level,
                                           self._title,
                                           self._level))
        for content in self._content:
            try:
                html.append(content.html())
            except AttributeError as ex:
                html.append("<p>%s</p>" % str(content))
        html.append('</div>')
        return '\n'.join(html)

class Table:
    """
    Utility class for constructing HTML tables

    Usage examples:

    - Table with two columns with titles
      "Key" and "Value", and with one row of
      data:

      >>> t = Table(('Key','Value'))
      >>> t.add_row(Key="Employee name",
      ...           Value="John Doe")
      >>> print(t.html())
      <table>
      <tr><th>Key</th><th>Value</th></tr>
      <tr><td>Employee name</td><td>John Doe</td></tr>
      </table>

    - Table with three columns, with aliases
      for column names:

      >>> t = Table(('name','date','result'),
      ...           name="Experiment",
      ...           date="Date run",
      ...           result="Final result")
      >>> t.add_row(name="Test setup",
      ...           date="10/10/2017",
      ...           result="Ok")
      >>> print(t.html())
      <table>
      <tr><th>Experiment</th><th>Date run</th><th>Final result</th></tr>
      <tr><td>Test setup</td><td>10/10/2017</td><td>Ok</td></tr>
      </table>

    - Add a column to the previous table, and
      set values:

      >>> t.append_columns('verified',verified="Is verified?")
      >>> t.set_value(0,"verified","no")
      >>> print(t.html())
      <table>
      <tr><th>Experiment</th><th>Date run</th><th>Final result</th><th>Is verified?</th></tr>
      <tr><td>Test setup</td><td>10/10/2017</td><td>Ok</td><td>no</td></tr>
      </table>

    The Table class provides the following properties:

    - nrows: number of rows in the table

    It provides the following methods:

    - add_css_classes: associate CSS classes with the table
    - no_header: turn off writing table header in output HTML
    - append_columns: add additional columns to the table
    - add_row: append a row to the table
    - set_value: set the value in a table cell

    Arguments:
      columns (list): list of column ids
      aliases (mapping): optional, mapping of
        column ids to aliases (i.e. actual
        names)
    """
    def __init__(self,columns,**aliases):
        """
        Create a new Table instance
        """
        self._columns = [x for x in columns]
        self._rows = []
        self._column_names = dict(aliases)
        self._css_classes = []
        self._css_classes_columns = {}
        self._output_header = True

    @property
    def nrows(self):
        """
        Return the number of rows in the table
        """
        return len(self._rows)

    def add_css_classes(self,*classes,**kws):
        """
        Associate CSS classes with the table

        By default the supplied CSS classes
        associated with the whole table. If
        the 'column' keyword is supplied then
        the classes will be associated only
        with the table cells in that column.

        Arguments:
          classes (list): one or more classes
            to associate with the table
          column (str): optional, if supplied
            then should specify a column name
            with which the classes will be
            associated
        """
        css_classes = self._css_classes
        for kw in kws:
            if kw == "column":
                col = kws[kw]
                if col not in self._columns:
                    raise KeyError("'%s'" % col)
                if col not in self._css_classes_columns:
                    self._css_classes_columns[col] = []
                css_classes = self._css_classes_columns[col]
            else:
                raise TypeError(
                    "add_css_classes() got an "
                    "unexpected argument '%s'"
                    % kw)
        for css_class in classes:
            css_classes.append(css_class)

    def no_header(self):
        """
        Don't write the table header on output
        """
        self._output_header = False

    def append_columns(self,*columns,**aliases):
        """
        Add a new columns to the table

        Arguments:
          columns (list): list of column ids
          aliases (mapping): optional, mapping of
            column ids to aliases (i.e. actual names)
        """
        for col in columns:
            if col in self._columns:
                raise KeyError("Column with id '%s' already defined"
                               % col)
            self._columns.append(col)
        for col in aliases:
            self._column_names[col] = aliases[col]

    def add_row(self,**kws):
        """
        Add (append) a row to the table

        Appends a new empty row to the end of
        the table, and returns its row index.

        Optionally values can also be assigned to
        columns in the new row, e.g.

        t.add_row(col1="value1",col2="value2",...)

        Arguments:
          kws (mapping): set of key=value
            assignments, defining the values to
            assign to columns in the new row

        Returns:
          Integer: the row index of the appended
            row (for use e.g. in 'set_value')
        """
        self._rows.append({})
        n = len(self._rows)-1
        for key in kws:
            self.set_value(n,key,kws[key])
        return n

    def set_value(self,row,key,value):
        """
        Set the value of a cell in the table

        Arguments:
          row (int): index of the row to update
          key (str): id of the column to set the
            value for in this row
          value (str): value to assign to the
            table cell indicated by row/column
        """
        if key not in self._columns:
            raise KeyError("Key '%s' not found" % key)
        self._rows[row][key] = value

    def html(self,css_id=None):
        """
        Generate HTML version of the table contents

        Arguments:
          css_id (str): optional, string to write
            into the CSS ``id`` attribute

        Returns:
          String: HTML representation of the table.
        """
        html = []
        # Opening tag
        table_tag = []
        table_tag.append("<table")
        if css_id is not None:
            table_tag.append(" id='%s'" % css_id)
        if self._css_classes:
            table_tag.append(" class='%s'" % ' '.join(self._css_classes))
        table_tag.append(">")
        html.append(''.join(table_tag))
        # Header
        if self._output_header:
            header = []
            header.append("<tr>")
            for col in self._columns:
                try:
                    css_classes = " class='%s'" % \
                                  ' '.join(self._css_classes_columns[col])
                except KeyError:
                    css_classes = ''
                try:
                    col_name = self._column_names[col]
                except KeyError:
                    col_name = col
                header.append("<th%s>%s</th>" % (css_classes,
                                                 str(col_name)))
            header.append("</tr>")
            html.append(''.join(header))
        # Body
        for row in self._rows:
            line = []
            line.append("<tr>")
            for col in self._columns:
                try:
                    css_classes = " class='%s'" % \
                                  ' '.join(self._css_classes_columns[col])
                except KeyError:
                    css_classes = ''
                try:
                    value = row[col].html()
                except KeyError:
                    value = '&nbsp;'
                except AttributeError:
                    value = row[col]
                line.append("<td%s>%s</td>" % (css_classes,
                                               value))
            line.append("</tr>")
            html.append(''.join(line))
        # Finish
        html.append("</table>")
        return '\n'.join(html)

class List:
    """
    Utility class for creating ordered and unordered lists

    Example usage:

    >>> lst = List(ordered=True)
    >>> lst.add_item("List item")
    >>> lst.html()
    "<ol><li>List item</li></ol>"

    """
    def __init__(self,name=None,ordered=False):
        """
        Create a new List instance

        Arguments:
          name (str): string to use as the 'id'
            attribute for the <ul>/<ol> tag
          ordered (boolean): if True then create
            an ordered (i.e. numbered) list;
            otherwise make an unnumbered list (the
            default)

        """
        self._name = name
        self._ordered = bool(ordered)
        self._items = []

    def add_item(self,*content):
        """
        Append an item to the list

        The item can consist of one or more
        objects (for example strings and other
        docwriter objects such as Lists, Links
        etc), which will be concatenated after
        rendering.

        Arguments:
          contents (sequence): one or more
            objects to add to the list item
        """
        item = [c for c in content]
        self._items.append(item)

    def html(self):
        """
        Generate HTML version of the list

        Returns:
          String: HTML representation of the list.
        """
        # Empty list?
        if not self._items:
            return ""
        # List type
        if self._ordered:
            tag = "ol"
        else:
            tag = "ul"
        # Build the html
        html = []
        html.append("<%s" % tag)
        if self._name:
            html.append(" id='%s'" % self._name)
        html.append(">")
        # Add items
        for item in self._items:
            html.append("<li>")
            for i in item:
                try:
                    html.append(i.html())
                except AttributeError:
                    html.append(str(i))
            html.append("</li>")
        # Close the list
        html.append("</%s>" % tag)
        return "".join(html)

class Img:
    """
    Utility class for embedding <img> tags

    Example usage:

    >>> img = Img('picture.png',width=150)
    >>> img.html()
    "<img src='picture.png' width='150' />"

    Arguments:
      src (str): string to use as the 'src'
        attribute for the <img.../> tag
      name (str): string to use as the 'id'
        attribute for the <img.../> tag
      height (int): optional height (pixels)
      width (int): optional width (pixels)
      href (str): if specified then the <img.../>
        will be wrapped in <a href..>...</a> with
        this used as the link target
      alt (str): if specified then used as the
        'alternative text' ('alt' attribute
        for <img.../> tag)
      title (str): if specified then assigned
        to the <img../> tag's 'title' attribute
    """
    def __init__(self,src,name=None,height=None,width=None,href=None,
                 alt=None,title=None):
        """
        Create a new Img instance
        """
        self._src = src
        self._height = height
        self._width = width
        self._name = name
        self._target = href
        self._alt = alt
        self._title = title

    @property
    def name(self):
        """
        Return the name (id) for the image
        """
        return self._name

    def html(self):
        """
        Generate HTML version of the image tag

        Returns:
          String: HTML representation of the image.
        """
        # Build the tag contents
        html = []
        html.append("<img")
        if self._name:
            html.append("id='%s'" % self._name)
        html.append("src='%s'" % self._src)
        # Optional height and width
        if self._height:
            html.append("height='%s'" % self._height)
        if self._width:
            html.append("width='%s'" % self._width)
        # Optional alt text
        if self._alt:
            html.append("alt='%s'" % self._alt)
        # Optional title
        if self._title:
            html.append("title='%s'" % self._title)
        # Close the tag
        html.append("/>")
        # Wrap in a hef
        if self._target:
            return Link(" ".join(html),self._target).html()
        else:
            return " ".join(html)

    def __repr__(self):
        return self.html()

class Link:
    """
    Utility class for embedding <a href=...> tags

    Example usage:

    >>> ahref = Link('My report','report.html')
    >>> ahref.html()
    "<a href='report.html'>My report</a>"

    The link target can be an object:

    >>> sect = Section("New section",name='new_section')
    >>> ahref = Link("New Section",sect)
    >>> ahref.html()
    "<a href='#new_section'>New Section</a>"

    The target can be seen via the 'href' property:

    >>> ahref.href
    '#new_section'

    Arguments:
      text (str): text to display for the link
      target (Object): target to link to; if not
        supplied then 'text' will be used as the
        link target
    """
    def __init__(self,text,target=None):
        """
        Create a new Link instance
        """
        self._text = text
        if target is None:
            self._target = text
        else:
            self._target = target

    @property
    def href(self):
        """
        Show the link target text
        """
        try:
            return '#%s' % self._target.name
        except AttributeError:
            return str(self._target)

    def html(self):
        """
        Generate HTML version of the link

        Returns:
          String: HTML representation of the link.
        """
        # Build the tag contents
        return "<a href='%s'>%s</a>" % (self.href,self._text)

    def __repr__(self):
        return self.html()

class Target:
    """
    Utility class for embedding <a id=... /> tags

    Example usage:

    >>> tgt = Target('my_target')
    >>> tgt.html()
    "<a id='my_target' />"

    The target name/id can be retrieved using the
    'name' property. It can also be supplied directly
    to a Link object:

    >>> ahref = Link("My target",tgt)
    >>> ahref.html()
    "<a href='#my_target'>My target</a>"

    Arguments:
      name (str): name (i.e. 'id') for the
        target
    """
    def __init__(self,name):
        """
        Create a new Target instance
        """
        self._name = name

    @property
    def name(self):
        """
        Return the name (id) for the target
        """
        return self._name

    def html(self):
        """
        Generate HTML version of the target

        Returns:
          String: HTML representation of the
            target.
        """
        # Build the anchor
        return "<a id='%s' />" % self._name

class Para:
    """
    Utility to wrap heterogeneous items into a single block

    The Para class provides a way to put several
    different items (e.g. plain text, links, images)
    into a single "block" (by default a '<p>'-delimited
    paragraph).

    Example usage:

    >>> s = Section("QC reports")
    >>> s.add(Para("%s" % project.name,
    ...            Link("qc_report.html")))

    The Para wrapper prevents each plain text (i.e.
    string-like) content element being wrapped in its
    own <p>...</p> pair; instead only the completely
    assembled Para HTML is enclosed in <p>...</p>.

    Arguments:
      items (sequence): optional, set of
        items to add to the block
      css_classes (list): list or iterable with
        names of CSS classes to associate with
        the section
    """
    def __init__(self,*items,**kws):
        """
        Create a new Para instance
        """
        self._content = [x for x in items]
        self._delimiter = " "
        self._tag = "p"
        self._css_classes = []
        for kw in kws:
            if kw == 'css_classes':
                css_classes = kws['css_classes']
                if css_classes:
                    self.add_css_classes(*css_classes)
            else:
                raise TypeError("__init__() got an unexpected "
                                "keyword argument '%s'" % kw)

    def add_css_classes(self,*classes):
        """
        Associate CSS classes with the paragraph

        Arguments:
          classes (str): the names of one or
            more CSS classes to associate with
            the section when it is rendered
        """
        for css_class in classes:
            self._css_classes.append(css_class)

    def add(self,*items):
        """
        Add (append) content to the block

        Arguments:
          items (sequence): optional, set of
            items to append to the bloc
        """
        for item in items:
            self._content.append(item)

    def html(self):
        """
        Generate HTML version of the block

        Returns:
          String: HTML representation of the block.
        """
        if not self._content:
            return ""
        html = []
        for content in self._content:
            try:
                html.append(content.html())
            except AttributeError as ex:
                html.append("%s" % str(content))
        if self._css_classes:
            classes = " class='%s'" % ' '.join(self._css_classes)
        else:
            classes = ''
        return "<%s%s>%s</%s>" % (self._tag,
                                  classes,
                                  self._delimiter.join(html),
                                  self._tag)

    def __bool__(self):
        """
        Para instance is True if has content, False otherwise
        """
        return bool(self._content)

    def __nonzero__(self):
        """
        Para instance is True if has content, False otherwise
        """
        return self.__bool__()

class WarningIcon(Img):
    """
    Create image with an inline warning icon

    Arguments:
      title (str): optional, content to mark with
        the warning
      size (int): optional height/width specifier for
        the icon (defaults to 25)
    """
    def __init__(self,title=None,size=25):
        """
        Create new WarningIcon instance
        """
        Img.__init__(self,
                     "data:image/png;base64,%s" % WARNING_ICON_BASE64,
                     title=title,
                     height=size,
                     width=size)

#######################################################################
# Functions
#######################################################################

def sanitize_css_string(s):
    """
    Remove or replace invalid (non-CSS) characters in a string

    Arguments:
      s (str): string to be sanitized

    Returns:
      String: version of the original string with whitespace
        converted to underscores, and all other invalid
        characters removed.
    """
    return ''.join(filter(lambda c: c in VALID_CSS_ID_CHARS,
                          str(s).replace(' ','_').replace('\t','_')))
