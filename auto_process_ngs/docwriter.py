#!/usr/bin/env python
#
# document generation library

from bcftbx.htmlpagewriter import HTMLPageWriter

VALID_CSS_ID_CHARS = "-_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

class Document:
    """
    Utility class for constructing documents

    """
    def __init__(self,title=None):
        self._title = title
        self._sections = []
        self._css_rules = []

    def add_section(self,title=None,section=None,name=None):
        """
        Add a new section

        """
        if section is None:
            section = Section(title=title,name=name)
        self._sections.append(section)
        return section

    def add_css_rule(self,css_rule):
        """
        """
        self._css_rules.append(css_rule)

    def html(self):
        """
        Generate HTML version of the document contents

        """
        html = []
        if self._title is not None:
            html.append("<h1>%s</h1>" % self._title)
        for section in self._sections:
            try:
                html.append(section.html())
            except AttributeError,ex:
                html.append("<p>Failed to render section: %s</p>" % ex)
        return '\n'.join(html)

    def write(self,outfile):
        """
        Write document contents to a file

        """
        html = HTMLPageWriter(self._title)
        for css_rule in self._css_rules:
            html.addCSSRule(css_rule)
        html.add(self.html())
        html.write("%s" % outfile)

class Section:
    """
    Class representing a generic document section

    """
    def __init__(self,title=None,name=None,level=2):
        """
        Create a new Section instance

        Arguments:
          title (str): title text
          name (str): name used for the 'id' of
            the section
          level (int): section heading level
            (defaults to 2)

        """
        self._title = title
        self._name = name
        self._content = []
        self._css_classes = []
        self._level = level

    @property
    def name(self):
        """
        Return the name (id) for the section

        """
        if self._name is not None:
            return self._name
        if self._title is not None:
            return filter(lambda c: c in VALID_CSS_ID_CHARS,
                          str(self._title).replace(' ','_').replace('\t','_'))
        else:
            return self._title

    @property
    def title(self):
        """
        Return the title for the section
        """
        return self._title

    def add_css_classes(self,*classes):
        """
        Associate CSS classes with the section

        """
        for css_class in classes:
            self._css_classes.append(css_class)

    def add(self,*args):
        """
        Add content to the section

        """
        for content in args:
            self._content.append(content)

    def add_subsection(self,title=None,section=None,name=None):
        """
        Add subsection within the section

        Arguments:
          title (str): title text
          section (Section): if supplied then must be
            a Section instance which is inserted into
            the current Section as-is
          name (str): name used for the 'id' of
            the section

        """
        if section is None:
            subsection = Section(title=title,name=name,level=self._level+1)
        else:
            subsection = section
        self.add(subsection)
        return subsection

    def html(self):
        """
        Generate HTML version of the section

        """
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
            except AttributeError,ex:
                html.append("<p>%s</p>" % str(content))
        html.append('</div>')
        return '\n'.join(html)

class Table:
    """
    Utility class for constructing tables for output

    Example usage:

    >>> t = Table('Key','Value')

    """
    def __init__(self,columns,**kws):
        """
        Create a new ReportTable instance

        Arguments:
          columns (list): list of column ids
          kws (mapping): optional, mapping of
            column ids to actual names

        """
        self._columns = [x for x in columns]
        self._rows = []
        self._column_names = dict(kws)
        self._css_classes = []
        self._output_header = True

    @property
    def nrows(self):
        """
        Return the number of rows in the table
        """
        return len(self._rows)

    def add_css_classes(self,*classes):
        """
        Associate CSS classes with the table

        """
        for css_class in classes:
            self._css_classes.append(css_class)

    def no_header(self):
        """
        Don't write the table header on output

        """
        self._output_header = False

    def append_columns(self,*columns,**kws):
        """
        Add a new columns to the table

        Arguments:
          columns (list): list of column ids
          kws (mapping): optional, mapping of
            column ids to actual names

        """
        for col in columns:
            if col in self._columns:
                raise KeyError("Column with id '%s' already defined"
                               % col)
            self._columns.append(col)
        for col in kws:
            self._column_names[col] = kws[col]

    def add_row(self,**kws):
        """
        Add a row to the table

        """
        self._rows.append({})
        n = len(self._rows)-1
        for key in kws:
            self.set_value(n,key,kws[key])
        return n

    def set_value(self,row,key,value):
        """
        Set the value of a field in a row

        """
        if key not in self._columns:
            raise KeyError("Key '%s' not found" % key)
        self._rows[row][key] = value

    def html(self,css_id=None):
        """
        Generate HTML version of the table contents

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
                    col_name = self._column_names[col]
                except KeyError:
                    col_name = col
                header.append("<th>%s</th>" % col_name)
            header.append("</tr>")
            html.append(''.join(header))
        # Body
        for row in self._rows:
            line = []
            line.append("<tr>")
            for col in self._columns:
                try:
                    value = row[col].html()
                except KeyError:
                    value = '&nbsp;'
                except AttributeError:
                    value = row[col]
                line.append("<td>%s</td>" % value)
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

        """
        item = [c for c in content]
        self._items.append(item)

    def html(self):
        """
        Generate HTML version of the list

        """
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

    """
    def __init__(self,src,name=None,height=None,width=None,href=None,
                 alt=None):
        """
        Create a new Img instance

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
            for <img.../> tag

        """
        self._src = src
        self._height = height
        self._width = width
        self._name = name
        self._target = href
        self._alt = alt

    @property
    def name(self):
        """
        Return the name (id) for the image

        """
        return self._name

    def html(self):
        """
        Generate HTML version of the image tag

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
        # Close the tag
        html.append("/>")
        # Wrap in a hef
        if self._target:
            return Link(" ".join(html),self._target).html()
        else:
            return " ".join(html)

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

    """
    def __init__(self,text,target=None):
        """
        Create a new Link instance

        Arguments:
          text (str): text to display for the link
          target (Object): target to link to; if not
            supplied then 'text' will be used as the
            link target

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

    """
    def __init__(self,name):
        """
        Create a new Target instance

        Arguments:
          name (str): name (i.e. 'id') for the
            target

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

        """
        # Build the anchor
        return "<a id='%s' />" % self._name
