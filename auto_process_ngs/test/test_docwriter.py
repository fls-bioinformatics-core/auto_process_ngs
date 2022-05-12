#######################################################################
# Tests for docwriter.py module
#######################################################################

import unittest
import os
import tempfile
import shutil
from auto_process_ngs.docwriter import Document
from auto_process_ngs.docwriter import Section
from auto_process_ngs.docwriter import Table
from auto_process_ngs.docwriter import List
from auto_process_ngs.docwriter import Img
from auto_process_ngs.docwriter import Link
from auto_process_ngs.docwriter import Target
from auto_process_ngs.docwriter import Para
from auto_process_ngs.docwriter import WarningIcon
from auto_process_ngs.docwriter import DocumentIcon
from auto_process_ngs.docwriter import sanitize_css_string
from auto_process_ngs.docwriter import WARNING_ICON_BASE64
from auto_process_ngs.docwriter import DOCUMENT_ICON_BASE64

# Unit tests

class TestDocument(unittest.TestCase):
    """
    Tests for the Document class
    """
    def setUp(self):
        # Create a temp working dir
        self.dirn = tempfile.mkdtemp(suffix='TestDocument')

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.dirn)

    def test_empty_document(self):
        d = Document()
        self.assertEqual(d.title,None)
        self.assertEqual(d.css_rules,[])
        self.assertEqual(d.html(),"")

    def test_document_with_title(self):
        d = Document("Document with title")
        self.assertEqual(d.title,"Document with title")
        self.assertEqual(d.css_rules,[])
        self.assertEqual(d.html(),
                         "<h1>Document with title</h1>")

    def test_document_with_sections(self):
        d = Document("Document with sections")
        s1 = d.add_section("First section")
        self.assertTrue(isinstance(s1,Section))
        self.assertEqual(d.html(),
                         "<h1>Document with sections</h1>\n"
                         "<div id='First_section'>\n"
                         "<h2>First section</h2>\n"
                         "</div>")
        s2 = d.add_section("Second section")
        self.assertTrue(isinstance(s2,Section))
        self.assertEqual(d.html(),
                         "<h1>Document with sections</h1>\n"
                         "<div id='First_section'>\n"
                         "<h2>First section</h2>\n"
                         "</div>\n"
                         "<div id='Second_section'>\n"
                         "<h2>Second section</h2>\n"
                         "</div>")

    def test_document_add_existing_section(self):
        d = Document("Document with existing section")
        s = Section("External section")
        d.add_section(section=s)
        self.assertEqual(d.html(),
                         "<h1>Document with existing section</h1>\n"
                         "<div id='External_section'>\n"
                         "<h2>External section</h2>\n"
                         "</div>")

    def test_document_add_section_with_css_classes(self):
        d = Document("Document with section with CSS classes")
        s = d.add_section("New section",css_classes=("summary","new",))
        self.assertEqual(d.html(),
                         "<h1>Document with section with CSS classes</h1>\n"
                         "<div id='New_section' class='summary new'>\n"
                         "<h2>New section</h2>\n"
                         "</div>")

    def test_document_with_css_rules(self):
        d  = Document("CSS rules")
        d.add_css_rule("h1 { color: black; }")
        d.add_css_rule("h2 { color: grey; }")
        self.assertEqual(d.css_rules,
                         ["h1 { color: black; }",
                          "h2 { color: grey; }"])
        self.assertEqual(d.html(),
                         "<h1>CSS rules</h1>")

    def test_document_write_to_file(self):
        d = Document("Test Document")
        outfile = os.path.join(self.dirn,"test.html")
        self.assertFalse(os.path.exists(outfile))
        d.write(outfile)
        self.assertTrue(os.path.exists(outfile))
        with open(outfile,'r') as fp:
            html = fp.read()
            self.assertEqual(html,
                             "<html>\n"
                             "<head>\n"
                             "<title>Test Document</title>\n"
                             "</head>\n"
                             "<body>\n"
                             "<h1>Test Document</h1></body>\n"
                             "</html>\n")

    def test_document_write_to_file_with_javascript(self):
        d = Document("Test Document")
        d.add_javascript("// Placeholder for script code")
        outfile = os.path.join(self.dirn,"test.html")
        self.assertFalse(os.path.exists(outfile))
        d.write(outfile)
        self.assertTrue(os.path.exists(outfile))
        with open(outfile,'r') as fp:
            html = fp.read()
            self.assertEqual(html,
                             "<html>\n"
                             "<head>\n"
                             "<title>Test Document</title>\n"
                             "<script language='javascript' type='text/javascript'><!--\n"
                             "// Placeholder for script code\n"
                             "--></script>\n"
                             "</head>\n"
                             "<body>\n"
                             "<h1>Test Document</h1></body>\n"
                             "</html>\n")

class TestSection(unittest.TestCase):
    """
    Tests for the Section class
    """
    def test_empty_section(self):
        s = Section()
        self.assertEqual(s.title,None)
        self.assertEqual(s.name,None)
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),"")

    def test_empty_section_with_css_class(self):
        s = Section(css_classes=("clear",))
        self.assertEqual(s.title,None)
        self.assertEqual(s.name,None)
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),
                         "<div class='clear'>\n"
                         "</div>")

    def test_empty_section_with_style_attribute(self):
        s = Section(style="display: block;")
        self.assertEqual(s.title,None)
        self.assertEqual(s.name,None)
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),
                         "<div style='display: block;'>\n"
                         "</div>")

    def test_section_no_content(self):
        s = Section("Empty section")
        self.assertEqual(s.title,"Empty section")
        self.assertEqual(s.name,"Empty_section")
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),
                         "<div id='Empty_section'>\n"
                         "<h2>Empty section</h2>\n"
                         "</div>")

    def test_section_with_content(self):
        s = Section("Section with content")
        s.add("Some stuff")
        self.assertEqual(s.title,"Section with content")
        self.assertEqual(s.name,"Section_with_content")
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),
                         "<div id='Section_with_content'>\n"
                         "<h2>Section with content</h2>\n"
                         "<p>Some stuff</p>\n"
                         "</div>")
        s.add("More stuff","and still more!")
        self.assertEqual(s.html(),
                         "<div id='Section_with_content'>\n"
                         "<h2>Section with content</h2>\n"
                         "<p>Some stuff</p>\n"
                         "<p>More stuff</p>\n"
                         "<p>and still more!</p>\n"
                         "</div>")

    def test_section_with_subsection(self):
        s = Section("Section with subsection")
        sub = s.add_subsection("Subsection")
        self.assertTrue(isinstance(sub,Section))
        self.assertEqual(s.title,"Section with subsection")
        self.assertEqual(s.name,"Section_with_subsection")
        self.assertEqual(s.level,2)
        self.assertEqual(sub.level,3)
        self.assertEqual(s.html(),
                         "<div id='Section_with_subsection'>\n"
                         "<h2>Section with subsection</h2>\n"
                         "<div id='Subsection'>\n"
                         "<h3>Subsection</h3>\n"
                         "</div>\n"
                         "</div>")

    def test_section_with_subsection_and_css_classes(self):
        s = Section("Section with subsection and CSS classes")
        sub = s.add_subsection("Subsection",css_classes=("subsection","new"))
        self.assertTrue(isinstance(sub,Section))
        self.assertEqual(s.title,"Section with subsection and CSS classes")
        self.assertEqual(s.name,"Section_with_subsection_and_CSS_classes")
        self.assertEqual(s.level,2)
        self.assertEqual(sub.level,3)
        self.assertEqual(s.html(),
                         "<div id='Section_with_subsection_and_CSS_classes'>\n"
                         "<h2>Section with subsection and CSS classes</h2>\n"
                         "<div id='Subsection' class='subsection new'>\n"
                         "<h3>Subsection</h3>\n"
                         "</div>\n"
                         "</div>")

    def test_section_with_subsection_and_style_attribute(self):
        s = Section("Section with subsection and style attribute")
        sub = s.add_subsection("Subsection",style="display: block;")
        self.assertTrue(isinstance(sub,Section))
        self.assertEqual(s.title,"Section with subsection and style attribute")
        self.assertEqual(s.name,"Section_with_subsection_and_style_attribute")
        self.assertEqual(s.level,2)
        self.assertEqual(sub.level,3)
        self.assertEqual(s.html(),
                         "<div id='Section_with_subsection_and_style_attribute'>\n"
                         "<h2>Section with subsection and style attribute</h2>\n"
                         "<div id='Subsection' style='display: block;'>\n"
                         "<h3>Subsection</h3>\n"
                         "</div>\n"
                         "</div>")

    def test_section_with_external_subsection(self):
        s = Section("Section with external subsection")
        sub = Section("External subsection")
        self.assertEqual(sub.level,2)
        s.add_subsection(section=sub)
        self.assertTrue(isinstance(sub,Section))
        self.assertEqual(s.title,"Section with external subsection")
        self.assertEqual(s.name,"Section_with_external_subsection")
        self.assertEqual(s.html(),
                         "<div id='Section_with_external_subsection'>\n"
                         "<h2>Section with external subsection</h2>\n"
                         "<div id='External_subsection'>\n"
                         "<h2>External subsection</h2>\n"
                         "</div>\n"
                         "</div>")

    def test_section_set_name(self):
        s = Section("Test section",name="awesome_section")
        self.assertEqual(s.title,"Test section")
        self.assertEqual(s.name,"awesome_section")
        self.assertEqual(s.level,2)
        self.assertEqual(s.html(),
                         "<div id='awesome_section'>\n"
                         "<h2>Test section</h2>\n"
                         "</div>")

    def test_section_set_level(self):
        s = Section("Test section",level=4)
        self.assertEqual(s.title,"Test section")
        self.assertEqual(s.name,"Test_section")
        self.assertEqual(s.level,4)
        self.assertEqual(s.html(),
                         "<div id='Test_section'>\n"
                         "<h4>Test section</h4>\n"
                         "</div>")

    def test_section_add_css_classes(self):
        s = Section("Test section")
        s.add_css_classes("section")
        self.assertEqual(s.title,"Test section")
        self.assertEqual(s.name,"Test_section")
        self.assertEqual(s.level,2)
        self.assertEqual(s.css_classes,["section",])
        self.assertEqual(s.html(),
                         "<div id='Test_section' class='section'>\n"
                         "<h2>Test section</h2>\n"
                         "</div>")
        s.add_css_classes("space-cadet","big")
        self.assertEqual(s.css_classes,["section",
                                        "space-cadet",
                                        "big",])
        self.assertEqual(s.html(),
                         "<div id='Test_section' "
                         "class='section space-cadet big'>\n"
                         "<h2>Test section</h2>\n"
                         "</div>")

class TestTable(unittest.TestCase):
    """
    Tests for the Table class
    """
    def test_basic_table_no_content(self):
        t = Table(('Key','Value'))
        self.assertEqual(t.nrows,0)
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "</table>")

    def test_basic_table(self):
        t = Table(('Key','Value'))
        idx = t.add_row(Key="Employee name",Value="John Doe")
        self.assertEqual(idx,0)
        self.assertEqual(t.nrows,1)
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")

    def test_table_with_column_aliases(self):
        t = Table(('name','date','result'),
                  name="Experiment",
                  date="Date run",
                  result="Final result")
        idx = t.add_row(name="Test setup",date="10/10/2017",result="Ok")
        self.assertEqual(idx,0)
        self.assertEqual(t.nrows,1)
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Experiment</th><th>Date run</th><th>Final result</th></tr>\n"
                         "<tr><td>Test setup</td><td>10/10/2017</td><td>Ok</td></tr>\n"
                         "</table>")

    def test_table_multiple_rows(self):
        t = Table(('Key','Value'))
        idx = t.add_row(Key="Name",Value="John Doe")
        self.assertEqual(idx,0)
        self.assertEqual(t.nrows,1)
        idx = t.add_row(Key="D.O.B",Value="12/06/1982")
        self.assertEqual(idx,1)
        self.assertEqual(t.nrows,2)
        idx = t.add_row(Key="Mobile",Value="+44 1726254")
        self.assertEqual(idx,2)
        self.assertEqual(t.nrows,3)
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "<tr><td>Mobile</td><td>+44 1726254</td></tr>\n"
                         "</table>")

    def test_table_append_columns(self):
        t = Table(('Key','Person1'))
        t.add_row(Key="Name",Person1="John Doe")
        t.add_row(Key="D.O.B",Person1="12/06/1982")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Person1</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "</table>")
        t.append_columns("Person2","Person3")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Person1</th><th>Person2</th><th>Person3</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n"
                         "</table>")

    def test_table_append_columns_with_aliases(self):
        t = Table(('Key','Person1'),
                  Person1="Person #1")
        t.add_row(Key="Name",Person1="John Doe")
        t.add_row(Key="D.O.B",Person1="12/06/1982")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Person #1</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "</table>")
        t.append_columns("Person2","Person3",
                         Person2="Person #2",
                         Person3="Person #3")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Person #1</th><th>Person #2</th><th>Person #3</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n"
                         "</table>")

    def test_table_append_columns_fails_with_multiple_names(self):
        t = Table(('Key','Person1'))
        t.add_row(Key="Name",Person1="John Doe")
        t.add_row(Key="D.O.B",Person1="12/06/1982")
        self.assertRaises(KeyError,
                          t.append_columns,
                          "Person1")

    def test_table_set_value(self):
        t = Table(('Key','Value'))
        t.add_row(Key="Name",Value="John Doe")
        t.add_row(Key="D.O.B",Value="12/06/1982")
        t.add_row(Key="Mobile",Value="+44 1726254")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "<tr><td>Mobile</td><td>+44 1726254</td></tr>\n"
                         "</table>")
        t.append_columns("Value2")
        t.set_value(0,"Value2","Jane Doe")
        t.set_value(1,"Value2","19/04/1979")
        t.set_value(2,"Value2","+44 1745262")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th><th>Value2</th></tr>\n"
                         "<tr><td>Name</td><td>John Doe</td><td>Jane Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td><td>19/04/1979</td></tr>\n"
                         "<tr><td>Mobile</td><td>+44 1726254</td><td>+44 1745262</td></tr>\n"
                         "</table>")

    def test_table_add_css_classes(self):
        t = Table(('Key','Value'))
        t.add_row(Key="Employee name",Value="John Doe")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")
        t.add_css_classes("employee-data")
        self.assertEqual(t.html(),
                         "<table class='employee-data'>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")
        t.add_css_classes("summary","sunny")
        self.assertEqual(t.html(),
                         "<table class='employee-data summary sunny'>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")

    def test_table_add_css_classes_to_column(self):
        t = Table(('Key','Value'))
        t.add_row(Key="Employee name",Value="John Doe")
        t.add_row(Key="D.O.B",Value="12/06/1982")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "<tr><td>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "</table>")
        t.add_css_classes("key-column",column="Key")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th class='key-column'>Key</th><th>Value</th></tr>\n"
                         "<tr><td class='key-column'>Employee name</td><td>John Doe</td></tr>\n"
                         "<tr><td class='key-column'>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "</table>")
        t.add_css_classes("summary","sunny",column="Key")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th class='key-column summary sunny'>Key</th><th>Value</th></tr>\n"
                         "<tr><td class='key-column summary sunny'>Employee name</td><td>John Doe</td></tr>\n"
                         "<tr><td class='key-column summary sunny'>D.O.B</td><td>12/06/1982</td></tr>\n"
                         "</table>")
        self.assertRaises(KeyError,
                          t.add_css_classes,
                          "blah",
                          column="doesntexist")

    def test_table_no_header(self):
        t = Table(('Key','Value'))
        t.add_row(Key="Employee name",Value="John Doe")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")
        t.no_header()
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")

    def test_table_render_with_css_id(self):
        t = Table(('Key','Value'))
        t.add_row(Key="Employee name",Value="John Doe")
        self.assertEqual(t.html(),
                         "<table>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")
        self.assertEqual(t.html(css_id="employee-tbl"),
                         "<table id='employee-tbl'>\n"
                         "<tr><th>Key</th><th>Value</th></tr>\n"
                         "<tr><td>Employee name</td><td>John Doe</td></tr>\n"
                         "</table>")

class TestList(unittest.TestCase):
    """
    Tests for the List class
    """
    def test_empty_list(self):
        lst = List()
        self.assertEqual(lst.html(),"")

    def test_unordered_list(self):
        lst = List()
        lst.add_item("First item")
        self.assertEqual(lst.html(),
                         "<ul>"
                         "<li>First item</li>"
                         "</ul>")
        lst.add_item("Second item")
        self.assertEqual(lst.html(),
                         "<ul>"
                         "<li>First item</li>"
                         "<li>Second item</li>"
                         "</ul>")

    def test_ordered_list(self):
        lst = List(ordered=True)
        lst.add_item("First item")
        self.assertEqual(lst.html(),
                         "<ol>"
                         "<li>First item</li>"
                         "</ol>")
        lst.add_item("Second item")
        self.assertEqual(lst.html(),
                         "<ol>"
                         "<li>First item</li>"
                         "<li>Second item</li>"
                         "</ol>")

    def test_nested_lists(self):
        inner_lst = List()
        inner_lst.add_item("First item")
        inner_lst.add_item("Second item")
        self.assertEqual(inner_lst.html(),
                         "<ul>"
                         "<li>First item</li>"
                         "<li>Second item</li>"
                         "</ul>")
        outer_lst = List()
        outer_lst.add_item(inner_lst)
        self.assertEqual(outer_lst.html(),
                         "<ul>"
                         "<li><ul>"
                         "<li>First item</li>"
                         "<li>Second item</li>"
                         "</ul></li>"
                         "</ul>")

    def test_list_set_name(self):
        lst = List(name="my-list")
        lst.add_item("First item")
        lst.add_item("Second item")
        self.assertEqual(lst.html(),
                         "<ul id='my-list'>"
                         "<li>First item</li>"
                         "<li>Second item</li>"
                         "</ul>")

class TestImg(unittest.TestCase):
    """
    Tests for the Img class
    """
    def test_img(self):
        img = Img('picture.png')
        self.assertEqual(img.html(),
                         "<img src='picture.png' />")

    def test_img_with_height_and_width(self):
        img = Img('picture.png',height=100,width=150)
        self.assertEqual(img.html(),
                         "<img src='picture.png' "
                         "height='100' width='150' />")

    def test_img_with_name(self):
        img = Img('picture.png',name="awesome-picture")
        self.assertEqual(img.html(),
                         "<img id='awesome-picture' "
                         "src='picture.png' />")

    def test_img_with_alt_text(self):
        img = Img('picture.png',alt="An awesome picture")
        self.assertEqual(img.html(),
                         "<img src='picture.png' "
                         "alt='An awesome picture' />")

    def test_img_with_href(self):
        img = Img('picture.png',href="http://awesome.pix.com/")
        self.assertEqual(img.html(),
                         "<a href='http://awesome.pix.com/'>"
                         "<img src='picture.png' /></a>")

    def test_img_with_title(self):
        img = Img('picture.png',title="An awesome picture")
        self.assertEqual(img.html(),
                         "<img src='picture.png' "
                         "title='An awesome picture' />")
        
class TestLink(unittest.TestCase):
    """
    Tests for the Link class
    """
    def test_link(self):
        ahref = Link("My report","report.html")
        self.assertEqual(ahref.href,"report.html")
        self.assertEqual(ahref.html(),
                         "<a href='report.html'>My report</a>")

    def test_link_no_target(self):
        ahref = Link("http://example.com/report.html")
        self.assertEqual(ahref.href,"http://example.com/report.html")
        self.assertEqual(ahref.html(),
                         "<a href='http://example.com/report.html'>"
                         "http://example.com/report.html</a>")

    def test_link_to_section(self):
        s = Section("Introduction",name='intro_section')
        ahref = Link("Introduction",s)
        self.assertEqual(ahref.href,"#intro_section")
        self.assertEqual(ahref.html(),
                         "<a href='#intro_section'>Introduction</a>")

class TestTarget(unittest.TestCase):
    """
    Tests for the Target class
    """
    def test_target(self):
        tgt = Target("my_target")
        self.assertEqual(tgt.name,"my_target")
        self.assertEqual(tgt.html(),"<a id='my_target' />")

    def test_target_with_link(self):
        tgt = Target("my_target")
        ahref = Link("My target",tgt)
        self.assertEqual(ahref.href,"#my_target")
        self.assertEqual(ahref.html(),
                         "<a href='#my_target'>My target</a>")

class TestPara(unittest.TestCase):
    """
    Tests for the Para class
    """
    def test_empty_para(self):
        p = Para()
        self.assertEqual(p.html(),"")

    def test_para(self):
        p = Para("some text")
        self.assertEqual(p.html(),"<p>some text</p>")
        p.add("that was added")
        self.assertEqual(p.html(),"<p>some text that was added</p>")
        p = Para("some text","that was added")
        self.assertEqual(p.html(),"<p>some text that was added</p>")
        p = Para()
        p.add("some text","that was added")
        self.assertEqual(p.html(),"<p>some text that was added</p>")

    def test_para_with_mixture_of_types(self):
        img = Img("picture.png")
        ahref= Link("http://example.com")
        p = Para("Beautiful picture:",img,"See more at",ahref)
        self.assertEqual(p.html(),
                         "<p>Beautiful picture: "
                         "<img src='picture.png' /> "
                         "See more at "
                         "<a href='http://example.com'>"
                         "http://example.com</a></p>")

    def test_para_with_css_classes(self):
        p = Para("some text",css_classes=('cls1','cls2'))
        self.assertEqual(p.html(),"<p class='cls1 cls2'>some text</p>")
        p.add_css_classes("cls3","cls4")
        self.assertEqual(p.html(),
                         "<p class='cls1 cls2 cls3 cls4'>some text</p>")

    def test_para_non_zero(self):
        p = Para()
        self.assertFalse(p)
        p.add("some text")
        self.assertTrue(p)

class TestWarningIcon(unittest.TestCase):
    """
    Tests for the WarningIcon class
    """
    def test_warningicon(self):
        w = WarningIcon()
        self.assertEqual(w.html(),
                         "<img src='data:image/png;base64,%s' height='25' width='25' />" % WARNING_ICON_BASE64)

    def test_warningicon_with_title(self):
        w = WarningIcon(title="This is a warning")
        self.assertEqual(w.html(),
                         "<img src='data:image/png;base64,%s' height='25' width='25' title='This is a warning' />" % WARNING_ICON_BASE64)

    def test_warningicon_change_size(self):
        w = WarningIcon(size=12)
        self.assertEqual(w.html(),
                         "<img src='data:image/png;base64,%s' height='12' width='12' />" % WARNING_ICON_BASE64)

class TestDocumentIcon(unittest.TestCase):
    """
    Tests for the DocumentIcon class
    """
    def test_documenticon(self):
        d = DocumentIcon()
        self.assertEqual(d.html(),
                         "<img src='data:image/png;base64,%s' height='35' width='35' />" % DOCUMENT_ICON_BASE64)

    def test_documenticon_with_title(self):
        d = DocumentIcon(title="This is a document")
        self.assertEqual(d.html(),
                         "<img src='data:image/png;base64,%s' height='35' width='35' title='This is a document' />" % DOCUMENT_ICON_BASE64)

    def test_documenticon_change_size(self):
        d = DocumentIcon(size=25)
        self.assertEqual(d.html(),
                         "<img src='data:image/png;base64,%s' height='25' width='25' />" % DOCUMENT_ICON_BASE64)

class TestSanitizeCssStringFunction(unittest.TestCase):
    """
    Tests for the sanitize_css_string function
    """
    def test_sanitize_css_string(self):
        self.assertEqual(sanitize_css_string("test"),"test")
        self.assertEqual(sanitize_css_string("test123"),"test123")
        self.assertEqual(sanitize_css_string("test 123"),"test_123")
        self.assertEqual(sanitize_css_string("test\t 123"),"test__123")
        self.assertEqual(sanitize_css_string("test123.gz"),"test123gz")
        self.assertEqual(sanitize_css_string("test@123.90.10.3"),"test12390103")
