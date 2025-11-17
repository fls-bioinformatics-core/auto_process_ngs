#######################################################################
# Tests for docs.py module
#######################################################################

from auto_process_ngs.docs import RstSimpleTab
from auto_process_ngs.docs import RstGridTablele
import unittest


class TestRstSimpleTable(unittest.TestCase):
    """
    Tests for RstSimpleTable class.
    """

    def test_construct_table_basic(self):
        """
        RstSimpleTable: construct basic table
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        header = ['Column 1', 'Column 2', 'Column 3']
        table = RstSimpleTable(table_data)
        expected_output = [
            '=========  =========  =========',
            'Column 1   Column 2   Column 3 ',
            '=========  =========  =========',
            'col1_row1  col2_row1  col3_row1',
            'col1_row2  col2_row2  col3_row2',
            'col1_row3  col2_row3  col3_row3',
            '=========  =========  =========',
        ]
        output = list(table.construct_table(header=header))
        self.assertEqual(output, expected_output)

    def test_construct_table_no_header(self):
        """
        RstSimpleTable: construct table without header
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        table = RstSimpleTable(table_data)
        expected_output = [
            '=========  =========  =========',
            'col1_row1  col2_row1  col3_row1',
            'col1_row2  col2_row2  col3_row2',
            'col1_row3  col2_row3  col3_row3',
            '=========  =========  =========',
        ]
        output = list(table.construct_table())
        self.assertEqual(output, expected_output)

    def test_construct_table_with_indent(self):
        """
        RstSimpleTable: make table with indentation
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        header = ['Column 1', 'Column 2', 'Column 3']
        table = RstSimpleTable(table_data)
        expected_output = [
            '   =========  =========  =========',
            '   Column 1   Column 2   Column 3 ',
            '   =========  =========  =========',
            '   col1_row1  col2_row1  col3_row1',
            '   col1_row2  col2_row2  col3_row2',
            '   col1_row3  col2_row3  col3_row3',
            '   =========  =========  =========',
        ]
        output = list(table.construct_table(header=header, indent="   "))
        self.assertEqual(output, expected_output)


class TestRstGridTable(unittest.TestCase):
    """
    Tests for RstGridTable class.
    """

    def test_construct_table_basic(self):
        """
        RstTable: construct basic table
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        header = ['Column 1', 'Column 2', 'Column 3']
        table = RstGridTable(table_data)
        expected_output = [
            '+-----------+-----------+-----------+',
            '| Column 1  | Column 2  | Column 3  |',
            '+===========+===========+===========+',
            '| col1_row1 | col2_row1 | col3_row1 |',
            '+-----------+-----------+-----------+',
            '| col1_row2 | col2_row2 | col3_row2 |',
            '+-----------+-----------+-----------+',
            '| col1_row3 | col2_row3 | col3_row3 |',
            '+-----------+-----------+-----------+',
        ]
        output = list(table.construct_table(header=header))
        self.assertEqual(output, expected_output)

    def test_construct_table_merge(self):
        """
        RstGridTable: construct table construction with merged cells
        """
        table_data = [
            ['block1', 'col2_row1', 'col3_row1'],
            ['block1', 'col2_row2', 'col3_row2'],
            ['block2', 'col2_row3', 'col3_row3'],
        ]
        header = ['Block', 'Column 2', 'Column 3']
        table = RstGridTable(table_data)
        expected_output = [
            '+--------+-----------+-----------+',
            '| Block  | Column 2  | Column 3  |',
            '+========+===========+===========+',
            '| block1 | col2_row1 | col3_row1 |',
            '|        +-----------+-----------+',
            '|        | col2_row2 | col3_row2 |',
            '+--------+-----------+-----------+',
            '| block2 | col2_row3 | col3_row3 |',
            '+--------+-----------+-----------+',
        ]
        output = list(table.construct_table(header=header))
        self.assertEqual(output, expected_output)

    def test_construct_table_no_header(self):
        """
        RstGridTable: construct table without header
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        table = RstGridTable(table_data)
        expected_output = [
            '+-----------+-----------+-----------+',
            '| col1_row1 | col2_row1 | col3_row1 |',
            '+-----------+-----------+-----------+',
            '| col1_row2 | col2_row2 | col3_row2 |',
            '+-----------+-----------+-----------+',
            '| col1_row3 | col2_row3 | col3_row3 |',
            '+-----------+-----------+-----------+',
        ]
        output = list(table.construct_table())
        self.assertEqual(output, expected_output)

    def test_construct_table_with_indent(self):
        """
        RstGridTable: make table with indentation
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        header = ['Column 1', 'Column 2', 'Column 3']
        table = RstGridTable(table_data)
        expected_output = [
            '   +-----------+-----------+-----------+',
            '   | Column 1  | Column 2  | Column 3  |',
            '   +===========+===========+===========+',
            '   | col1_row1 | col2_row1 | col3_row1 |',
            '   +-----------+-----------+-----------+',
            '   | col1_row2 | col2_row2 | col3_row2 |',
            '   +-----------+-----------+-----------+',
            '   | col1_row3 | col2_row3 | col3_row3 |',
            '   +-----------+-----------+-----------+',
        ]
        output = list(table.construct_table(header=header, indent="   "))
        self.assertEqual(output, expected_output)


