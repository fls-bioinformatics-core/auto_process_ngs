#######################################################################
# Tests for docs.py module
#######################################################################

from auto_process_ngs.docs import RstTable
import unittest


class TestRstTable(unittest.TestCase):
    """
    Tests for RstTable class.
    """

    def test_construct_table_basic(self):
        """
        Test basic table construction.
        """
        table_data = [
            ['col1_row1', 'col2_row1', 'col3_row1'],
            ['col1_row2', 'col2_row2', 'col3_row2'],
            ['col1_row3', 'col2_row3', 'col3_row3'],
        ]
        header = ['Column 1', 'Column 2', 'Column 3']
        table = RstTable(table_data)
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
        Test table construction with merged cells.
        """
        table_data = [
            ['block1', 'col2_row1', 'col3_row1'],
            ['block1', 'col2_row2', 'col3_row2'],
            ['block2', 'col2_row3', 'col3_row3'],
        ]
        header = ['Block', 'Column 2', 'Column 3']
        table = RstTable(table_data)
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

