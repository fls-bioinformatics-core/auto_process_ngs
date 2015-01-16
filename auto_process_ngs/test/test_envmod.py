#######################################################################
# Tests for envmod.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
import auto_process_ngs.envmod as envmod

class TestEnvMod(unittest.TestCase):
    """Tests for the envmod module

    """

    def setUp(self):
        # Create set of fake modulefiles
        self.modulesdir = tempfile.mkdtemp()
        try:
            os.environ['MODULEPATH'] += ":%s" % self.modulesdir
            for mod in ('app1','app2'):
                filen = os.path.join(self.modulesdir,mod)
                open(filen,'w').write('#%%Module1.0\nprepend-path PATH /opt/%s\n' % mod)
        except KeyError:
            shutil.rmtree(self.wd)
            self.fail("Unable to set up test environment")

    def tearDown(self):
        # Remove the temporary test directory
        shutil.rmtree(self.wd)

    def test_loaded(self):
        """envmod.loaded function returns list of loaded modules
        """
        self.assertEqual(envmod.loaded(),[])

    def test_load(self):
        """envmod.load function loads specified module
        """
        envmod.load('app1')
        self.assertEqual(envmod.loaded(),['app1'])
        envmod.load('app2')
        self.assertEqual(envmod.loaded(),['app1','app2'])

    def test_load_missing_module(self):
        """envmod.load function raises exception if module is missing
        """
        self.assertRaises(Exception,envmod.load,'app3')

    def test_unload(self):
        """envmod.unload function unloads specified module
        """
        envmod.load('app1','app2')
        self.assertEqual(envmod.loaded(),['app1','app2'])
        envmod.unload('app1')
        self.assertEqual(envmod.loaded(),['app2'])

    def test_unload_all(self):
        """envmod.unload_all function unloads all modules
        """
        envmod.load('app1','app2')
        self.assertEqual(envmod.loaded(),['app1','app2'])
        envmod.unload_all()
        self.assertEqual(envmod.loaded(),[])

