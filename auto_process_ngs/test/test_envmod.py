#######################################################################
# Tests for envmod.py module
#######################################################################

import unittest
import tempfile
import os
import shutil
import auto_process_ngs.envmod as envmod

def _setup_example_modulefiles(*apps):
    # Create set of fake modulefiles
    modulesdir = tempfile.mkdtemp()
    try:
        os.environ['MODULEPATH'] += ":%s" % modulesdir
        for mod in apps:
            filen = os.path.join(modulesdir,mod)
            open(filen,'w').write('#%%Module1.0\nprepend-path PATH /opt/%s\n' % mod)
    except KeyError:
        shutil.rmtree(modulesdir)
        raise Exception("No MODULEPATH environment variable")
    # Wipe existing environment
    try:
        os.environ["LOADEDMODULES"] = ''
    except KeyError:
        pass
    return modulesdir

def _teardown_example_modulefiles(modulesdir):
    # Remove the temporary test directory
    shutil.rmtree(modulesdir)

class TestEnvModLoaded(unittest.TestCase):
    """Tests for the envmod module

    """

    def setUp(self):
        self.modulesdir = _setup_example_modulefiles('app1','app2')

    def tearDown(self):
        _teardown_example_modulefiles(self.modulesdir)

    def test_loaded(self):
        """envmod.loaded function returns list of loaded modules
        """
        self.assertEqual(envmod.loaded(),[])

class TestEnvModLoad(unittest.TestCase):
    """Tests for the envmod module

    """

    def setUp(self):
        self.modulesdir = _setup_example_modulefiles('app1','app2')

    def tearDown(self):
        _teardown_example_modulefiles(self.modulesdir)

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

class TestEnvModUnload(unittest.TestCase):
    """Tests for the envmod module

    """

    def setUp(self):
        self.modulesdir = _setup_example_modulefiles('app1','app2')

    def tearDown(self):
        _teardown_example_modulefiles(self.modulesdir)

    def test_unload(self):
        """envmod.unload function unloads specified module
        """
        envmod.load('app1','app2')
        loaded = envmod.loaded()
        loaded.sort()
        self.assertEqual(loaded,['app1','app2'])
        envmod.unload('app1')
        loaded = envmod.loaded()
        loaded.sort()
        self.assertEqual(loaded,['app2'])

class TestEnvModUnloadAll(unittest.TestCase):
    """Tests for the envmod module

    """

    def setUp(self):
        self.modulesdir = _setup_example_modulefiles('app1','app2')

    def tearDown(self):
        _teardown_example_modulefiles(self.modulesdir)

    def test_unload_all(self):
        """envmod.unload_all function unloads all modules
        """
        envmod.load('app1','app2')
        self.assertEqual(envmod.loaded(),['app1','app2'])
        envmod.unload_all()
        self.assertEqual(envmod.loaded(),[])

