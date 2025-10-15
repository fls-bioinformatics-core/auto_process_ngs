# In order to run tests without installing the auto_process_ngs package
# we need to add the package "bin" directory

import os
__AUTO_PROCESS_BIN__ = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "bin"))
os.environ["PATH"] = f"{__AUTO_PROCESS_BIN__}{os.pathsep}{os.environ['PATH']}"