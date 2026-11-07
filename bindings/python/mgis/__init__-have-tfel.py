import os
import warnings

if os.name == 'nt' and hasattr(os, 'add_dll_directory'):
    # since Python 3.8, the directories listed in PATH are not searched for
    # the DLLs needed by the extension modules: the bin directory of MGIS,
    # where MFrontGenericInterface.dll is installed, is registered
    if "MGISHOME" in os.environ:
        mgishome = os.environ["MGISHOME"]
    else:
        # remove \bin\python...\site-packages\mgis\__init__.py
        mgishome = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
    if os.path.isdir(os.path.join(mgishome, "bin")):
        os.add_dll_directory(os.path.join(mgishome, "bin"))

try:
    import mfront.database
except ImportError:
    warnings.warn("loading mfront.database model failed")

from ._mgis import *
