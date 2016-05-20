"""
importing all of the sub modules
"""
import hyak
import compare
import pyrfile
import design
import fold
import analyze
import utilities
from foldingsub import FoldingSubData

try:
    from pyrfold._version import __version__
except ImportError:
    pass


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
