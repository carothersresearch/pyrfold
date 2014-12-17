"""
importing all of the sub modules
"""
from . import hyak
from . import compare
from . import pyrfile
from . import design
from . import fold
from . import analyze
from . foldingsub import FoldingSubData

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
