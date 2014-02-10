"""importing all of the sub modules
"""
from . import convert
from . import hyak
from . import kineprocess
from . import process
from . import pyrfile
from . import design
from . import submit

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
