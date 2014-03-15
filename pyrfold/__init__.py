"""importing all of the sub modules
"""
from . import hyak
from . import compare
from . import pyrfile
from . import design

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class KineSubData:
    def __init__(self):
        self.sequence =""
        self.windowstart = 0
        self.windowstop = 0
        self.partstartstoplist = []
        self.partnamelist = []
        self.polrate = 30
        self.foldtimeafter = 1
        self.experimenttype = 2
        self.pseudoknots = 0
        self.entanglements = 0
