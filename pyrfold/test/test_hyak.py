"""
Testing package specifically for all hyak functions
"""
import unittest
import os
from pyrfold.hyak import create
from pyrfold import pyrfile

class TestCreateFunctions(unittest.TestCase):

    def setUp(self):
        # Need to load in a test case
        filepath = os.path.join(os.path.dirname(__file__),
                                                 'testdata/submissiondata.csv')
        self.devicenametosubobj = pyrfile.sub_file(filepath)
        self.nodes = [1, 2, 3, 4, 5, 6, 7, 'auto']

    def test_nodes_to_devivces(self):
        #Now confirm that all of the devicies are in nodestodevices
        numberofdevices = len(self.devicenametosubobj.keys())
        deviceset = set(self.devicenametosubobj.keys())
        for nodetotest in self.nodes:
            nodestodevices = create.nodes_to_devivces(
            self.devicenametosubobj, nodetotest)
            countofdevices = 0
            for node in nodestodevices:
                countofdevices += len(nodestodevices[node])
            self.assertEqual(countofdevices, numberofdevices)
        for nodestotest in self.nodes:
            nodestodevices = create.nodes_to_devivces(
            self.devicenametosubobj, nodetotest)
            for node in nodestodevices:
                tempdeviceslist = nodestodevices[node]
                tempdevicesset = set(tempdeviceslist)
                self.assertEqual(len(tempdeviceslist), len(tempdevicesset))

if __name__ == '__main__':
    unittest.main()
