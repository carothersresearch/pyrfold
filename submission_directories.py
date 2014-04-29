"""2014-02-10 11:15 WEV
This is a framework for a submission form.

It makes the following assumptions:
your submission .csv file exists within a directory that also contains
a copy of pyrfold/
"""
import os
from pyrfold.hyak import submission

#Need to fill in the list of submission files
LISTOFDIRECTORIES = ['ASBV-1', 'ASBV-3', 'Sman', 'sTRSV-3']
subfilename = 'round-1_sub.csv'
EMAIL = 'wmvoje@uw.edu'
NUMBEROFNODES = '2'

ROOT = os.getcwd()
for direct in LISTOFDIRECTORIES:
    temppath = os.path.join(ROOT, direct, subfilename)
    SUBROOT = os.path.join(ROOT, direct)
    submission.submit_file(temppath, SUBROOT, EMAIL, nodes=NUMBEROFNODES)

