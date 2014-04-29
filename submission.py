"""2014-02-10 11:15 WEV
This is a framework for a submission form.

It makes the following assumptions:
your submission .csv file exists within a directory that also contains
a copy of pyrfold/
"""
import os
from pyrfold.hyak import submission

#Need to fill in the list of submission files
LISTOFSUBMISSIONFILES = []
EMAIL = ''
NUMBEROFNODES = 2

ROOT = os.getcwd()
for submissionfile in LISTOFSUBMISSIONFILES:
    submission.submit_file(os.path.join(ROOT, submissionfile), ROOT, EMAIL,
        cores=12, backfill=True, nameofexperiment='auto', nodes=NUMBEROFNODES)
