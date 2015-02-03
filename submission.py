"""2014-02-10 11:15 WEV
This is a framework for a submission form.

It makes the following assumptions:
your submission .csv file exists within a directory that also contains
a copy of pyrfold/
"""
import os
from pyrfold.hyak import submission
import shutil

LISTOFSUBMISSIONFILES = []
EMAIL = ''
NUMBEROFNODES = 2


ROOT = os.getcwd()
PROCESSINGSCRIPTNAME = 'kinefold_simulation_and_processing.py'
PATHTOPROCESSINGSCRIPT = os.path.join(ROOT, 'pyrfold', 'hyak', 'scripts',
                                      PROCESSINGSCRIPTNAME)
# Copy the processing script
shutil.copy(PATHTOPROCESSINGSCRIPT, ROOT)
PATHTOPROCESSINGSCRIPT = os.path.join(ROOT, PROCESSINGSCRIPTNAME)
# Need to fill in the list of submission files

ROOT = os.getcwd()
for submissionfile in LISTOFSUBMISSIONFILES:
    submission.submit_file(os.path.join(ROOT, submissionfile), ROOT, EMAIL,
                           processing_script_path=PATHTOPROCESSINGSCRIPT,
                           cores=12, backfill=True, nameofexperiment='auto',
                           nodes=NUMBEROFNODES)
