"""
This is a framework for a submission form.

It makes the following assumptions:
your submission .csv file exists within a directory that also contains
a copy of pyrfold/

So you dir should look like this
>dir
>>pyrfold
>>submissionname.csv
"""

import os
from pyrfold.hyak import submission
import shutil

# MUST FILL IN!
# Should be in the form of  ['*.csv', *.csv']
LISTOFSUBMISSIONFILES = []
# Email address for hyak updates
EMAIL = ''

# HYAK LOGISTICS
# Backfill
BACKFILL = True
# Can use 8, 12, or 16
# NOTE if NOT using backfill must be 16
CORES = 12
# Number of nodes to subimt to
# NOTE if NOT using backfill we have 3 nodes total
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
                           cores=CORES, backfill=BACKFILL,
                           nameofexperiment='auto', nodes=NUMBEROFNODES)
