"""  WEV
This script will go through a number of rounds pull out the best
performers and resubmit them through a final trial
"""
#Import Modules
import shutil
import os
from collections import Counter
import pyrfold.hyak.create as hyakcreate
import pyrfold.hyak.submission as hyaks
import csv

#Define the experiment sub directories to be processed
#MUST FULL IN
NAMEOFDIRECTORYLIST = ['ASBV-1', 'ASBV-3', 'sTRSV-2', 'Sman']
REFERENCEFILE = 'processed_reference_summary.csv'
NAMEOFNEXTROUND = 'final-round'
PREVIOUSUBFORM = 'round-1_sub.csv'
POLRATES = [10, 20, 40]
FIVEPOS = [50, 60]
THREEPOS = [30, 25, 40]
DWELLTIME = [5]
LISTSOFCONDITIONS = []
for polrate in POLRATES:
    for five in FIVEPOS:
        for three in THREEPOS:
            for dwell in DWELLTIME:
                LISTSOFCONDITIONS.append([polrate, dwell, five, three])

EMAIL = 'wmvoje@uw.edu'
NUMBEROFSIMULATIONS = 50
NUMBEROFDEVICESTOTEST = 20

#Do the work
ROOT = os.getcwd()
PROCESSINGSCRIPTNAME = 'delete_nodes_process_timecourse.py'

for directoryname in NAMEOFDIRECTORYLIST:
    PARTPATH = os.path.join(ROOT, directoryname)
    pathtofirstsub = os.path.join(PARTPATH, PREVIOUSUBFORM)
    pathtonextsub = os.path.join(PARTPATH, NAMEOFNEXTROUND)
    #List the directors
    directorylist = [ls for ls in os.listdir(PARTPATH)
                            if os.path.isdir(os.path.join(PARTPATH, ls))]
    for dirpath in directorylist:
        #Find the summary directory
        if 'summary' not in dirpath:
            continue
        SUMPATH = os.path.join(PARTPATH, dirpath)
        #Now we have the summary path
        summaryfiles = [ls for ls in os.listdir(SUMPATH)
                            if os.path.isfile(os.path.join(SUMPATH, ls))]
        dictofsummaryfiles = {}
        for sumfilename in summaryfiles:
            if 'round' not in sumfilename:
                continue
            tempnum = int(sumfilename.split('-')[1].split('_')[0])
            dictofsummaryfiles[tempnum] = os.path.join(SUMPATH, sumfilename)
        tempkeys = dictofsummaryfiles.keys()
        tempkeys.sort(reverse=True)
        dictofwinners = Counter()
        for key in tempkeys:
            with open(dictofsummaryfiles[key], 'r') as csvfile:
                reader = csv.reader(csvfile)
                next(reader, None)
                for row in reader:
                    if row[0] in dictofwinners:
                        dictofwinners[row[0]] += float(row[1])
                    else:
                        dictofwinners[row[0]] = float(row[1])
        #From this it's easy to select the winners
        dicttosubmit = {}
        for winner in dictofwinners.most_common(NUMBEROFDEVICESTOTEST):
            dicttosubmit[winner[0]] = winner[1]
        hyaks.additional_round_submission(dicttosubmit, pathtofirstsub,
            pathtonextsub, LISTSOFCONDITIONS, [EMAIL, NUMBEROFSIMULATIONS])
