"""
This is a python script that is designed to operate with the hyak
framework for parallel processing
"""
import shutil
import fileinput
import os
import pyrfold.hyak.process as hyakp
import pyrfold.hyak.submission as hyaks
import pyrfold.pyrfile as pyrfile
import pyrfold.hyak
import cPickle as pickle
from pyrfold import compare
import csv

#from pyrfold import hyak.process as haykp

############################################################
################# Grabbing data from input  ################
############################################################
INPUTFILE = fileinput.input()
PROCESSINGTYPE = INPUTFILE.next().strip()

if PROCESSINGTYPE == 'deletenodes':
    #Gather the rest of the information
    NODEDIRECTORY = INPUTFILE.next().strip()
    shutil.rmtree(NODEDIRECTORY)
if PROCESSINGTYPE == 'timecourse':
    #Pull all of the boolean decision variables
    FINALSTRUCTURES = int(INPUTFILE.next().strip())
    SINGLEDIRECTORY = int(INPUTFILE.next().strip())
    #ADDITIONALROUNDOFPROCESS = int(INPUTFILE.next().strip())
    #Process the data to get timecourse data
    FOLDERTOPROCESS = INPUTFILE.next().strip()
    PROCESSEDDIRECTORY = INPUTFILE.next().strip()
    NAMEOFFOLDERTOPROCESS = os.path.basename(FOLDERTOPROCESS)
    #Do processing specific based on the single or multiple directory thing
    if SINGLEDIRECTORY:
        #create the timecourse files
        hyakp.timecourse(FOLDERTOPROCESS, PROCESSEDDIRECTORY,
            singledirectory=True)
        #remove unwanted files
        hyakp.clean_up_output_data(FOLDERTOPROCESS, singlefile=True)
    else:
        #create the timecourse files
        hyakp.timecourse(FOLDERTOPROCESS, PROCESSEDDIRECTORY)
        #remove unwanted files
        hyakp.clean_up_output_data(FOLDERTOPROCESS)
    if FINALSTRUCTURES:
        #Consolidate the data into a final structure file
        ##Get the reference files
        SUBSUMMARYFILE = INPUTFILE.next().strip()
        SUBSUMMARYDATA = pyrfile.load_pickled_sub_summary(SUBSUMMARYFILE)
        ##process to single files
        NAMEOFFOLDER = os.path.basename(FOLDERTOPROCESS)
        if SINGLEDIRECTORY:
            hyakp.finalstructure(PROCESSEDDIRECTORY, SUBSUMMARYDATA,
                                    singledirectoryname=NAMEOFFOLDER + '.p')
        else:
            #NEED TO BUILD
            pass
if PROCESSINGTYPE == 'additionalround':
    #This will be linked to a job file that is submitted after
    #all other files have been submitted
    #just summary data
    JUSTSUMMARY = int(INPUTFILE.next().strip())
    #Path to timecourse and final data:
    PROCESSEDDIRECTORY = INPUTFILE.next().strip()
    #Path to summary files directory
    SUMMARYPATH = INPUTFILE.next().strip()
    #Path to reference part summary .csv file
    REFERENCESUMMARYFILE = INPUTFILE.next().strip()
    #path (to be made) that holds the additional round of processing
    NEXTROUNDPATH = INPUTFILE.next().strip()
    #Name of previous round
    CURRENTROUNDNAME = INPUTFILE.next().strip()
    #path to currentround submission file:
    CURRENTSUBMISSIONPATH = INPUTFILE.next().strip()
    #Cutoff frequency
    FOLDINGCUTOFF = float(INPUTFILE.next().strip())
    #Next positions of the simulation
    POLRATE = float(INPUTFILE.next().strip())
    DWELLTIME = float(INPUTFILE.next().strip())
    FIVEPOS = INPUTFILE.next().strip()
    THREEPOS = INPUTFILE.next().strip()
    if FIVEPOS == 'None':
        FIVEPOS = None
    else:
        FIVEPOS = int(FIVEPOS)
    if THREEPOS == 'None':
        THREEPOS = None
    else:
        THREEPOS = int(THREEPOS)
    #Number of simulations
    NUMBEROFSIMULATIONS = int(INPUTFILE.next().strip())
    EMAIL = INPUTFILE.next().strip()
    ####################################################################
    # Make the comparisions
    ####################################################################
    #Get the reference structures to make comaprisons to
    REFSTRUCT = pyrfile.summary_refpart_to_sequence(REFERENCESUMMARYFILE)
    #Grab all of the final structure data
    PATHTOFINALSTRUCTS = os.path.join(PROCESSEDDIRECTORY, 'finalstructure')
    #Grab all of the .p files for iteration
    PICKLELIST = [ls for ls in os.listdir(PATHTOFINALSTRUCTS) if
                    os.path.isfile(os.path.join(PATHTOFINALSTRUCTS, ls))]
    #Make the comparison for all of the devices
    STATDICT = {}
    for pickfile in PICKLELIST:
        pickpath = os.path.join(PATHTOFINALSTRUCTS, pickfile)
        with open(pickpath, 'r') as pick:
            tempfinaldict = pickle.load(pick)
            devicename = pickfile.split('.')[0]
            # devicename = devicename.split('.')[0]
            #Make the comparison
            tempcomparedict = compare.folding_frequency(REFSTRUCT,
                                                                tempfinaldict)
            STATDICT[devicename] = tempcomparedict

    #Build the roundsummary folder
    if not os.path.exists(SUMMARYPATH):
        os.mkdir(SUMMARYPATH)
    #Build the experimental condition.csv and the round-summary
    PATHTOROUND = os.path.join(SUMMARYPATH, CURRENTROUNDNAME + '.csv')
    with open(PATHTOROUND, 'wb') as roundfile:
        writer = csv.writer(roundfile)
        #build header
        HEADER = ['device_name']
        HEADER.extend(STATDICT[STATDICT.keys()[0]].keys())
        writer.writerow(HEADER)
        for device in STATDICT:
            linetowrite = [device]
            for part in STATDICT[device]:
                linetowrite.append(STATDICT[device][part])
            writer.writerow(linetowrite)

    PATHTOEXPCOND = os.path.join(SUMMARYPATH, 'exp_cond_summary.csv')
    if os.path.isfile(PATHTOEXPCOND):
        #open this file and append the information
        with open(PATHTOEXPCOND, 'a') as expcondfile:
            writer = csv.writer(expcondfile)
            linetowrite = [CURRENTROUNDNAME]
            linetowrite.extend(pyrfile.sub_file(CURRENTSUBMISSIONPATH,
                                             justexperimentalconditions=True))
            writer.writerow(linetowrite)
    else:
        with open(PATHTOEXPCOND, 'w') as expcondfile:
            writer = csv.writer(expcondfile)
            HEADER = ['roundname', 'polrate', 'dwelltime', 'fiveprime',
                                                                'threeprime']
            writer.writerow(HEADER)
            linetowrite = [CURRENTROUNDNAME]
            linetowrite.extend(pyrfile.sub_file(CURRENTSUBMISSIONPATH,
                                             justexperimentalconditions=True))
            writer.writerow(linetowrite)
    #Write hyakp.get_exp_conditions
    #BUILD OTHER SUBMISSION
    hyaks.additional_round_submission(STATDICT, CURRENTSUBMISSIONPATH,
         FOLDINGCUTOFF, NEXTROUNDPATH, [POLRATE, DWELLTIME, FIVEPOS, THREEPOS],
         [EMAIL, NUMBEROFSIMULATIONS])














