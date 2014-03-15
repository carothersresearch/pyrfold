"""2014-02-10 WEV
This is a template script for processing a directory of hyak submissions
It makes the following assumptions:
- your directory should look something like this
>Experiment Super Directory
>>hyak_processing_script.py
>>Experiment_Sub_directory1
>>>node-1
>>>node-2
>>>output
>>Experiment_Sub_directory1
>>>node-1
>>>node-2
>>>output

In brief this script makes a file framework which enables parallel
computational processing of large kinefold siumulation submissions.
It acomplishes this by created two types of jobfiles:
- A node deletion script for every node file
   This deletes the node file completely (a task that takes a surprising
    amount of computational time)
- A device processing script
    Deletes most of the raw output of kinefold
    Crates a proc_output folder which contains all of the desired data in
        a pickled data structure
These job files contain a series of booleans and directories
that tell a general processing script what to do
"""
#Import Modules
import shutil
import os
import pyrfold.hyak.create as hyakcreate
import pyrfold.hyak.submission as sub

#Define the experiment sub directories to be processed
#MUST FULL IN
NAMEOFDIRECTORYLIST = ['ASBV-1']
NAMEOFROUND = 'round-1_sub'
EMAIL = 'wmvoje@uw.edu'
NAMEOFNEXTROUND = 'round-2'
REFERENCEFILE = 'processed_reference_summary.csv'
ADDITIONALROUNDOFPROCESS = True
JUSTSUMMARY = True
CUTOFFFREQ = 0.1
POLRATE = 20
DWELLTIME = 5
FIVEPOS = 60
THREEPOS = 40
NUMBEROFSIMULATIONS = 10

#Do the work
ROOT = os.getcwd()
PROCESSINGSCRIPTNAME = 'delete_nodes_process_timecourse.py'
NAMEOFPROCESSINGFOLDER = 'hyakprocessing' #all processing will be done this dir
PATHTOPROCESSINGSCRIPT = os.path.join(ROOT, 'pyrfold', 'hyak', 'scripts',
                         PROCESSINGSCRIPTNAME)
shutil.copy(PATHTOPROCESSINGSCRIPT, ROOT)
NUMBEROFCORES = int(raw_input("Number of cores (8, 12, 16): "))
NUMBEROFNODES = 1
PATHTOPARMS = hyakcreate.framework_shell(os.path.join(ROOT,
    NAMEOFPROCESSINGFOLDER))

#This assumes that directory contains a set of round files
for directoryname in NAMEOFDIRECTORYLIST:
    PARTPATH = os.path.join(ROOT, directoryname)
    #List the directors
    rounddirectorylist = [ls for ls in os.listdir(PARTPATH)
                            if os.path.isdir(os.path.join(PARTPATH, ls))]
    for roundpath in rounddirectorylist:
        if NAMEOFROUND not in roundpath:
            continue
        EXPPATH = os.path.join(PARTPATH, roundpath)
        #Assumes that the node and output directories are immediately below this
        dirls = [ls for ls in os.listdir(EXPPATH)
                            if os.path.isdir(os.path.join(EXPPATH, ls))]
        for dirname in dirls:
            NODEOUTPATH = os.path.join(EXPPATH, dirname)
            if 'node' in dirname:
                #write node job files
                with open(os.path.join(PATHTOPARMS, 'job-' + directoryname +
                                           dirname + '.txt'), 'wb') as txtfile:
                    txtfile.write('deletenodes')
                    txtfile.write('\n')
                    txtfile.write(NODEOUTPATH)
                    txtfile.write('\n')
            if 'output' in dirname:
                simulationdirs = [ls for ls in os.listdir(NODEOUTPATH)
                    if os.path.isdir(os.path.join(NODEOUTPATH, ls))]
                for sim in simulationdirs:
                    SIMDIR = os.path.join(NODEOUTPATH, sim)
                    #Now to build the file
                    with open(os.path.join(PATHTOPARMS, 'job-' + sim + '.txt')
                                                             , 'wb') as txtfile:
                        txtfile.write('timecourse')
                        txtfile.write('\n')
                        #Final Structure
                        txtfile.write('1')
                        txtfile.write('\n')
                        #Single directory(as opposed to all output)
                        txtfile.write('1')
                        txtfile.write('\n')
                        #Folder to process
                        txtfile.write(SIMDIR)
                        txtfile.write('\n')
                        #path to the output
                        txtfile.write(os.path.join(EXPPATH, 'proc_output'))
                        txtfile.write('\n')
                        #Get the reference files
                        txtfile.write(os.path.join(EXPPATH, 'sub_summary.csv'))
                        txtfile.write('\n')

        if ADDITIONALROUNDOFPROCESS:
            with open(os.path.join(PATHTOPARMS, 'finaljob-' + directoryname
                                            + '.txt'), 'wb') as txtfile:
                txtfile.write('additionalround')
                txtfile.write('\n')
                #Just processe or not
                if JUSTSUMMARY:
                    txtfile.write('1')
                else:
                    txtfile.write('0')
                txtfile.write('\n')
                #Processed Directory
                txtfile.write(os.path.join(EXPPATH, 'proc_output'))
                txtfile.write('\n')
                #Path to summary files directory
                txtfile.write(os.path.join(PARTPATH, 'summary'))
                txtfile.write('\n')
                #Path to reference part summary .csv file
                txtfile.write(os.path.join(ROOT, REFERENCEFILE))
                txtfile.write('\n')
                #path (to be made) that holds the additional round of processing
                txtfile.write(os.path.join(ROOT, PARTPATH, NAMEOFNEXTROUND))
                txtfile.write('\n')
                #Name of previous round
                txtfile.write(NAMEOFROUND)
                txtfile.write('\n')
                #path to currentround submission file:
                txtfile.write(os.path.join(ROOT, PARTPATH, NAMEOFROUND +
                                                                      '.csv'))
                txtfile.write('\n')
                #Cutoff frequency
                txtfile.write(str(CUTOFFFREQ))
                txtfile.write('\n')
                #Next positions of the simulation
                #POLRATE
                txtfile.write(str(POLRATE))
                txtfile.write('\n')
                #DWELLTIME
                txtfile.write(str(DWELLTIME))
                txtfile.write('\n')
                #FIVEPOS
                txtfile.write(str(FIVEPOS))
                txtfile.write('\n')
                #THREEPOS
                txtfile.write(str(THREEPOS))
                txtfile.write('\n')
                #Number of simulations
                txtfile.write(str(NUMBEROFSIMULATIONS))
                txtfile.write('\n')
                #Email to send the updates to
                txtfile.write(EMAIL)
                txtfile.write('\n')


CALLCOMMAND = 'python ' + os.path.join(ROOT, PROCESSINGSCRIPTNAME)
os.chdir(NAMEOFPROCESSINGFOLDER)
hyakcreate.mybundle_sub(EMAIL, NUMBEROFCORES, 1,
                    '15:00:00', 'highthroughputproc', finaljob=True)
hyakcreate.general_myscript_sub(CALLCOMMAND)
ROOT = os.getcwd()
hyakcreate.symlinks('MyScript.sh', os.path.join(ROOT,'myscript-parms'),
                 os.path.join(ROOT,'myscript-links'), '*.txt')

#Submit the file
sub.submit_basic_hyak_framework(os.path.join(ROOT,'myscript-links'))
