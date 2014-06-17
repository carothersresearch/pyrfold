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
########################################################################
########################################################################
####################  FILL IN ##########################################

NAMEOFDIRECTORYLIST = []
EMAIL = ''
FINALSTRUCTURE = 1

### COMPRESSION stuff ###
# ENTER 1 or 0 (yes or no)
COMPRESSOUTPUT = 1
COMPRESSPROCESSEDOUTPUT = 1

####################  FILL IN ##########################################
########################################################################
########################################################################

#####Script does the Remainder of the work FOR YOU!!!!!!!!!!!!!!!!
ROOT = os.getcwd()
PROCESSINGSCRIPTNAME = 'delete_nodes_process_timecourse.py'
NAMEOFPROCESSINGFOLDER = 'hyakprocessing' #all processing will be done this dir
PATHTOPROCESSINGSCRIPT = os.path.join(ROOT, 'pyrfold', 'hyak', 'scripts',
                         PROCESSINGSCRIPTNAME)
shutil.copy(PATHTOPROCESSINGSCRIPT, ROOT)
NUMBEROFCORES = 16 #int(raw_input("Number of cores (8, 12, 16): "))
NUMBEROFNODES = 1
PATHTOPARMS = hyakcreate.framework_shell(os.path.join(ROOT,
    NAMEOFPROCESSINGFOLDER))

#This assumes that the directory contains only node
for directoryname in NAMEOFDIRECTORYLIST:
    #Assumes that the node and output directories are immediately below this
    EXPPATH = os.path.join(ROOT, directoryname)
    dirls = [ls for ls in os.listdir(EXPPATH)
                        if os.path.isdir(os.path.join(EXPPATH, ls))]
    for dirname in dirls:
        NODEOUTPATH = os.path.join(EXPPATH, dirname)
        if 'node' in dirname:
            #write node job files
            with open(os.path.join(PATHTOPARMS, 'job-' + dirname + '.txt')
                                                        , 'wb') as txtfile:
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
                    txtfile.write(str(FINALSTRUCTURE))
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
                    txtfile.write(os.path.join(EXPPATH, 'sub_summary.p'))
                    txtfile.write('\n')

if COMPRESSOUTPUT:
    with open(os.path.join(PATHTOPARMS, 'finaljob-compression-out' + '.txt'), 'wb') as txtfile:
        txtfile.write('compression')
        txtfile.write('\n')
        #path to standard output
        txtfile.write(os.path.join(EXPPATH, 'output'))
        txtfile.write('\n')

if COMPRESSPROCESSEDOUTPUT:
    for path in ['timecourse', 'compressedtime', 'finalstructure']:
        with open(os.path.join(PATHTOPARMS, 'finaljob-compression-' + path + '.txt'), 'wb') as txtfile:
            txtfile.write('compression')
            txtfile.write('\n')
            #path to standard output
            txtfile.write(os.path.join(EXPPATH, 'proc_output', path))
            txtfile.write('\n')

CALLCOMMAND = 'python ' + os.path.join(ROOT, PROCESSINGSCRIPTNAME)
#os.chdir(NAMEOFPROCESSINGFOLDER)
PATHTOFRAMEWORK = os.path.join(ROOT, NAMEOFPROCESSINGFOLDER)
hyakcreate.mybundle_sub(PATHTOFRAMEWORK, EMAIL, NUMBEROFCORES, 1,
                        '15:00:00', 'highthroughputproc', finaljob=True)
hyakcreate.general_myscript_sub(PATHTOFRAMEWORK, CALLCOMMAND)
hyakcreate.symlinks(os.path.join(PATHTOFRAMEWORK, 'MyScript.sh'),
                    os.path.join(PATHTOFRAMEWORK,'myscript-parms'),
                    os.path.join(PATHTOFRAMEWORK,'myscript-links'), '*.txt')

#Submit the file
sub.submit_basic_hyak_framework(os.path.join(PATHTOFRAMEWORK, 'myscript-links'))
