"""
This is a python script that is designed to operate with the hyak
framework for parallel processing
"""
import shutil
import fileinput
import os
from pyrfold.hyak import process as hyakp

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
elif PROCESSINGTYPE == 'timecourse':
    #Pull all of the boolean decision variables
    MAKEREFERENCESTRUCTURECOMP = int(INPUTFILE.next().strip())
    SINGLEDIRECTORY = int(INPUTFILE.next().strip())
    ADDITIONALROUNDOFPROCESS = int(INPUTFILE.next().strip())
    if MAKEREFERENCESTRUCTURECOMP:
        #NEED TO BUILD THIS
        if ADDITIONALROUNDOFPROCESS:
            pass
        else:
            pass
        pass
    else:
        FOLDERTOPROCESS = INPUTFILE.next().strip()
        PROCESSEDDIRECTORY = INPUTFILE.next().strip()
        NAMEOFFOLDERTOPROCESS = os.path.basename(FOLDERTOPROCESS)
        if SINGLEDIRECTORY:
            #remove unwanted files
            hyakp.clean_up_output_data(FOLDERTOPROCESS, singlefile=True)
            #create the timecourse files
            hyakp.timecourse(FOLDERTOPROCESS, PROCESSEDDIRECTORY,
                singledirectory=True)
            if ADDITIONALROUNDOFPROCESS:
                #Need to build this
                pass
        else:
            #remove unwanted files
            hyakp.clean_up_output_data(FOLDERTOPROCESS)
            #create the timecourse files
            hyakp.timecourse(FOLDERTOPROCESS, PROCESSEDDIRECTORY)
            if ADDITIONALROUNDOFPROCESS:
                #Need to build this
                pass
