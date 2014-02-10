import fileinput
import os
from pyrfold import hyak
from pyrfold import pyrfile


#The strip command is needed to remove the newline character
INPUTFILE = fileinput.input()
#Name of sub file
EXPERIMENTALDIRECTORY = INPUTFILE.next().strip()
#The big thing that I need is the experiment name
NAMEOFEXPERIMENT = INPUTFILE.next().strip()
#The path of the directory which contains the part
DIRECTORYPATH = INPUTFILE.next().strip()
#The name of the sub file
SUBFILE = INPUTFILE.next().strip()
#Additionally the parameters
NUMBEROFCORES = int(INPUTFILE.next().strip())
NUMBEROFNODES = int(INPUTFILE.next().strip())
NUMBEROFSIMULATIONS = int(INPUTFILE.next().strip())


os.chdir(DIRECTORYPATH)
# Gather all of the information
devicetosequence, devicetopart, devicetokinefoldparms, \
     devicetoexperimentalparms  = pyrfile.sub_file(SUBFILE)
ROOT = os.getcwd()
if not os.path.exists(os.path.join(DIRECTORYPATH, EXPERIMENTALDIRECTORY)):
    os.makedirs(os.path.join(DIRECTORYPATH, EXPERIMENTALDIRECTORY))
os.chdir(os.path.join(DIRECTORYPATH, EXPERIMENTALDIRECTORY))
pyrfile.submission_summary(devicetosequence, devicetopart,
                    devicetokinefoldparms, devicetoexperimentalparms)
hyak.framework(devicetosequence, devicetopart, devicetokinefoldparms,
    devicetoexperimentalparms, NUMBEROFCORES, NUMBEROFNODES,
    NUMBEROFSIMULATIONS, NAMEOFEXPERIMENT)

hyak.submit_multiple_nodes(os.path.join(DIRECTORYPATH, EXPERIMENTALDIRECTORY))
