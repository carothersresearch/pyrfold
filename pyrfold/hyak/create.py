

############################ Modules ###################################
import os
import glob
import random
import math
#import shutil

############################ Helper Functions ##########################
def dna_to_rna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq

def kine(devicetosequence, devicetokinefoldparms):
    """({str,tuple}, {str, tuple}) -> {str,tuple}

    The required input of kinefold is different than what was put in the
    input sheet, consequently it must be converted
    """
    #First reprocess devicetokinefoldparms to have the requested information
    #for kinefold (rate [ms/nt], wanted time [ms])
    for device in devicetokinefoldparms:
        #stop position of device - start position of device + 1
        windowsize = (devicetosequence[device][2] -
                        devicetosequence[device][1] + 1)
        #original is tuple (pol (nt/s), extratime)
        original = devicetokinefoldparms[device]
        if float(original[0]) == 0.0: #no folding
            rate = 0
        else:
            rate = float(1)/float(original[0])*float(1000) #ms/nt
        #orignal[1] is additional time of folding
        requestedtime = rate * windowsize + original[1] * 1000 #ms
        devicetokinefoldparms[device] = (rate, requestedtime)
    return devicetokinefoldparms

def number_of_nodes(numberofsimulations, numberofexperiments):
    """(int, int) --> int
    This will calculate the number of nodes to use for a given experiment
    it takes ~4 hours for a 12 core node to do 5000 exp
    """
    numberofthingspernode = 5000.0 #This is an educated guess
    totalnumberofthings = numberofsimulations * numberofexperiments
    return int(math.ceil(totalnumberofthings / numberofthingspernode))

def new_directory(directory):
    """(str)->directory
    simple code to make a directory in the root if it does not exist
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

############################ Core Functions   ##########################
def framework(directorypath, devicetosequence, devicetokinefoldparms,
    devicetoexperimentalparms, cores, numberofsimulations,
    nameofexperiment, forcedhelixes, email, nodes='auto'):
    """(dict, dict, dict) -> .dat file containing all of the information
    edited 2013-10-02 WEV
    """
    # import local dependencies
    #import dependencies
    root = directorypath
    numberofexperiments = len(devicetosequence.keys())
    if nodes == 'auto':
        nodes = number_of_nodes(numberofsimulations, numberofexperiments)
    #Now we break the work up to number of nodes
    numexppernode = int(math.ceil(numberofexperiments/float(nodes)))
    #write all of the dat files (sequence files)
    #Convert the devicetokinefoldparms to the required units
    #see dockstring of convert_kine_input for more info
    devicetokinefoldparms = kine(devicetosequence,
                                         devicetokinefoldparms)
    nodecount = 0
    for counter, device in enumerate(devicetosequence):
        #Check to see the number of directories
        if counter % numexppernode == 0:
            os.chdir(root)
            nodecount += 1
            #Make the new directory
            tempnodename = 'node-' + str(nodecount)
            new_directory(tempnodename)
            os.chdir(tempnodename)
             #writing directy system
            new_directory('dat')
            #new_directory('output')
            new_directory('myscript-parms')
            new_directory('myscript-links')
            new_directory('completed-jobs')
            #wrtie MyBundle master submission script
            mybundle_sub(email, int(cores), 1, '15:00:00',
                            nameofexperiment + '-' + tempnodename)
            #write_MyScript.sh
            myscript_sub()
        #create $device.dat file for every device being tested
        f = open(os.path.join(root, tempnodename,'dat',(device + '.dat'))
                                                                    , 'wb')
        f.write("<" + device + '\n')
        #truncating device from tuple (seq,winStart,winStop)
        seq = devicetosequence[device][0]
        start = int(devicetosequence[device][1])
        stop = int(devicetosequence[device][2])
        #accounting for shift in frame str[n:c] doesn't actually read
        #through to c it will stop at c-1
        seq = seq[start - 1: stop]
        seq = dna_to_rna(seq)
        f.write(seq + '\n')
        f.close()
        #Submit req files
        req_files(device, int(numberofsimulations),
                        devicetokinefoldparms[device],
                        devicetoexperimentalparms[device],
                        os.path.join(root, tempnodename, 'myscript-parms'),
                        os.path.join(root, 'output'), forcedhelixes[device])


    os.chdir(root)
    #create symlinks between myscript-parms and MyScript.sh to be placed in
    #myscirpt-links
    nodelist = [ls for ls in os.listdir(root) \
                        if os.path.isdir(os.path.join(root,ls))]
    for node in nodelist:
        node = os.path.join(root, node)
        os.chdir(node)
        symlinks('MyScript.sh', os.path.join(node,'myscript-parms'),
                     os.path.join(node,'myscript-links'), '*.req')
        os.chdir(root)
    #make the output directory
    new_directory('output')
    #create a output directory for every device
    for device in devicetosequence:
        #create directory for every device
        new_directory(os.path.join(root, 'output', device))

def symlinks(exefile, filepath, linkdest, ext):
    """(str, str(path), str(path)) -> symlinks
    This script connects all files filepath to exefile and places the
    symlink in linkdest

    This is specific to *.req files
    """
    # collect all of the file names
    root = os.getcwd()
    os.chdir(filepath)
    filelist = glob.glob(ext)
    for fi in filelist:
        os.symlink(('../' + exefile), os.path.join(linkdest, fi))
    os.chdir(root)

def mybundle_sub(email, cores, nodes, walltime, nameofexperiment):
    """(str,num,num,str) -> .sh script for submission
    nameofexperiment <- reference for emails
    This script is needed for hyak checkpointing
    See the comments for the complete description
    """
    root = os.getcwd()
    #Writing file
    f = open('MyBundle.sh','wb')
    f.write("#!/bin/bash")
    ## RENAME FOR YOUR JOB
    f.write("\n#PBS -N " + nameofexperiment)
    ## EDIT FOR YOUR JOB
    ## For 16 core nodes
    f.write("\n#PBS -l nodes=" + str(nodes) + ":ppn=" + str(cores) + ",mem=22gb,feature=" + str(cores) + "core -q bf")
    ## WALLTIME DEFAULTS TO ONE HOUR - ALWAYS SPECIFY FOR LONGER JOBS
    ## If the job doesn't finish in 10 minutes, cancel it
    f.write("\n#PBS -l walltime=:" + walltime)
    ## EDIT FOR YOUR JOB
    ## Put the STDOUT and STDERR from jobs into the below directory
    f.write("\n#PBS -o " + root)
    ## Put both the stderr and stdout into a single file
    f.write("\n#PBS -j oe")
    f.write("\n#PBS -M " + email)
    ## a mail is sent when the job is aborted by the batch system.
    ## b mail is sent when the job begins execution.
    ## e mail is sent when the job terminates.
    f.write("\n#PBS -m abe")
    ## EDIT FOR YOUR JOB
    ## Sepcify the working directory for this job bundle
    ## NOTE: we run from the LINKS directory - NOT the dir
    ## containing the actual script files.  This is CRITICAL
    ## for our checkpoinitng scheme
    f.write("\n#PBS -d " + os.path.join(root, "myscript-links"))
    ## If you can't run as many tasks as there are cores due to
    ## memory constraints you can simply set HYAK_SLOTS to a
    ## number instead.
    ## HYAK_SLOTS=4
    f.write("\nHYAK_SLOTS=`wc -l < $PBS_NODEFILE`")
    ## Prevent tasks from exceeding the total RAM of the node
    ## Requires HYAK_SLOTS to be set to number of tasks started.
    f.write("\nNODEMEM=`grep MemTotal /proc/meminfo | awk '{print $2}'`")
    f.write("\nNODEFREE=$((NODEMEM-2097152))")
    f.write("\nMEMPERTASK=$((NODEFREE/HYAK_SLOTS))")
    f.write("\nulimit -v $MEMPERTASK")
    ## need a file in the completed-jobs dir in order for the 'for'
    ## expression below to evaluate correctly - stupid *nix quirk
    f.write("\ntouch ../completed-jobs/null")
    f.write('\nmodule load epd_7.3_2')
    # remove sym links to scripts that have completed
    # this is the "checkpointing" operation.
    f.write("\nfor i in ../completed-jobs/*; do F=`echo $i | cut -d / -f 3`; rm $F; done")
    # run the rest of them
    f.write("\nfind . -name job\* | parallel -j $HYAK_SLOTS")
    f.write("\nexit 0\n")
    f.close()
    os.chmod('MyBundle.sh', 0777)

def framework_shell(pathtoframework):
    """2013-12-16 16:02 WEV
    This will create all of the directories needed for the hyak
    framework
    This is being built to support high throughput processing functionality
    """
    new_directory(pathtoframework)
    new_directory(os.path.join(pathtoframework, 'output'))
    new_directory(os.path.join(pathtoframework,
        'myscript-links'))
    new_directory(os.path.join(pathtoframework,
        'myscript-parms'))
    new_directory(os.path.join(pathtoframework,
        'completed-jobs'))
    return os.path.join(pathtoframework,
        'myscript-parms')

def general_myscript_sub(callfilecommand):
    """ -> .sh file for submission
    The call file command should be a string of txt that oulines
    what exactly is being called and how it is being modife
    This script writes a .sh file that will actually do the work of
    submitting a job to hyak
    """
    root = os.getcwd()
    #open file
    f = open('MyScript.sh', 'wb')
    f.write("#!/bin/bash")
    ##################################################################
    # Set some parameters here
    ##################################################################
    f.write("\nTopDir=" + root)
    # All of the parameter files should be in one directory
    f.write("\nparmdir=$TopDir/myscript-parms")
    # We write output into a standard location
    # output file has the same name as the script
    f.write("\noutputdir=$TopDir/output")
    # To support our primitive checkpointing scheme, provide a location
    # for storing semaphores to indicate which jobs have completed
    f.write("\ncompleted=$TopDir/completed-jobs")
    ##################################################################
    # $0 expands to the name of the shell of script being run
    f.write("\nParmFile=`echo $0 | cut -d / -f 2`")
    f.write("\nParmFile=$parmdir/$ParmFile")
    #Checkpointing
    f.write("\nrm $outputdir/$0 2>/dev/null")
    # Do your work
    f.write("\n" + callfilecommand + ' $ParmFile')
    # We're done!
    # CHECKPOINTING: Create a semaphore to prevent this script instance
    # from rerunning if the job bundle is resubmitted
    f.write("\ntouch $completed/$0\n")
    f.close()
    os.chmod('MyScript.sh', 0777)

def myscript_sub():
    """ -> .sh file for submission

    This script writes a .sh file that will actually do the work of
    submitting a job to hyak
    """
    root = os.getcwd()
    #open file
    f = open('MyScript.sh', 'wb')
    f.write("#!/bin/bash")
    ##################################################################
    # Set some parameters here
    ##################################################################
    f.write("\nTopDir=" + root)
    # All of the parameter files should be in one directory
    f.write("\nparmdir=$TopDir/myscript-parms")
    # We write output into a standard location
    # output file has the same name as the script
    f.write("\noutputdir=$TopDir/output")
    # To support our primitive checkpointing scheme, provide a location
    # for storing semaphores to indicate which jobs have completed
    f.write("\ncompleted=$TopDir/completed-jobs")
    ##################################################################

    # $0 expands to the name of the shell of script being run
    f.write("\nParmFile=`echo $0 | cut -d / -f 2`")
    f.write("\nParmFile=$parmdir/$ParmFile")
    #Checkpointing
    #f.write("\nrm $outputdir/$0 2>/dev/null")
    # Do your work
    f.write("\n/gscratch/esci/carothers/kinefold/kinefold_long_static $ParmFile -noprint >/dev/null")
    # We're done!
    # CHECKPOINTING: Create a semaphore to prevent this script instance
    # from rerunning if the job bundle is resubmitted
    f.write("\ntouch $completed/$0\n")
    f.close()
    os.chmod('MyScript.sh', 0777)

def req_files(devicename, numberofsimulations, parms, expparms, filepath,
                outputpath, forcedhelixes):
    """(str, int, tuple) -> write .req files for device
    tuple contains (rate ms/nt, total time)
    expparms tuple contains (cotrans/renaturation, psudoknots, entanglements)
    """
    root = os.getcwd()
    random.seed() #uses system time to initialize random generator
    #EXPERIMENTAL VARIABLES
    #psudoknots 1 = yes, 0 = no
    psudoknots = expparms[1]
    #entanglements 1 = yes, 0 = no
    entanglements = expparms[2]
    for i in range(numberofsimulations):
        #create name for req files and outputs
        reqname = str(i + 1).zfill(3) + devicename
        #open unique req file
        f = open(os.path.join(filepath, 'job.' + reqname
                                        + '.req'), 'wb')
        #Write random number seed
        f.write(str(int(round(random.random() * 10000))).zfill(4))
        #Write all output directories
        fileext = ['.p', '.e' , '.rnm', '.rnms', '.rnml', '.rnm2']
        for ext in fileext:
            f.write('\n')
            f.write(os.path.join(outputpath, devicename, reqname + ext))
        #Write .dat directory
        f.write('\n' + os.path.join(root, 'dat', devicename + '.dat'))
        #Write 0 RNA 1 for DNA
        f.write('\n' + str(0))
        # helix minimum free energy in kcal/mol: 6.3460741=10kT
        f.write('\n' + str(6.3460741))
        #just something needed
        f.write('\n' + str(10000000))
        #requested folding time
        f.write('\n' + str(parms[1]))
        #psudoknots 1 = yes, 0 = no
        f.write('\n' + str(int(psudoknots)))
        #entanglements 1 = yes, 0 = no
        f.write('\n' + str(int(entanglements)))
        # simulation type: 1=renaturation; 2 20 =cotrans. @ 20msec/nt
        if int(expparms[0]) == 2:
            f.write('\n' + '2 ' + str(parms[0]))
        else:
            f.write('\n' + '1')
        for forc in forcedhelixes:
            f.write('\n' + 'F ' + str(forc[0]) + ' ' + str(forc[1]) +
                                                 ' ' + str(forc[2]))
        #filename and filename.zip
        f.write('\n' + devicename + '\n' + devicename + '.zip' + '\n')
        f.close()
