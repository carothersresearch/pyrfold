

############################ Modules ###################################
import os
import glob
import random
import math
import cPickle as pickle
#import shutil

############################ Helper Functions ##########################
def dna_to_rna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq

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
def framework(directorypath, devicenametosubobj, cores, nameofexperiment,
            email, backfill=True, nodes=1, finaljob=False):
    """ Builds a directory framework for the highthroughput processing
    of kinefold folding simulations.
    :param directorypath:
    :type directorypath:
    :param devicenametosubobj:
    :type devicenametosubobj:
    :param numberofsimulations: The total number of simulations to run for a
        specific simulations.
    :type numberofsimulations: int
    :param cores: number of cores a node will have (this is important for
        submission like things)
    :type cores: int
    :param nameofexperiment: Name of the experiment that will be displayed in
        hyak display and referenced in emails.
    :type nameofexperiment: str
    :param nodes: The number of nodes to run the simulations on. If 'auto' it
        will make an approximation for the nodes based on simulations requested
    :type nodes: int, str
    """
    #First write the output framework
    create_output_directory(directorypath, devicenametosubobj)
    #Write the submission summary file:
    write_sub_summary(directorypath, devicenametosubobj)
    #Need a sophisicated function to determine where jobs go
    nodestodevices = calculate_number_of_nodes(devicenametosubobj, nodes)
    #write all of the dat files (sequence files)
    for node in nodestodevices:
        tempnodename = 'node-' + str(node)
        nodedirectory = os.path.join(directorypath, tempnodename)
        #make the new directory
        new_directory(nodedirectory)
        datdirectory = os.path.join(nodedirectory, 'dat')
        new_directory(datdirectory)
        parmdirectory = os.path.join(nodedirectory, 'myscript-parms')
        new_directory(parmdirectory)
        linkdirectory = os.path.join(nodedirectory, 'myscript-links')
        new_directory(linkdirectory)
        mybundle_sub(nodedirectory, email, cores, 1, '15:00:00',
            nameofexperiment + '-' + tempnodename, backfill)
        #write_MyScript.sh
        myscript_sub_kinefold(nodedirectory)
        #Write all of the dat files
        write_dat_files(datdirectory, nodestodevices[node],
                    devicenametosubobj)
        #Write the req files that kinefold requires
        outputdirectory = os.path.join(directorypath, 'output')
        write_req_files(parmdirectory, outputdirectory, datdirectory,
                                    nodestodevices[node], devicenametosubobj)
        #Make symbolic links for devices
        myscriptpath = os.path.join(nodedirectory, 'MyScript.sh')
        symlinks(myscriptpath, parmdirectory, linkdirectory, '*.req')

def write_sub_summary(directory, devicenametosubobj):
    """This will create a pickle structure for the devicenametosubobj
    This will serve both as a condensed record for what was submitted as well
    as something that can be used for part references"""
    pathtopicklefile = os.path.join(directory, 'sub_summary.p')
    with open(pathtopicklefile, 'wb') as picklefile:
        pickle.dump(devicenametosubobj, picklefile, protocol=2)

def write_req_files(parmdirectory, outputpath, datdirectory, listofdevices,
                                                        devicenametosubobj):
    """(str, int, tuple) -> write .req files for device
    This will write all of the neccessary req files for kinefold to process
    """
    random.seed() #uses system time to initialize random generator
    #EXPERIMENTAL VARIABLES
    for devicename in listofdevices:
        numberofsimulations = devicenametosubobj[devicename].numberofsimulations
        polrate, requestedtime = devicenametosubobj[devicename].kine_folding_data()
        psudoknots = devicenametosubobj[devicename].pseudoknots
        entanglements = devicenametosubobj[devicename].entanglements
        forcedhelixes = devicenametosubobj[devicename].forcedhelixes
        exptype = devicenametosubobj[devicename].experimenttype
        for i in range(numberofsimulations):
            #create name for req files and outputs
            reqname = str(i + 1).zfill(3) + devicename
            #open unique req file
            with open(os.path.join(parmdirectory, 'job.' + reqname
                                                        + '.req'), 'wb') as f:
                #Write random number seed
                f.write(str(int(round(random.random() * 10000))).zfill(4))
                #Write all output directories
                fileext = ['.p', '.e' , '.rnm', '.rnms', '.rnml', '.rnm2']
                for ext in fileext:
                    f.write('\n')
                    f.write(os.path.join(outputpath, devicename, reqname + ext))
                #Write .dat directory
                f.write('\n' + os.path.join(datdirectory, devicename + '.dat'))
                #Write 0 RNA 1 for DNA
                f.write('\n' + str(0))
                # helix minimum free energy in kcal/mol: 6.3460741=10kT
                f.write('\n' + str(6.3460741))
                #just something needed
                f.write('\n' + str(10000000))
                #requested folding time
                f.write('\n' + str(requestedtime))
                #psudoknots 1 = yes, 0 = no
                f.write('\n' + str(int(psudoknots)))
                #entanglements 1 = yes, 0 = no
                f.write('\n' + str(int(entanglements)))
                # simulation type: 1=renaturation; 2 20 =cotrans. @ 20msec/nt
                if exptype == 2:
                    f.write('\n' + '2 ' + str(polrate))
                else:
                    f.write('\n' + '1')
                for forc in forcedhelixes:
                    f.write('\n' + 'F ' + str(forc[0]) + ' ' + str(forc[1]) +
                                                         ' ' + str(forc[2]))
                #filename and filename.zip
                f.write('\n' + devicename + '\n' + devicename + '.zip' + '\n')

def symlinks(exefilepath, filepath, linkdirectory, ext):
    """(str, str(path), str(path)) -> symlinks
    This script connects all files filepath to exefile and places the
    symlink in linkdirectory

    This is specific to *.req files
    """
    # collect all of the file names
    filelist = glob.glob(os.path.join(filepath, ext))
    for fi in filelist:
        fi = os.path.basename(fi)
        os.symlink(exefilepath, os.path.join(linkdirectory, fi))

def create_output_directory(directorypath, devicenametosubobj):
    """Writes the framework for the output directory"""
    for devicename in devicenametosubobj:
        tempath = os.path.join(directorypath, 'output', devicename)
        new_directory(tempath)

def write_dat_files(datdirectory, listofdevices, devicenametosubobj):
    """This will fill in the tempnode/dat directory with the neccessary
    sequence information for processing"""
    for device in listofdevices:
        with open(os.path.join(datdirectory,(device + '.dat')), 'wb') as f:
            f.write("<" + device + '\n')
            #truncating device from tuple (seq,winStart,winStop)
            seq = devicenametosubobj[device].sequence
            start = devicenametosubobj[device].windowstart
            stop = devicenametosubobj[device].windowstop
            #accounting for shift in frame str[n:c] doesn't actually read
            #through to c it will stop at c-1
            seq = seq[start - 1: stop]
            seq = dna_to_rna(seq)
            f.write(seq + '\n')

def calculate_number_of_nodes(devicenametosubobj, nodes):
    """Since the number of simulations done for every device is dicated now
    in the spreadsheet we have to come up with a clever way of spreading out
    the jobs"""
    outdict = {}
    #calculate the total number of jobs
    numberofjobs = 0
    for device in devicenametosubobj:
        numberofjobs += devicenametosubobj[device].numberofsimulations
    #Calculate the number of jobs to go per node
    if nodes == 'auto':
        nodes = 2
    jobspernode = math.ceil(numberofjobs/nodes)
    for node in range(nodes):
        outdict[node] = []
    #Distribute the devices tot he dictionaries
    countforasinglenode = 0
    node = 0
    for device in devicenametosubobj:
        outdict[node].append(device)
        countforasinglenode += devicenametosubobj[device].numberofsimulations
        if countforasinglenode >= jobspernode:
            countforasinglenode = 0
            node += 1
            #this is a sloppy way to make sure we don't intoduce to many nodes
            if node > nodes:
                node = nodes
    return outdict

def mybundle_sub(directorypath, email, cores, nodes, walltime,
            nameofexperiment, backfill=False, finaljob=False):
    """(str,num,num,str) -> .sh script for submission
    nameofexperiment <- reference for emails
    This script is needed for hyak checkpointing
    See the comments for the complete description
    """
    #Writing file
    pathtomybundle = os.path.join(directorypath, 'MyBundle.sh')
    with open(pathtomybundle,'wb') as f:
        f.write("#!/bin/bash")
        ## RENAME FOR YOUR JOB
        f.write("\n#PBS -N " + nameofexperiment)
        ## EDIT FOR YOUR JOB
        ## For 16 core nodes
        if backfill:
            f.write("\n#PBS -l nodes=" + str(nodes) + ":ppn=" + str(cores) +
                ",mem=22gb,feature=" + str(cores) + "core -q bf")
        else:
            f.write("\n#PBS -l nodes=" + str(nodes) + ":ppn=" + str(cores) +
                ",mem=22gb,feature=" + str(cores) + "core")
        ## WALLTIME DEFAULTS TO ONE HOUR - ALWAYS SPECIFY FOR LONGER JOBS
        ## If the job doesn't finish in 10 minutes, cancel it
        f.write("\n#PBS -l walltime=:" + walltime)
        ## EDIT FOR YOUR JOB
        ## Put the STDOUT and STDERR from jobs into the below directory
        f.write("\n#PBS -o " + directorypath)
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
        f.write("\n#PBS -d " + os.path.join(directorypath, "myscript-links"))
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
        #f.write("\ntouch ../completed-jobs/null")
        f.write('\nmodule load epd_7.3_2')
        # remove sym links to scripts that have completed
        # this is the "checkpointing" operation.
        #f.write("\nfor i in ../completed-jobs/*; do F=`echo $i | cut -d / -f 3`; rm $F; done")
        # run the rest of them
        f.write("\nfind . -name job\* | parallel -j $HYAK_SLOTS --joblog paralleljobs.log --resume")
        if finaljob:
            f.write("\nfind . -name finaljob\* | parallel -j $HYAK_SLOTS --joblog paralleljobs.log --resume")
        f.write("\nexit 0\n")
    os.chmod(pathtomybundle, 0777)

def framework_shell(pathtoframework, finaljob=False):
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
    # new_directory(os.path.join(pathtoframework,
    #     'completed-jobs'))
    return os.path.join(pathtoframework,
        'myscript-parms')

def general_myscript_sub(directory, callfilecommand):
    """ -> .sh file for submission
    The call file command should be a string of txt that oulines
    what exactly is being called and how it is being modife
    This script writes a .sh file that will actually do the work of
    submitting a job to hyak
    """
    filepath = os.path.join(directory, 'MyScript.sh')
    #open file
    f = open(filepath, 'wb')
    f.write("#!/bin/bash")
    ##################################################################
    # Set some parameters here
    ##################################################################
    f.write("\nTopDir=" + directory)
    # All of the parameter files should be in one directory
    f.write("\nparmdir=$TopDir/myscript-parms")
    # We write output into a standard location
    # output file has the same name as the script
    f.write("\noutputdir=$TopDir/output")
    # To support our primitive checkpointing scheme, provide a location
    # for storing semaphores to indicate which jobs have completed
    # f.write("\ncompleted=$TopDir/completed-jobs")
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
    #f.write("\ntouch $completed/$0\n")
    f.close()
    os.chmod(filepath, 0777)

def myscript_sub_kinefold(parentdirectory):
    """ -> .sh file for submission

    This script writes a .sh file that will actually do the work of
    submitting a job to hyak
    """
    myscriptpath = os.path.join(parentdirectory, 'MyScript.sh')
    #open file
    with open(myscriptpath, 'wb') as f:
        f.write("#!/bin/bash")
        ##################################################################
        # Set some parameters here
        ##################################################################
        f.write("\nTopDir=" + parentdirectory)
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
        f.write("\n/gscratch/rna/compiled_binaries/kinefold/kinefold_long_static $ParmFile -noprint >/dev/null")
    os.chmod(myscriptpath, 0777)


