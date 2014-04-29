"""
This should handle every submission step for hyak submissions

TODO - Find a better place for these submission files
"""
import os
import csv
from . import create
from .. import pyrfile

def submit_file(filepath, writedirectory, email, cores=16, backfill=True,
                        nameofexperiment='auto', nodes='auto'):
    """2014-01-24 09:52 WEV
    :param filepath: path to the hyak submission file
    :type filepath: str
    :param writedirectory: directory which the experiment files will be written
    :type writedirectory: str
    :param email: Email address which will recieve hyak updates on experiment
    :type email: str
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
    #Define all of the needed variables
    experimentname = os.path.basename(filepath).split('.')[0]
    experimentpath = os.path.join(writedirectory, experimentname)
    devicenametosubobj = pyrfile.sub_file(filepath)
    if nameofexperiment == 'auto':
        nameofexperiment = experimentname
    #Make the directory
    if not os.path.exists(experimentpath):
        os.makedirs(experimentpath)
    #Make the sub summary
    # pyrfile.submission(devicetosequence, devicetopart,
    #                 devicetokinefoldparms, devicetoexperimentalparms,
    #                 experimentpath)
    create.framework(experimentpath, devicenametosubobj, cores,
        nameofexperiment, email, backfill=backfill, nodes=nodes)
    #submit files for simulation
    submit_multiple_nodes(experimentpath)

################# SUBMISSION STUFF ####################################
def submit_multiple_nodes(parentdirectory):
    """2014-01-04 16:08 WEV
    Parent directory is the directory which contains all of all of
    the node files
    Note that the files must be named nodes
    """
    from subprocess import call
    callcommand = 'qsub ../MyBundle.sh'
    nodelist = [ls for ls in os.listdir(parentdirectory)
                    if os.path.isdir(os.path.join(parentdirectory, ls))
                    and 'node' in ls]
    for node in nodelist:
        node = os.path.join(parentdirectory, node)
        os.chdir(os.path.join(node, 'myscript-links'))
        call(callcommand, shell = True)

def submit_basic_hyak_framework(myscriptlinksdirectory):
    """ This should submit a given hyak framework
    :param myscriptlinksdirectory: This is the directory of the
        myscript-links folder that's created for all hyak frameworks
    :type myscriptlinksdirectory: str
    """
    from subprocess import call
    callcommand = 'qsub ../MyBundle.sh'
    os.chdir(myscriptlinksdirectory)
    call(callcommand, shell=True)

############## RESUBMISSION STUFF ######################################
def additional_round_submission(performdict, summarydata, foldcutoff,
            nextroundpath, poldwell, fivethreeshift, email, numsimulations):
    """ This is used to submit an additional round of simulations based on
    a greedy selection
    :param performdict:
    :type performdict:
    :param cursubname:
    :type cursubname:
    :param foldcutoff:
    :type foldcutoff:
    :param nextroundpath:
    :type nextroundpath:
    :param listofconditions:
    :type listofconditions:
    :param submissiondata: [email, numberofsimulations]
    :type submissiondata: list
    TODO - unhardcode this additional round sub
    """
    #Select the winners
    listofdevices = []
    for device in performdict:
        fail = False
        #Everypart has to fold better than the cutoff
        for part in performdict[device]:
            if performdict[device][part] < foldcutoff:
                fail = True
        if not fail:
            listofdevices.append(device)
    #The sub_summary data contains all of the parts of interest
    nextroundsubdict = {}
    for device in listofdevices:
        tempsubobj = summarydata[device]
        if fivethreeshift[0] or fivethreeshift[1]:
            relativepositionpart = tempsubobj.positionrefpart
            relpartstartstop = tempsubobj.part_start_stop(relativepositionpart)
            if fivethreeshift[0]:
                tempsubobj.winodwstart = relpartstartstop[0] - fivethreeshift[0]
            if fivethreeshift[1]:
                tempsubobj.windowstop = relpartstartstop[1] + fivethreeshift[1]
        tempsubobj.polrate = poldwell[0]
        tempsubobj.foldtimeafter = poldwell[1]
        tempsubobj.numberofsimulations = numsimulations
        nextroundsubdict[device] = tempsubobj
    #Get the name of the next round
    nextroundname = os.path.basename(nextroundpath)
    subfiledirectory = os.path.dirname(nextroundpath)
    #Print the next sub file
    pathtonewsub = os.path.join(subfiledirectory, nextroundname)
    pyrfile.filled_in_form(pathtonewsub, nextroundsubdict)
    #Submit the file that was just written
    submit_file(pathtonewsub, os.path.dirname(pathtonewsub),
        email, cores=16, backfill=False, nodes=2)

def combine_sequences_calculate_windowranges(listofsequencecomponents, partindextopartname, windowsize, shiftfraqthree):
    """should simply stitch together components and deterimine the proper
    window sequence as well as the part placement

    partindextopartname states whether something is a part or not
    This will have to be expanded to deal with things that have more than
    a single part [[index,part]]

    list of sequences contains, in order, the sequence components that will make up the entire submission
    """
    combinedsequence = ""
    halfwayindexpartlist = []
    outkinedata = KineSubData()
    for index, sequencecomponent in enumerate(listofsequencecomponents):
        if index in partindextopartname:
            #This is a part need to store part position
            #the 6 is added in a hacky fashion in order to accont for insulated
            #regions
            tempstartstop = []
            tempstartstop.append(len(combinedsequence) + 1 + 6)
            halfwayindexpartlist.append(len(sequencecomponent)/2 + len(
                                                        combinedsequence))
            combinedsequence += sequencecomponent
            tempstartstop.append(len(combinedsequence) - 6)
            outkinedata.partstartstoplist.append(tempstartstop)
            outkinedata.partnamelist.append(partindextopartname[index])
        else:
            combinedsequence += sequencecomponent
    outkinedata.windowstart = halfwayindexpartlist[0] - \
                      int(windowsize/2*shiftfraqthree)
    outkinedata.windowstop = halfwayindexpartlist[0] + \
                                int(windowsize/2*(2-shiftfraqthree))
    outkinedata.sequence = combinedsequence
    return outkinedata

