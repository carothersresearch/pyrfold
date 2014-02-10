"""
This should handle every submission step for hyak submissions

TODO - Find a better place for these submission files
"""
import os
from . import create
from .. import pyrfile

def submit_file(filepath, writedirectory, numberofsimulations=100,
                            cores=12, nameofexperiment='auto', nodes='auto'):
    """2014-01-24 09:52 WEV
    :param filepath: path to the hyak submission file
    :type filepath: str
    """
    #Define all of the needed variables
    experimentname = os.path.basename(filepath).split('.')[0]
    experimentpath = os.path.join(writedirectory, experimentname)
    (devicetosequence, devicetopart, devicetokinefoldparms,
     devicetoexperimentalparms,
     devicetoforced)  = pyrfile.sub_file(filepath)
    if nameofexperiment == 'auto':
        nameofexperiment = experimentname
    #Make the directory
    if not os.path.exists(experimentpath):
        os.makedirs(experimentpath)
    #Make the sub summary
    pyrfile.submission(devicetosequence, devicetopart,
                    devicetokinefoldparms, devicetoexperimentalparms,
                    experimentpath)
    if nodes == 'auto':
        create.framework(experimentpath, devicetosequence,
            devicetokinefoldparms, devicetoexperimentalparms, cores,
            numberofsimulations, nameofexperiment, devicetoforced)
    else:
        create.framework(experimentpath, devicetosequence,
            devicetokinefoldparms, devicetoexperimentalparms, cores,
            numberofsimulations, nameofexperiment, devicetoforced, nodes=nodes)
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


#holdovers yet to be sorted
class KineSubData:
    def __init__(self):
        self.sequence =""
        self.windowstart = 0
        self.windowstop = 0
        self.partstartstoplist = []
        self.partnamelist = []
        self.polrate = 30
        self.foldtimeafter = 1
        self.experimenttype = 2
        self.pseudoknots = 0
        self.entanglements = 0

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

