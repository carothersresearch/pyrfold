"""
This should handle every submission step for hyak submissions

TODO - Find a better place for these submission files
"""
import os
import csv
from . import create
from .. import pyrfile

def submit_file(filepath, writedirectory, email, numberofsimulations=100,
                            cores=12, nameofexperiment='auto', nodes='auto'):
    """2014-01-24 09:52 WEV
    :param filepath: path to the hyak submission file
    :type filepath: str
    :param writedirectory: directory which the experiment files will be written
    :type writedirectory: str
    :param email: Email address which will recieve hyak updates on experiment
    :type email: str
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
    create.framework(experimentpath, devicetosequence,
        devicetokinefoldparms, devicetoexperimentalparms, cores,
        numberofsimulations, nameofexperiment, devicetoforced, email,
         nodes=nodes)
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
def additional_round_submission(performdict, cursubpath, nextroundpath,
        listsofconditions, submissiondata, foldcutoff=None):
    """ This is used to submit an additional round of simulations based on
    a greedy selection
    :param performdict: A dictionary that contains the device names:foldingfreq
    :type performdict: dict
    :param cursubpath: Path to the submission.csv of the previous round
    :type cursubpath: str
    :param foldcutoff: the minimum folding frequency allowed
    :type foldcutoff: float
    :param nextroundpath: directory path to the next round
    :type nextroundpath: str
    :param listsofconditions: list of lists of all of the exp conditions
    :type listsofconditions: list [[polrate, dwelltime, fivepos, threepos],...]
    :param submissiondata: [email, numberofsimulations]
    :type submissiondata: list

    TODO - unhardcode this additional round sub
    """
    #Select the winners
    listofdevices = []
    if foldcutoff:
        for device in performdict:
            fail = False
            #Everypart has to fold better than the cutoff
            for part in performdict[device]:
                if performdict[device][part] < foldcutoff:
                    fail = True
            if not fail:
                listofdevices.append(device)
    else:
        for device in performdict:
            listofdevices.append(device)
    #Pass the devices with the former spreadsheet to make a new spreadsheet
    nextsubfilelist = []
    with open(cursubpath, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        nextsubfilelist.append(next(reader))
        for row in reader:
            if row[1] in listofdevices:
                for count, listofconditions in enumerate(listsofconditions):
                    temprow = row[:]
                    #Now have to change the values 18, 19, 20
                    if count == 0:
                        pass
                    else:
                        temprow[1] += '--' + str(count)
                    #windowstart = partstart - shift
                    temprow[2] = int(temprow[19]) - listofconditions[2]
                    #Right window
                    temprow[3] = listofconditions[3] + int(temprow[20])
                    #Pol Rate
                    temprow[4] = listofconditions[0]
                    #dwell time
                    temprow[5] = listofconditions[1]
                    nextsubfilelist.append(temprow)
    #Get the name of the next round
    nextroundname = os.path.basename(nextroundpath)
    subfiledirectory = os.path.dirname(cursubpath)
    #Print the next sub file
    pathtonewsub = os.path.join(subfiledirectory, nextroundname + '_sub.csv')
    with open(pathtonewsub , 'w') as csvfile:
        writer = csv.writer(csvfile)
        for row in nextsubfilelist:
            writer.writerow(row)
    #Submit the file that was just written
    submit_file(pathtonewsub, os.path.dirname(pathtonewsub),
        submissiondata[0], numberofsimulations=submissiondata[1])

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

