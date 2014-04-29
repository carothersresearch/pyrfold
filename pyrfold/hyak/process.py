"""This module exists to process all of the kinefold output data
into a form which is actually useful

Things TODO
- ADD NO PARTLIST CONDITON TO FINALSTRUCTURES
"""

########################################################################
############################ Modules ###################################
########################################################################
import os
import glob
import cPickle as pickle
from collections import Counter
from . _convert import KineRmnOutput
import errno

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

########################################################################
############################ Core Process     ##########################
########################################################################
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def timecourse(rawoutputfolder, processeddirectory,
                                        singledirectory=False):
    """2014-01-26 15:31 WEV
    :param rawoutputfolder: name of the folder which contians kinefold raw
        output
    :type rawoutputfolder: str
    :param processeddirectory: name of the directory to print the data
    :param singledirectory: Determines whether it should dig in a directory or
                            if it should use the rawoutputfolder excusively
    output of this will be a pickled data structure which will be a
    dict structure with the following classes:
    'time' : list
    'dotbracket' : list
    'kine' : list
    'total bases' : int
    'base count' : list
    'helix list' : list
    'energy' [kcal/mol] : list
    'type' : str (this is the type of experimnet renaturation)
    """
    #Make the directory to write the files to
    outputfilename = os.path.join(processeddirectory, 'timecourse')
    if not os.path.exists(outputfilename):
        make_sure_path_exists(outputfilename)
    #collect all of the experiment files
    if not singledirectory:
        experimentdirectorylist = os.walk(rawoutputfolder).next()[1]
    else:
        experimentdirectorylist = [rawoutputfolder]
    for directoryname in experimentdirectorylist:
        if not singledirectory:
            temprnmfilelist = glob.glob(os.path.join(rawoutputfolder,
                                                directoryname, '*.rnm'))
        else:
            temprnmfilelist = glob.glob(os.path.join(rawoutputfolder,
                                                            '*.rnm'))
            directoryname = os.path.basename(directoryname)
        dicttopickle = {}
        for rnmfile in temprnmfilelist:
            rnmname = os.path.basename(rnmfile).split('.')[0]
            rnmnumber = int(rnmname[0:3]) #should pull just the runnumber
            #open the file for reading
            tempdictionary = {}
            tempdictionary['time'] = []
            tempdictionary['dotbracket'] = []
            # tempdictionary['kine'] = []
            sequencelist = []
            # tempdictionary['base count'] = []
            # tempdictionary['kine helix'] = []
            tempdictionary['energy'] = []
            with open(rnmfile, 'r') as rnmopenfile:
                for i, line in enumerate(rnmopenfile):
                    if i in [0, 1]:
                        continue
                    elif i % 2 == 0:
                        kineobj = KineRmnOutput(line)
                    elif (i+1) % 2 == 0: #should account for odds
                        kineobj.addHelix(line)
                        kineobj.generateDotBracket()
                        tempdictionary['time'].append(kineobj.time)
                        tempdictionary['dotbracket'].append(
                            kineobj.dotbracket)
                        # tempdictionary['kine'].append(kineobj.structure)
                        # tempdictionary['base count'].append(
                        #     kineobj.basenumber)
                        sequencelist.append(kineobj.sequence)
                        # tempdictionary['kine helix'].append(kineobj.helix)
                        tempdictionary['energy'].append(kineobj.freeenergy)
                        tempdictionary['type'] = kineobj.type
                        tempdictionary['total bases'] = kineobj.basetotal
            tempdictionary['sequence'] = sequencelist[-1]
            dicttopickle[rnmnumber] = tempdictionary
        #Now to write the output for a single experiment
        temppickledest = os.path.join(outputfilename, directoryname + '.p')
        pickle.dump(dicttopickle, open(temppickledest, 'wb'), protocol=2)

def finalstructure(processeddirectory, subsummarydict=None,
                 singledirectoryname=None):
    """This will go through the previously picked timecoures data
    and pull out part structures for all of the parts that are outlined
    in the summart dictionary
    :param processeddirectory: path to the processed output data
    :type processeddirectory: str
    :param summardict: dictionary with key being part name and value being
        [referencepart name, part start, part stop] these are adjusted to the
        submitted files
    :param singledirectoryname: this is for highthroughput processing tells
        the function to either dig through all of the completed pickles or to
        process a single pickle if not none this will be the name of the
    :return: builds a a system of pickle files
    """
    #Sort the paths
    pathtotimecourse = os.path.join(processeddirectory, 'timecourse')
    pathtofinalstructures = os.path.join(processeddirectory, 'finalstructure')
    if not os.path.exists(pathtofinalstructures):
        make_sure_path_exists(pathtofinalstructures)
    #collect all of the experiment files
    if not singledirectoryname:
        pickletimelist = [ls for ls in os.listdir(pathtotimecourse)
                    if os.path.isfile(os.path.join(pathtotimecourse, ls))]
    else:
        pickletimelist = [singledirectoryname]
    #Now to go through all of these files
    for pick in pickletimelist:
        if subsummarydict:
            partdict = subsummarydict[pick.split('.')[0]].adjustedpartstartstop
        else:
            partdict = None
        pathtopickle = os.path.join(pathtotimecourse, pick)
        with open(pathtopickle, 'rb') as picktemp:
            timecoursedict = pickle.load(picktemp)
        outpickdict = create_finalpart_dict(timecoursedict, partdict)
        pathtopickle = os.path.join(pathtofinalstructures, pick)
        with open(pathtopickle, 'wb') as topickle:
            pickle.dump(outpickdict, topickle, protocol=2)

def create_finalpart_dict(timecoursedict, dictofparts=None):
    """helper function to final structure
    :param timecoursedict: complete timecoursedict for a experiment
    :type timecoursedict: dict
    :param listsofparts: listoflists of parts [[partname, start, stop],...]
    :type listofparts: list
    """
    tempdict = {}
    listofdotbracket = []
    listofsequences = []
    dictofenergy = {}
    for runnumber in timecoursedict:
        if timecoursedict[runnumber]['dotbracket']:
            listofdotbracket.append(
                timecoursedict[runnumber]['dotbracket'][-1])
            listofsequences.append(
                timecoursedict[runnumber]['sequence'])
            dictofenergy[timecoursedict[runnumber]['dotbracket'][-1]] = timecoursedict[runnumber]['energy'][-1]
    dotbracketcount = Counter(listofdotbracket)
    sequencecount = Counter(listofsequences)
    totalcount = float(sum(dotbracketcount.viewvalues()))
    for dotbrac in dotbracketcount:
        dotbracketcount[dotbrac] /= totalcount
    tempdict['energy'] = dictofenergy
    tempdict['sequence'] = sequencecount.keys()[0]
    tempdict['dotbracket'] = dotbracketcount
    if dictofparts:
        tempdict['parts'] = {}
        for part in dictofparts:
            startstop = dictofparts[part]
            tempdict['parts'][part] = {}
            start = startstop[0] - 1 #shift to account for indexing
            stop = startstop[1]
            listofdotbracket = []
            listofsequences = []
            for runnumber in timecoursedict:
                listofdotbracket.append(
                    timecoursedict[runnumber]['dotbracket'][-1][start:stop])
                listofsequences.append(
                    timecoursedict[runnumber]['sequence'][start:stop])
            dotbracketcount = Counter(listofdotbracket)
            sequencecount = Counter(listofsequences)
            if len(sequencecount) > 1:
                print "Sequences DONOT match up for {}".format(part)
            #Convert count to average
            totalcount = float(sum(dotbracketcount.viewvalues()))
            for dotbrac in dotbracketcount:
                dotbracketcount[dotbrac] /= totalcount
            tempdict['parts'][part]['sequence'] = sequencecount.keys()[0]
            tempdict['parts'][part]['dotbracket'] = dotbracketcount
    return tempdict

########################################################################
################# Gen.  Destruction  ###################################
########################################################################
def clean_up_output_data(experimentfolder, singlefile=False):
    """2014-01-08 17:12 WEV
    Removes all of the undesired files from the output of kinefold
    """
    otherfilestodelete = ['*.rnms', '*.rnm2','*.e'
                                    ,'*.p', '*.i', '*.rnm']
    if not singlefile:
        experimentfolder = os.path.join(experimentfolder, 'output')
        #change to this path
        deletefileiterator = glob.iglob(os.path.join(experimentfolder,
             '*.req'))
        # delete req files
        for filename in deletefileiterator:
            os.remove(filename)
        for filetype in otherfilestodelete:
            filetype = '*/' + filetype
            deletefileiterator = glob.iglob(os.path.join(experimentfolder,
                                                filetype))
            for filename in deletefileiterator:
                os.remove(filename)
    if singlefile:
        for filetype in otherfilestodelete:
            deletefileiterator = glob.iglob(os.path.join(experimentfolder,
                                                filetype))
            for filename in deletefileiterator:
                os.remove(filename)

def remove_all_node_files(exppath):
    """2013-12-19 12:21 WEV
    This script is made to complement the hyak_processing script
    it will remove all of the node files that are generated for
    multiple node processing of data
    """
    nodelist = [ls for ls in os.listdir(exppath) \
                if os.path.isdir(os.path.join(exppath, ls)) and
                    'node' in ls]
    for node in nodelist:
        shutil.rmtree(os.path.join(exppath, node))
