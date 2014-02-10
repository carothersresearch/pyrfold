"""This module exists to process all of the kinefold output data
into a form which is actually useful

Things TODO
- Add psudeonot conversion tools
"""

########################################################################
############################ Modules ###################################
########################################################################
import os
import glob
import cPickle as pickle
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
                            singledirectory=False, csvoutput=False):
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
    if csvoutput == False:
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
                tempdictionary['kine'] = []
                tempdictionary['sequence'] = []
                tempdictionary['base count'] = []
                tempdictionary['kine helix'] = []
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
                            tempdictionary['kine'].append(kineobj.structure)
                            tempdictionary['base count'].append(
                                kineobj.basenumber)
                            tempdictionary['sequence'].append(kineobj.sequence)
                            tempdictionary['kine helix'].append(kineobj.helix)
                            tempdictionary['energy'].append(kineobj.freeenergy)
                            tempdictionary['type'] = kineobj.type
                            tempdictionary['total bases'] = kineobj.basetotal
                dicttopickle[rnmnumber] = tempdictionary
            #Now to write the output for a single experiment
            temppickledest = os.path.join(outputfilename, directoryname + '.p')
            pickle.dump(dicttopickle, open(temppickledest, 'wb'))
    """
    if csvoutput:
        root = os.getcwd()
        if not os.path.exists(os.path.join(
                                    processeddirectory, rawoutputfolder)):
            os.makedirs(os.path.join(processeddirectory, rawoutputfolder))
        files = os.walk(root).next()[1]
        #allows for rules used by the unix shell to find files
        for f in files:
            if not os.path.exists(os.path.join(processeddirectory,
                                                        rawoutputfolder, f)):
                # creates subdir
                os.mkdir(os.path.join(processeddirectory, rawoutputfolder, f))
            else:
                continue
            os.chdir(os.path.join(root, f)) #enters experiment specific folder
            # rnmslist = glob.glob('*.rnms')
            rnmlist = glob.glob('*.rnm')
            for rnm in rnmlist: #goes through every file
                rnmname = os.path.splitext(rnm)[0] #removes .extension

                #open rnm file
                rnmfile = open(rnm, 'r')
                # These are all of the things
                listofkineobjects = []
                for i, line in enumerate(rnmfile):
                    if i in [0, 1]:
                        continue
                    elif i % 2 == 0:
                        tempobject = KineRmnOutput(line)
                    elif (i+1) % 2 == 0: #should account for odds
                        tempobject.addHelix(line)
                        tempobject.generateDotBracket()
                        listofkineobjects.append(tempobject)
                rnmfile.close()

                output = open(os.path.join(processeddirectory,
                 rawoutputfolder, f, rnmname + '.csv'), 'wb')
                #write header
                headerlist = ['structure', 'helix numbering',
                'free energy (kcal/mol)', 'time (ms)', 'base #', 'Base total',
                 'sequence', 'dotbracket']
                for head in headerlist:
                    output.write(head)
                    output.write(',')
                output.write('\n')

                for kineob in listofkineobjects:
                    output.write(kineob.structure + ',')
                    output.write(kineob.helix + ',')
                    output.write(kineob.freeenergy + ',')
                    output.write(kineob.time + ',')
                    output.write(kineob.basenumber + ',')
                    output.write(kineob.basetotal + ',')
                    output.write(kineob.sequence + ',')
                    output.write(kineob.dotbracket + ',')
                    output.write('\n')
                output.close()
        os.chdir(root)
        """

########################################################################
################# Gen.  Destruction  ###################################
########################################################################
def clean_up_output_data(experimentfolder, singlefile=False):
    """2014-01-08 17:12 WEV
    Removes all of the undesired files from the output of kinefold
    """
    otherfilestodelete = ['*/*.rnms', '*/*.rnm2','*/*.e'
                                    ,'*/*.p', '*/*.i']
    if not singlefile:
        experimentfolder = os.path.join(experimentfolder, 'output')
        #change to this path
        deletefileiterator = glob.iglob(os.path.join(experimentfolder,
             '*.req'))
        # delete req files
        for filename in deletefileiterator:
            os.remove(filename)
        for filetype in otherfilestodelete:
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
