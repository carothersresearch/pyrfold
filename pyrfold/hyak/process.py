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
from collections import Counter, OrderedDict
from . _convert import KineRmnOutput
import errno
import numpy as np
import shutil


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


def timecourse(rawoutputfolder, processeddirectory, subsummarydict=None,
               singledirectory=False, return_dictionary=False):
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
    # Make the directory to write the files to
    TIMECOURSEDIR = os.path.join(processeddirectory, 'timecourse')
    COMPRESSEDTIMEDIR = os.path.join(processeddirectory, 'compressedtime')
    make_sure_path_exists(TIMECOURSEDIR)
    make_sure_path_exists(COMPRESSEDTIMEDIR)
    # collect all of the experiment files
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
            rnmnumber = int(rnmname[0:3])  # should pull just the runnumber
            # open the file for reading
            tempdictionary = {}
            tempdictionary['time'] = []
            tempdictionary['dotbracket'] = []
            sequencelist = []
            tempdictionary['energy'] = []
            with open(rnmfile, 'r') as rnmopenfile:
                for i, line in enumerate(rnmopenfile):
                    if i in [0, 1]:
                        continue
                    elif i % 2 == 0:
                        kineobj = KineRmnOutput(line)
                    elif (i+1) % 2 == 0:  # should account for odds
                        kineobj.addHelix(line)
                        kineobj.generateDotBracket()
                        tempdictionary['time'].append(kineobj.time)
                        tempdictionary['dotbracket'].append(
                            kineobj.dotbracket)
                        sequencelist.append(kineobj.sequence)
                        tempdictionary['energy'].append(kineobj.freeenergy)
                        tempdictionary['type'] = kineobj.type
                        tempdictionary['total bases'] = kineobj.basetotal
            try:
                tempdictionary['sequence'] = sequencelist[-1]
                dicttopickle[rnmnumber] = tempdictionary
            except KeyError:
                print rnmfile
                print 'KeyError, probably did not finish a simulation'
            except IndexError:
                print rnmfile
                print 'IndexError, probably did not finish a simulation'
        # Now to write the complete timecourse data
        temppickledest = os.path.join(TIMECOURSEDIR, directoryname + '.p')
        if not return_dictionary:
            with open(temppickledest, 'wb') as topick:
                pickle.dump(dicttopickle, topick, protocol=2)
        # Now to generate the compressed version of this dictionary
        if subsummarydict:
            polrate = subsummarydict[directoryname].kine_folding_data()[0]
            polrate = np.round(polrate, 1)
        else:
            polrate = None
        temprundict, baseadditiontime, completesequence, simulation_type = \
            consolidate_run_dictionary(dicttopickle, polrate=polrate)
        compresseddict = compress_run_dictionary(temprundict,
                                                 baseadditiontime,
                                                 completesequence,
                                                 simulation_type)
        temppickledest = os.path.join(COMPRESSEDTIMEDIR, directoryname + '.p')
        if return_dictionary:
            return dicttopickle, compresseddict
        with open(temppickledest, 'wb') as topick:
            pickle.dump(compresseddict, topick, protocol=2)


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
    # Sort the paths
    pathtotimecourse = os.path.join(processeddirectory, 'timecourse')
    pathtofinalstructures = os.path.join(processeddirectory, 'finalstructure')
    if not os.path.exists(pathtofinalstructures):
        make_sure_path_exists(pathtofinalstructures)
    # collect all of the experiment files
    if not singledirectoryname:
        pickletimelist = [ls for ls in os.listdir(pathtotimecourse)
                          if os.path.isfile(os.path.join(pathtotimecourse, ls))
                          ]
    else:
        pickletimelist = [singledirectoryname]
    # Now to go through all of these files
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
            # Convert count to average
            totalcount = float(sum(dotbracketcount.viewvalues()))
            for dotbrac in dotbracketcount:
                dotbracketcount[dotbrac] /= totalcount
            tempdict['parts'][part]['sequence'] = sequencecount.keys()[0]
            tempdict['parts'][part]['dotbracket'] = dotbracketcount
    return tempdict


###############################################################################
###############################################################################
# CONSOLIDATE RUN DICTIONARY BLOCK START
###############################################################################
###############################################################################
def consolidate_run_dictionary(rundictionary, polrate=None):
    """ This function will distill a run dictionary into a very basic dictionary
    which will be a single set of time points to all of the structures that
    are observed at various times
    :param rundictionary: Dictionary that hyak processing produces
    """
    #First calcualte the pol rate and base timeline
    timebofbaseadition, timeline, sequence, simulation_type = calculate_pol_rate(rundictionary,
                                                        polrate)
    #we do this because we want all of the structural shifts to be discrete
    rundictionary = change_duplicate_time_points(rundictionary)
    #Rescale all of the time vectors in the rundictionary to align with timelin
    if simulation_type == 'cotrans':
        rundictionary = rescale_time_vectors(rundictionary, timeline)
        rundictionary = change_duplicate_time_points(rundictionary)
    #Do a final polishing step to make sure that there are not duplicate time points
    return rundictionary, timebofbaseadition, sequence, simulation_type


def change_duplicate_time_points(dictionaryofruns):
    """This will scan through all of the runs and add a small
    amount of time to timepoints which are reported as being
    identical"""
    for runnumber in dictionaryofruns:
        timelist = dictionaryofruns[runnumber]['time']
        previoustime = -0.1
        for index, time in enumerate(timelist):
            currenttime = time
            if previoustime == currenttime:
                # need to shift current time up
                timelist[index] += 0.1
            previoustime = currenttime
        # dictionaryofruns[runnumber]['time'] = TIMELIST
    return dictionaryofruns


def calculate_pol_rate(dictionaryofruns, polrate=None):
    """ function scans through run dictionary to calcuate pol rate if this
    value is not contained in the pickled structure
    :param dictionaryofruns: Silly dictionary run:dotrbacket etc....
    :type dictionaryofruns: dict
    :param polrate: The polrate in ms/nt
    :type polrate: float
    :return timeofbaseaddition: ms/nt
    :return timeline: Timeline of base addition. EVery time point a base was added
    :return timeofbaseaddition: list
    """
    simulation_type = 'cotrans'

    if polrate is not None:
        timeofbaseaddition = np.round(polrate, 1)
        for runnumber in dictionaryofruns:
            if dictionaryofruns[runnumber]['sequence']:
                sequence = dictionaryofruns[runnumber]['sequence']
                break
    else:
        timeofbaseaddition = []
        #Polymerization rate
        #print dictionaryofruns
        for runnumber in dictionaryofruns:
            #Check to see if the simulation is meltanneal
            if dictionaryofruns[runnumber]['type'] == 'meltanneal':
                simulation_type = 'meltanneal'
                sequence = dictionaryofruns[runnumber]['sequence']
                polrate = None
                break
            tempbaseadditionlist = []
            dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
            timelist = dictionaryofruns[runnumber]['time']
            if dictionaryofruns[runnumber]['sequence']:
                sequence = dictionaryofruns[runnumber]['sequence']
            first_iteration = True
            for counter, dotbracket in enumerate(dotbracketlist):
                if first_iteration:
                    first_iteration = False
                    currentlength = len(dotbracket)
                else:
                    if len(dotbracket) == currentlength:
                        continue
                    currentlength = len(dotbracket)
                shift = 1
                while counter-shift > 0:
                    if len(dotbracketlist[counter - shift]) == currentlength:
                        shift += 1
                    else:
                        tempbaseadditionlist.append(timelist[counter] - timelist[counter-shift])
                        break
                if len(tempbaseadditionlist) == 8:
                    timeofbaseaddition.append(np.mean(tempbaseadditionlist))
                    break
        timeofbaseaddition = np.round(np.mean(timeofbaseaddition), 1)
        #Make the timeline

    if simulation_type == 'cotrans':
        timeline = []
        for basecount in range(len(sequence)):
            timeline.append(basecount * timeofbaseaddition)
    else:
        timeline = None
    return timeofbaseaddition, timeline, sequence, simulation_type


def rescale_time_vectors(dictionaryofruns, timeline):
    """This function needs to rescale all time vectors to align with the
    new timeline that was generated
    """
    #Cycle through all of the time vectors
    for runnumber in dictionaryofruns:
        timelist = dictionaryofruns[runnumber]['time']
        dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
        dictionaryofruns[runnumber]['time'] = rescale_time_vector(
                                        timelist, dotbracketlist, timeline)
        #length of the sequence at first base addition
    return dictionaryofruns


def rescale_time_vector(timelist, dotbracketlist, timeline):
    totalsize = len(timeline)
    sizeoftimelist = len(timelist)
    countstart = 0
    presize = 0
    for count, dotbracket in enumerate(dotbracketlist):
        cursize = len(dotbracket)
        #Second condition accounts for first pass
#         print 'presize - ' + str(presize)
#         print 'pretime - ' + str(timelist[countstart])
#         print 'cursize - ' + str(cursize)
#         print 'curtime - ' + str(timelist[count])
#         raw_input()
        if cursize > presize:
            if presize == 0:
                presize = cursize
                countstart = count
                continue
            #update presize
            adjust_time_window(timelist, countstart, count, presize, timeline)
            presize = cursize
            countstart = count
        #check to see if there are any more addition steps
        if presize == totalsize:

            diff_last_add = timeline[-1] - timelist[count]
#             print diff_last_add
#             print timelist[count - 1:]
            timelist[count] = timeline[-1]
            for new_count in range(count + 1, sizeoftimelist):
#                 print timelist[new_count]
                timelist[new_count] += diff_last_add
#                 print timelist[new_count]
            break
    return np.around(timelist, decimals=1)


def adjust_time_window(timelist, indexstart, indexstop, size, timeline):
    basetime = timeline[size - 1]
    #If there is only one step between no interpoliation is neccesary
    if indexstart + 1 == indexstop:
        timelist[indexstart] = basetime
        return timelist
    ortinalstarttime = timelist[indexstart]
    orginalstoptime = timelist[indexstop]
    timelist[indexstart] = basetime
    newstart = timeline[size -1]
    newstop = timeline[size]
    for index in range(indexstart + 1, indexstop):
        timelist[index] = adjust_time_point(timelist[index], ortinalstarttime, orginalstoptime, newstart, newstop)


def adjust_time_point(orginaltime, orgstart, orgstop, newstart, newstop):
    """simple lever arm method to determine the new value of the time point
    """
    a = (orginaltime - orgstart)
    b = (orgstop - orgstart)
    value = a/b*(newstop - newstart) + newstart
    return value

###############################################################################
###############################################################################
# CONSOLIDATE RUN DICTIONARY BLOCK STOP
###############################################################################
###############################################################################

###############################################################################
###############################################################################
# COMPRESS RUN DICTIONARY BLOCK START
###############################################################################
###############################################################################


def compress_run_dictionary(rundictionary, baseadditiontime, completesequence,
                            simulation_type):
    dotdict = OrderedDict()
    energydict = {}
    timearray = calculate_time_list(rundictionary, baseadditiontime,
                                    simulation_type)
    #Initialize the dictionary
    for time in timearray:
        dotdict[time] = Counter()
    maxtime = timearray[-1]
    for runnumber in rundictionary:
        #print runnumber
        temptimelist = rundictionary[runnumber]['time']
        totalnumber = len(temptimelist)
        #add unstructured sequences structures where appropriate
        mintime = temptimelist[0]
        # numberofsinglebasestoadd = int(mintime/baseadditiontime)
        previousbasecountindex = 0
        if simulation_type == 'cotrans':
            numberofsinglebasestoadd = len(rundictionary[runnumber]['dotbracket'][0])

            for basecount in range(numberofsinglebasestoadd):
                temptime1 = find_nearest(timearray, basecount * baseadditiontime)
                temptime2 = find_nearest(timearray,
                                        (basecount + 1) * baseadditiontime)
                structure = '.'*(basecount + 1)
                dotdict, previousbasecountindex = add_structure_timedictionary(dotdict, structure,
                            [temptime1, temptime2], previousbasecountindex)
                add_energy_data(energydict, structure, energy=0)

        #print 'initial stuff added'
        #add all of the other structures
        for counter, time in enumerate(temptimelist):
            structure = rundictionary[runnumber]['dotbracket'][counter]
            currenttime = time
            energy = rundictionary[runnumber]['energy'][counter]
            if counter < totalnumber - 1:
                nexttime = temptimelist[counter + 1]
            else:
                nexttime = maxtime + 1
            dotdict, previousbasecountindex  = add_structure_timedictionary(dotdict, structure,
                            [currenttime, nexttime], previousbasecountindex)
            energydict = add_energy_data(energydict, structure, energy)
    outdict  = {}
    dotdict = normailize_orderddict_counters(dotdict)
    outdict['dotbracket'] = dotdict
    outdict['energy'] = energydict
    outdict['sequence'] = completesequence
    outdict['baseadditiontime'] = baseadditiontime
    return outdict


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def normailize_orderddict_counters(ordereddict):
    """adding hacky removal of not used bases"""
    for time in ordereddict:
        if not ordereddict[time]:
            print time
        tempcount = sum(ordereddict[time].values())
        for item in ordereddict[time]:
            ordereddict[time][item] /= float(tempcount)
    return ordereddict


def add_structure_timedictionary(ordereddict, structure, times, previousstart):
    """Previous start is to increase the efficency of this process"""
    timelist = ordereddict.keys()
    for index, time in enumerate(timelist[previousstart:]):
        if time >= times[0]:
            if time >= times[1]:
                break
            try:
                ordereddict[time][structure] += 1
            except KeyError:
                ordereddict[time][structure] = 1
    previousstart += index - 1
    return ordereddict, previousstart


def add_energy_data(energydict, structure, energy):
    if structure not in energydict:
        energydict[structure] = energy
    return energydict


def calculate_time_list(dictofruns, baseaddition, simulation_type):
    templist = []
    for run in dictofruns:
        templist.extend(dictofruns[run]['time'])
        #print dictofruns[run]['time']
    templist = set(templist)
    templist = list(templist)
    templist.sort()
    #print templist
    maxtime = max(templist)
    mintime = min(templist)
    #Add timesteps for initial bass addition
    if simulation_type == 'cotrans':
        timepointstoadd = np.arange(0, mintime, baseaddition)
        outlist = list(timepointstoadd)
    else:
        outlist = []
    previousstep = -0.1
    for time in templist:
        if time > previousstep:
            outlist.append(time)\
            previousstep = time
    return np.array(outlist)

###############################################################################
###############################################################################
# COMPRESS RUN DICTIONARY BLOCK stop
###############################################################################
###############################################################################

########################################################################
################# Gen.  Destruction  ###################################
########################################################################


def clean_up_output_data(experimentfolder, singlefile=False):
    """2014-01-08 17:12 WEV
    Removes all of the undesired files from the output of kinefold
    """
    otherfilestodelete = ['*.rnms', '*.rnm2', '*.e',
                          '*.p']  # , '*.i']  #, '*.rnm']
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


def delete_files_in_folder(directory, name_to_look_for):
    """This function is delete all of the left over dat files from processing
    """
    files_to_delete = [os.path.join(directory, ls) for ls in os.listdir(directory)
                       if name_to_look_for in ls]
    for file_name in files_to_delete:
        os.remove(file_name)


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

########################################################################
################# Compression        ###################################
########################################################################


def compress_and_delete_directory(directorypath, tarpath='auto'):
    import tarfile
    import os
    import shutil
    if tarpath == 'auto':
        tarpath = os.path.dirname(os.path.abspath(directorypath))
    nameofdir = os.path.basename(os.path.abspath(directorypath))
    nameoftar = nameofdir + '.tar.gz'
    tardir = os.path.join(tarpath, nameoftar)
    with tarfile.open(tardir, "w:gz") as tar:
        tar.add(directorypath, arcname=os.path.basename(directorypath))
    shutil.rmtree(directorypath)

