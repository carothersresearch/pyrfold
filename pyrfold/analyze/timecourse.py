"""
This module exists to do basic processing of timecourse data that is output
from kinefld simulations


Possible things TODO:
add a wildcardstructure thing like ****** when structures not yet exposed

"""

import math
import numpy as np
from scipy.signal import argrelextrema
from copy import deepcopy
from collections import OrderedDict
from collections import Counter

class TimeCourseStructure:
    """
    class will contain the information of a given set of experimental runs on
    Kinefold

    Requires an output dictionay structure that is output from hyak processing
    """
    def __init__(self, dictionaryofruns, structurewindow=None, timewindow=None,
                        rescale=False, maxlength=False, cutoff=0.0):
        """
        This function initializes the object - the only thing required for
        this step is the dictionaryofruns
        :param dictionaryofruns: This is the dictionary containing all of the
            simulation imformation
        :type dictionaryofruns: dict
        :param structurewindow: The start and stop window that you want to look
            at - (this is 1 indexed!)
        :type structurewindow: list
        :param timewindow: The start and stop of the timewindow you want to
            consider
        :type timewindow: list
        :param cutoff: The smallest observed max frequency of structure to
            consider. This is great for reducing the size of things that you
            produce
        :type cutoff: float
        :param rescale: This will extend the max time of a simulation if True.
            this is important for considering a diverse set of simulations
        :type rescale: boolean
        :param maxlength: If true the maximum length of the sequence will be
            used as the end of the indexstomine
        :type maxlength: boolean
        """
        (self.dictionary, self.baseadditiontime, self.__sequence) = \
                         self._compensate_for_intial_translation(dictionaryofruns)
        #A list of structures that have been found thus far and their max
        #frequency
        self.structures = {}
        self.timedata = {}
        self.stats = {}
        if maxlength:
            structurewindow = _calculate_indexs_from_time(
                                        self.dictionary, timewindow)
            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                structurewindow, timewindow, cutoff, rescale)
            self.sequence = self.__sequence[
                        structurewindow[0] - 1: structurewindow[1]]
        elif structurewindow:
            if structurewindow[1] > len(self.__sequence):
                structurewindow[1] = len(self.__sequence)
                print 'Requested max structure index was readjusted'
            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                structurewindow, timewindow, cutoff, rescale)
            self.sequence = self.__sequence[
                            structurewindow[0] - 1: structurewindow[1]]

    def generate_data(self, structurewindow, timewindow=None, cutoff=0.0,
                        rescale=False, maxlength=None):
        if maxlength:
            self.structurewindow = _calculate_indexs_from_time(
                                        self.dictionary, timewindow)

            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                self.structurewindow, timewindow, cutoff, rescale)
            self.sequence = self.__sequence[
                          self.structurewindow[0] - 1: self.structurewindow[1]]
        else:
            (self.timedata, self.structures, self.stats) = \
                    _structure_evaluation(self.dictionary,
                    structurewindow, timewindow, cutoff, rescale)
            self.sequence = self.__sequence[
                                    structurewindow[0] - 1: structurewindow[1]]

    def generate_time_to_structure_dict(self):
        outdict = {}
        for structure in self.timedata:
            for timefreq in self.timedata[structure]:
                if timefreq[0] in outdict:
                    outdict[timefreq[0]][structure] = timefreq[1]
                else:
                    outdict[timefreq[0]] = {}
                    outdict[timefreq[0]][structure] = timefreq[1]
        self.timetostructures = outdict
        #return outdict

    def __repr__(self):
        #Crude Sorting system
        sortlist = []
        for structure in self.structures:
            sortlist.append([self.structures[structure][0], structure])
        sortlist.sort()
        output = ''
        output += '======================================================'
        output += '\n'
        output += self.sequence + '   ' + 'Time (ms)' + '   ' + 'Max Freq'
        for structure in sortlist:
            output += '\n'
            structure = structure[1]
            linetoprint = structure
            linetoprint += '   '
            linetoprint += str(self.structures[structure][0]).zfill(8)
            linetoprint += '         '
            linetoprint += str(self.structures[structure][1])
            output += linetoprint
        return output

    def generate_compressed_dictionary(self):
        dotdict = OrderedDict()
        energydict = {}
        timelist = _calculate_time_list(self.dictionary)
        baseaddition = self.baseadditiontime
        #Rescale all of the time vectors

        #Try to account of numerical disagreement
        #Tolerance is the cutoff change in time for investigating if sizes of
        #The constrcuts are the same
        # tolerance = baseaddition/5.
        # pasttime = -50
        # for time in timelist:
        #     if abs(pasttime - time) < tolerance:
        #         #Need to investigate if this additional time point is due to a single base addition

        #Initialize the dictionary
        for time in timelist:
            dotdict[time] = Counter()
        maxtime = max(timelist)
        for runnumber in self.dictionary:
            totalnumber = len(self.dictionary[runnumber]['time'])
            for counter, time in enumerate(self.dictionary[runnumber]['time']):
                structure = self.dictionary[runnumber]['dotbracket'][counter]
                currenttime = time
                energy = self.dictionary[runnumber]['energy'][counter]
                if counter < totalnumber - 1:
                    nexttime = self.dictionary[runnumber]['time'][counter + 1]
                else:
                    nextime = maxtime + 1
                dotdict = _add_structure_timedictionary(dotdict, structure,
                                [currenttime, nexttime])
                energydict = _add_energy_data(energydict, structure, energy)
        outdict  = {}
        outdict['dotbracket'] = dotdict
        outdict['energy'] = energydict
        outdict['sequence'] = self.__sequence
        return outdict

    def _compensate_for_intial_translation(self, dictionaryofruns):
        """Kinefold takes timepoint 0 to be the time of the first helix
        formation. This function accounts for the intial polymerization that
        has to occur beforehand
         """
        #We need to identify the polymerization rate
        newdict = {}
        timeofbaseaddition = []
        #Polymerization rate
        for runnumber in dictionaryofruns:
            tempbaseadditionlist = []
            dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
            timelist = dictionaryofruns[runnumber]['time']
            for counter, dotbracket in enumerate(dotbracketlist):
                currentlength = len(dotbracket)
                shift = 1
                while counter-shift >= 0:
                    if len(dotbracketlist[counter - shift]) == currentlength:
                        shift += 1
                    else:
                        tempbaseadditionlist.append(timelist[counter] - timelist[counter-shift])
                        break
                if len(timeofbaseaddition) == 8:
                    timeofbaseaddition.append(tempbaseadditionlist)
                    break
        timeofbaseaddition = np.round(np.mean(tempbaseadditionlist), 2)
        for runnumber in dictionaryofruns:
            dotbracketlist = deepcopy(dictionaryofruns[runnumber]['dotbracket'])
            timelist = deepcopy(dictionaryofruns[runnumber]['time'])
            if not dotbracketlist:
                continue
            elif len(dotbracketlist[0]) == 1:
                sequence = dictionaryofruns[runnumber]['sequence']
                continue
            newdict[runnumber] = {}
            #Now have to add time and dotbrackets to this system
            newtimelist = []
            newdotbracketlist = []
            newenegylist = []
            for basenumber in range(len(dotbracketlist[0]) - 1):
                newtimelist.append(timeofbaseaddition * basenumber)
                newdotbracketlist.append('.' + '.'*basenumber)
                newenegylist.append(0)
            basetime = max(newtimelist) + timeofbaseaddition
            for counter, dotbracket in enumerate(dotbracketlist):
                newdotbracketlist.append(dotbracket)
                newtimelist.append(basetime + timelist[counter])
            newenegylist.extend(dictionaryofruns[runnumber]['energy'])
            newdict[runnumber]['dotbracket'] = deepcopy(newdotbracketlist)
            newdict[runnumber]['time'] = \
                                    np.round(np.array(deepcopy(newtimelist)), 3)
            newdict[runnumber]['energy'] = np.array(deepcopy(newenegylist))
            newdict[runnumber]['sequence'] = dictionaryofruns
            sequence = dictionaryofruns[runnumber]['sequence']
        #print newdict
        if not newdict:
            newdict = dictionaryofruns
        return newdict, timeofbaseaddition, sequence

def _calculate_indexs_from_time(dictionaryofruns, timewindow):
    """ This serves to find the last possible base of the window
    sequence """
    output = [1, 0]
    maxtime = timewindow[1]
    possiblemaxindexs = []
    for run in dictionaryofruns:
        for counter, timepoint in enumerate(dictionaryofruns[run]['time']):
            maxlength = len(dictionaryofruns[run]['dotbracket'][-1])
            if maxtime < timepoint:
                possiblemaxindexs.append(
                    len(dictionaryofruns[run]['dotbracket'][counter - 1]))
                break
    if not possiblemaxindexs:
        possiblemaxindexs.append(maxlength)
    maxbase = min(possiblemaxindexs)
    output[1] = maxbase
    return output

def _add_energy_data(energydict, structure, energy):
    if structure not in energydict:
        energydict[structure] = energy
    return energydict

def _add_structure_timedictionary(ordereddict, structure, times):
    for time in ordereddict:
        if (time >= times[0] and time < times[1]):
            if structure in ordereddict[time]:
                ordereddict[time][structure] += 1
            else:
                ordereddict[time][structure] = 1
    return ordereddict

def _structure_evaluation(dictionaryofruns, structurewindow,
                                    timewindow=None, cutoff=0., rescale=False):
    """
    :param dictionaryofruns: This is the standard dictionary that is output by
        Kinefold runnumber:['dotbracket', etc, etc,]:list
    :type dictionaryofruns: Dictionary
    :param structurewindow: The start and stop positions of the dotbrackets to
        mine over time (1 indexed)
    :type inddexestomine: list of lists of ints
        NOTE: This is currently broken and only supports a single index
    :param timewindow: If used this outlines the timepoints of interest
    :type timewindow: list of lists of floats
    :param cutoff: This is a way to remove not common structures - anystructure
        which doesn't appear more than the cutoff is removed from the output
    :type cutoff: float
    """
    #Probably best to have this set up as a dictionary output?
    tempdict = {} #keys: structures values: [count, time]
    maxtime = 0
    indexs = structurewindow
    for runnumber in dictionaryofruns:
        if not dictionaryofruns[runnumber]['time'].any():
            continue
        tempdict[runnumber] = {}
        tempdict[runnumber]['struct'] = []
        tempdict[runnumber]['time'] = []
        timelist = dictionaryofruns[runnumber]['time']
        if maxtime < timelist[-1]:
            maxtime = timelist[-1]
        dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
        #sequencelist = dictionaryofruns[runnumber]['sequence']
        #sequence = self.sequence[indexstomine[0] - 1: indexstomine[1]]
        for counter, dotbracket in enumerate(dotbracketlist):
            if len(dotbracket) >= indexs[1]:
                #start recording structure
                tempdict[runnumber]['struct'].append(
                                      dotbracket[indexs[0] - 1: indexs[1]])
                tempdict[runnumber]['time'].append(timelist[counter])
    timelist = _calculate_time_list(tempdict, timewindow, rescale)
    maxtime = max(timelist)
    mintime = min(timelist)
    #Find the total parts
    totalparts = len(tempdict)
    outdict = {}
    for runnumber in tempdict:
        # previousstructures = None
        for counter, structure in enumerate(tempdict[runnumber]['struct']):
            if counter == len(tempdict[runnumber]['time']) - 1:
                timepoints = [tempdict[runnumber]['time'][counter], maxtime]
            else:
                timepoints = tempdict[runnumber]['time'][counter:counter + 2]
            if timepoints[0] > maxtime:
                #print 'breaking'
                break
            elif timepoints[0] == timepoints[1]:
                #we don't want to double count these
                continue
            elif timepoints[0] < mintime:
                #We don't want to consider earlier timepoints
                #Check to see if this structure exists into the actual
                #Window
                if timepoints[1] > mintime:
                    timepoints = [mintime, timepoints[1]]
                else:
                    continue
            if structure in outdict:
                outdict = add_structures_to_timecountlist(
                                                outdict, structure, timepoints)
            else:
                outdict[structure] = [[ls, 0] for ls in timelist]
                # print structure
                # print timepoints
                outdict = add_structures_to_timecountlist(
                                                outdict, structure, timepoints)
    #Reduce the size of the dictionary
    #outdict = reduce_size_of_dict(outdict)
    structstoouput = []
    for structure in outdict:
        maxfreq = 0
        for count, timecount in enumerate(outdict[structure]):
            tempfreq = timecount[1]/float(totalparts)
            outdict[structure][count][1] = tempfreq
            if tempfreq > maxfreq:
                maxfreq = tempfreq
        if maxfreq >= cutoff:
            structstoouput.append(structure)
    timecoursedatadict = {x:np.array(y) for (x, y) in
                            outdict.items() if x in structstoouput}
    # print timecoursedatadict
    structures, stats, timecoursedatadict = calculate_stats(timecoursedatadict)
    return timecoursedatadict, structures, stats

def _calculate_time_list(dictofruns, timewindow=None, rescale=False,
        minstep=100):
    templist = [0]
    for run in dictofruns:
        templist.extend(dictofruns[run]['time'])
        #print dictofruns[run]['time']
    templist = set(templist)
    templist = list(templist)
    templist.sort()
    #print templist
    if timewindow:
        maxtime = timewindow[1]
        mintime = timewindow[0]
    else:
        maxtime = max(templist)
        mintime = min(templist)
    if rescale and maxtime > max(templist):
        #asking for data outside of timewindow
        entry = max(templist)
        while entry < maxtime:
            entry += minstep
            templist.append(entry)
    outlist = [mintime]
    previousstep = mintime
    for time in templist:
        if time >= maxtime:
            #outlist.append(maxtime - 0.005)
            outlist.append(maxtime)
            break
        elif time > previousstep:
            #outlist.append(time - 0.005)
            outlist.append(time)
            previousstep = time
        elif time > previousstep + minstep:
            outlist.append(time)
            previousstep = time

    return np.array(outlist)

def calculate_stats(timecoursedictionary):
    structs = {}
    stats = {}
    for structure in timecoursedictionary:
        timecountarray = timecoursedictionary[structure]
        #Find the first instance of it's max
        maxindex = np.argmax(timecountarray[:, 1])
        maxfeq = timecountarray[maxindex, 1]
        maxtime = timecountarray[maxindex, 0]
        structs[structure] = [maxtime, maxfeq]
        #finding the first appearance
        for counter, countofstru in enumerate(timecountarray[:, 1]):
            if countofstru > 0:
                timeoffirst = timecountarray[counter, 0]
                break
        #Finding all of the extrema
        indexofmaxes = local_maxima_index(timecountarray[:, 1])
        maxcountvalues = timecountarray[:, 1][indexofmaxes]
        timeofvalues = timecountarray[:, 0][indexofmaxes]
        stats[structure] = [timeoffirst, maxcountvalues, timeofvalues]
        #sorttime dict by time of first appearance
        structkeys = stats.keys()
        sortlist = zip([stats[struct][0] for struct in structkeys], structkeys)
        sortlist.sort()
        outtimedict = OrderedDict()
        for time, struct in sortlist:
            outtimedict[struct] = timecoursedictionary[struct]

    return structs, stats, outtimedict

def reduce_size_of_dict(timecoursedictionary):
    """
    This function will go through the time and frequency data and remove points
    Which provide no insight into how the structure changes with time.
    Basically it will only have time and structure values if there is a dynamic
    change.
    """
    outdict = {}
    minimumtimestep = 100
    previousstep = 0
    for structure in timecoursedictionary:
        outdict[structure] = []
        firstentry = False
        for counter, timefreq in enumerate(timecoursedictionary[structure]):
            if not firstentry:
                firstentry = True
                outdict[structure].append(timefreq)
                previousstep = timefreq[0]
                previousfreq = timefreq[1]
            #If there is a change in the value we need to add them to the list
            elif timefreq[1] != previousfreq:
                #append the previous value and the current value to appreciate
                #the curve
                outdict[structure].append(timecoursedictionary[structure][counter - 1])
                outdict[structure].append(timefreq)
                previousstep = timefreq[0]
                previousfreq = timefreq[1]
            #If we are past the minimum time step we need to add it to the list
            elif timefreq[0] > previousstep + minimumtimestep:
                outdict[structure].append(timefreq)
                previousstep = timefreq[0]
                previousfreq = timefreq[1]
    return outdict

def local_maxima_index(array):
    gradients = np.diff(array)
    leftgradeindex = 0
    listofmaxes = []
    foundincrease = True
    for count, grade in enumerate(gradients):
        if grade == 0:
            pass
        if grade < 0 and foundincrease:
            #Time to restart
            listofmaxes.append(leftgradeindex)
            foundincrease = False
        if grade > 0:
            leftgradeindex = count + 1
            foundincrease = True
    if foundincrease:
        listofmaxes.append(len(gradients))
    return listofmaxes


def add_structures_to_timecountlist(dictionary, structure, timepoints):
    """
    This will scan through the timepoints of the dictionary and add the
    structure and adjust the counts moving forward
    """
    sizeofdictionary = len(dictionary[structure])
    for counter, timecount in enumerate(dictionary[structure]):
        if timepoints[0] < timecount[0] < timepoints[1]:
            dictionary[structure][counter][1] += 1
        elif timepoints[0] == timecount[0]:
            dictionary[structure][counter][1] += 1
        elif timecount[0] > timepoints[1]:
            break
    return dictionary
