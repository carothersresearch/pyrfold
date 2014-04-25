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

class TimeCourseStructure:
    """
    class will contain the information of a given set of experimental runs on
    Kinefold

    Requires an output dictionay structure that is output from hyak processing
    """
    def __init__(self, dictionaryofruns=None, indexestomine=None, timescale=None,
                    cutoff=0.0, rescale=False, maxlength=None):

        (self.dictionary, self.baseadditiontime, self.__sequence) = \
                         self._compensate_for_intial_translation(dictionaryofruns)
        #A list of structures that have been found thus far and their max
        #frequency
        self.structures = {}
        self.timedata = {}
        self.stats = {}
        if maxlength:
            self.indexestomine = _calculate_indexs_from_time(
                                        self.dictionary, timescale)
            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                self.indexestomine, timescale, cutoff, rescale)
            self.sequence = self.__sequence[
                                        indexestomine[0] - 1: indexestomine[1]]

        elif indexestomine:
            if indexestomine[1] > len(self.__sequence):
                indexestomine[1] = len(self.__sequence)
                print 'readjusted max index'
            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                indexestomine, timescale, cutoff, rescale)
            self.sequence = self.__sequence[
                            self.indexestomine[0] - 1: self.indexestomine[1]]

    def generate_data(self, indexestomine=None, timescale=None, cutoff=0.0,
                        rescale=False, maxlength=None):
        # if indexestomine[1] > len(self.__sequence):
        #     indexestomine[1] = len(self.__sequence)
        #     print 'readjusted max index'
        if maxlength:
            self.indexestomine = _calculate_indexs_from_time(
                                        self.dictionary, timescale)

            (self.timedata, self.structures, self.stats) = \
                _structure_evaluation(self.dictionary,
                self.indexestomine, timescale, cutoff, rescale)
            self.sequence = self.__sequence[
                                self.indexestomine[0] - 1: self.indexestomine[1]]
        else:
            (self.timedata, self.structures, self.stats) = \
                    _structure_evaluation(self.dictionary,
                    indexestomine, timescale, cutoff, rescale)
            self.sequence = self.__sequence[
                                        indexestomine[0] - 1: indexestomine[1]]

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

    def _compensate_for_intial_translation(self, dictionaryofruns):
        "All of these experiments are not condering the initial time of eleongation"
        #We need to identify the polymerization rate
        newdict = {}
        timeofbaseaddition = None
        for runnumber in dictionaryofruns:
            dotbracketlist = deepcopy(dictionaryofruns[runnumber]['dotbracket'])
            timelist = deepcopy(dictionaryofruns[runnumber]['time'])
            if not dotbracketlist:
                continue
            elif len(dotbracketlist[0]) == 1:
                continue
            newdict[runnumber] = {}
            #Polymerization rate
            timeofbaseaddition = []
            for counter, dotbracket in enumerate(dotbracketlist):
                currentlength = len(dotbracket)
                shift = 1
                numberfound = 0
                while counter-shift >= 0:
                    if len(dotbracketlist[counter - shift]) == currentlength:
                        shift += 1
                    else:
                        timeofbaseaddition.append(timelist[counter] - timelist[counter-shift])
                        numberfound += 1
                        break
                if len(timeofbaseaddition) == 7:
                    timeofbaseaddition = np.mean(timeofbaseaddition)
                    break
            #Now have to add time and dotbrackets to this system
            newtimelist = []
            newdotbracketlist = []
            for basenumber in range(len(dotbracketlist[0]) - 1):
                newtimelist.append(timeofbaseaddition * basenumber)
                newdotbracketlist.append('.' + '.'*basenumber)
            basetime = max(newtimelist) + timeofbaseaddition
            for counter, dotbracket in enumerate(dotbracketlist):
                newdotbracketlist.append(dotbracket)
                newtimelist.append(basetime + timelist[counter])
            newdict[runnumber]['dotbracket'] = deepcopy(newdotbracketlist)
            newdict[runnumber]['time'] = deepcopy(newtimelist)
            sequence = dictionaryofruns[runnumber]['sequence'][-1]
        #print newdict
        return newdict, timeofbaseaddition, sequence

def _calculate_indexs_from_time(dictionaryofruns, timescale):
    """ This serves to find the last possible base of the window
    sequence """
    output = [1, 0]
    maxtime = timescale[1]
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

def _structure_evaluation(dictionaryofruns, indexstomine,
                                    timescale=None, cutoff=0., rescale=False):
    """
    :param dictionaryofruns: This is the standard dictionary that is output by
        Kinefold runnumber:['dotbracket', etc, etc,]:list
    :type dictionaryofruns: Dictionary
    :param indexestomine: The start and stop positions of the dotbrackets to
        mine over time (1 indexed)
    :type inddexestomine: list of lists of ints
        NOTE: This is currently broken and only supports a single index
    :param timescale: If used this outlines the timepoints of interest
    :type timescale: list of lists of floats
    :param cutoff: This is a way to remove not common structures - anystructure
        which doesn't appear more than the cutoff is removed from the output
    :type cutoff: float
    """
    #Probably best to have this set up as a dictionary output?
    tempdict = {} #keys: structures values: [count, time]
    maxtime = 0
    indexs = indexstomine
    #print dictionaryofruns[1]
    for runnumber in dictionaryofruns:
        #print 'run number'
        if not dictionaryofruns[runnumber]['time']:
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
    #print tempdict
    timelist = _calculate_time_list(tempdict, timescale, rescale)
    # print timelist
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
    #             #print outdict
    # print outdict
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
    timecoursedatadict = {x:y for (x, y) in
                            outdict.items() if x in structstoouput}
    # print timecoursedatadict
    structures, stats = calculate_stats(timecoursedatadict)
    return timecoursedatadict, structures, stats

def _calculate_time_list(dictofruns, timescale, rescale=False, minstep=75):
    templist = [0]
    for run in dictofruns:
        templist.extend(dictofruns[run]['time'])
        #print dictofruns[run]['time']
    templist = set(templist)
    templist = list(templist)
    templist.sort()
    #print templist
    if timescale:
        maxtime = timescale[1]
        mintime = timescale[0]
    else:
        maxtime = max(templist)
        mintime = min(templist)
    if rescale and maxtime > max(templist):
        #asking for data outside of timescale
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

    return outlist

def calculate_stats(timecoursedictionary):
    structs = {}
    stats = {}
    for structure in timecoursedictionary:
        timecountarray = np.array(timecoursedictionary[structure])
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
    return structs, stats

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
        #To account for structures that fall in the cracks
        # elif counter < sizeofdictionary - 1:
        #     if (timecount[0] < timepoints[0] < dictionary[structure][counter + 1][0]
        #         and timecount[0] < timepoints[1] < dictionary[structure][counter + 1][0]):
        #         dictionary[structure][counter][1] += 1
    return dictionary
