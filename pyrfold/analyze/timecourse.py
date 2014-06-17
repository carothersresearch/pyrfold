"""
This module exists to do basic processing of timecourse data that is output
from kinefld simulations


Possible things TODO:
add a wildcardstructure thing like ****** when structures not yet exposed

"""

import math
import numpy as np
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
from copy import deepcopy
from collections import OrderedDict
from collections import Counter

class TimeCourseStructure(object):
    """
    class will contain the information of a given set of experimental runs on
    Kinefold

    Requires an output dictionay structure that is output from hyak processing
    """
    def __init__(self, compresseddict, structurewindow=None, timewindow=None,
              rescale=False, firstexposure=False, maxlength=False, cutoff=0.0):
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
        self.structures = {}
        self.structuredataframe = None
        self.stats = {}
        self.timewindow = timewindow
        self.structurewindow = structurewindow
        #self.sequence = self.completesequence
        self.dictionary = compresseddict
        self.completesequence = compresseddict['sequence']

        if structurewindow or (timewindow and maxlength):
            self.generate_data(structurewindow, timewindow, cutoff, rescale,
                                                    maxlength, firstexposure)

    @classmethod
    def init_from_dictofruns(cls, dictionaryofruns, structurewindow=None, timewindow=None,
              rescale=False, firstexposure=False, maxlength=False, cutoff=0.0):
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
        temprundict, baseadditiontime, completesequence = \
                                   consolidate_run_dictionary(dictionaryofruns)
        #This step is time intensive and needs to be optimized
        dictionary = compress_run_dictionary(temprundict,
                                  baseadditiontime, completesequence)
        #A list of structures that have been found thus far and their max
        #frequency
        output = cls(dictionary, structurewindow, timewindow, rescale,
                                        firstexposure, maxlength, cutoff)
        return output


    def generate_data(self, structurewindow, timewindow=None, cutoff=0.0,
                        rescale=False, maxlength=None, firstexposure=False):

        self.structurewindow = structurewindow
        self.timewindow = timewindow
        if maxlength:
            self.structurewindow = calculate_indexs_from_time(
                                        self.dictionary, self.timewindow)
        self.structuredataframe = \
                structure_evaluation(self.dictionary, self.structurewindow,
                               self.timewindow, cutoff, rescale, firstexposure)
        self.sequence = self.completesequence[
                      self.structurewindow[0] - 1: self.structurewindow[1]]

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

    # def generate_compressed_dictionary(self):
    #     dotdict = OrderedDict()
    #     energydict = {}
    #     timelist = _calculate_time_list(self.dictionary)
    #     baseaddition = self.baseadditiontime
    #     #Initialize the dictionary
    #     for time in timelist:
    #         dotdict[time] = Counter()
    #     maxtime = max(timelist)
    #     for runnumber in self.dictionary:
    #         totalnumber = len(self.dictionary[runnumber]['time'])
    #         for counter, time in enumerate(self.dictionary[runnumber]['time']):
    #             structure = self.dictionary[runnumber]['dotbracket'][counter]
    #             currenttime = time
    #             energy = self.dictionary[runnumber]['energy'][counter]
    #             if counter < totalnumber - 1:
    #                 nexttime = self.dictionary[runnumber]['time'][counter + 1]
    #             else:
    #                 nextime = maxtime + 1
    #             dotdict = _add_structure_timedictionary(dotdict, structure,
    #                             [currenttime, nexttime])
    #             energydict = _add_energy_data(energydict, structure, energy)
    #     outdict  = {}
    #     outdict['dotbracket'] = dotdict
    #     outdict['energy'] = energydict
    #     outdict['sequence'] = self.completesequence
    #     return outdict

    def find_sequence_index(self, sequence):
        """this function will look through the parts and identify the start
        and stop location of the part if it exists in the part
        :param sequence: the sequence of the part that you are looking for
        :type sequence: str
        """
        #First convert the string to RNA
        sequence = dna_to_rna(sequence)
        try:
            lowindex = self.completesequence.index(sequence) + 1
            highindex = lowindex -1 + len(sequence)
        except ValueError:
            return 'Target sequence is not in the complete sequence'
        return [lowindex, highindex]

    def calculate_folding_rate(self, structurelist):
        """This function will fit an exponential function to fit folding data
        with the goal of finding a rate of folding a max value of folding
        sum the contribution of multiple parts for this folding rate
        :param structurelist: list of all of the structures to be considered
        :type structurelist: list of str
        :return: returns popt and pcov
        popt = [shift, , maxvalue, frequency]
        pcov an array of the variance of these values
        """
        #First consolidate all of the structures
        timearray = np.array(self.structuredataframe.index.tolist())
        timearray = timearray - min(timearray)
        structurearray = np.zeros(len(timearray))
        for structure in structurelist:
            structurearray += self.structuredataframe[structure].values
        # timearray = timearray[:1000]
        # structurearray = structurearray[:1000]
        #shift = self.baseadditiontime * self.structurewindow[1]
        popt, pcov = curve_fit(rate_of_folding_func, timearray, structurearray,
            p0=[0.5, 0.0001])
        # popt, pcov = curve_fit(
        #                        lambda t, maxvalue, rate: rate_of_folding_func(
        #                                              t, shift, maxvalue, rate),
                               # temptime, freq)
        return popt, pcov

    # def calculate_folding_rate(self, structurelist):
    #     """This function will fit an exponential function to fit folding data
    #     with the goal of finding a rate of folding a max value of folding
    #     sum the contribution of multiple parts for this folding rate
    #     :param structurelist: list of all of the structures to be considered
    #     :type structurelist: list of str
    #     :return: returns popt and pcov
    #     popt = [shift, , maxvalue, frequency]
    #     pcov an array of the variance of these values
    #     """
    #     #First consolidate all of the structures
    #     for structure in structurelist:
    #         temptime = self.timedata[structure][1:,0]
    #         try:
    #             freq += self.timedata[structure][1:,1]
    #         except NameError:
    #             freq = self.timedata[structure][1:,1]
    #     shift = self.baseadditiontime * self.structurewindow[1]

    #     popt, pcov = curve_fit(
    #                            lambda t, maxvalue, rate: rate_of_folding_func(
    #                                                  t, shift, maxvalue, rate),
    #                            temptime, freq)
    #     return popt, pcov, shift, freq, temptime


###############################################################################
###############################################################################
# CONSOLIDATE RUN DICTIONARY BLOCK START
###############################################################################
###############################################################################
def consolidate_run_dictionary(rundictionary):
    """ This function will distill a run dictionary into a very basic dictionary
    which will be a single set of time points to all of the structures that
    are observed at various times
    :param rundictionary: Dictionary that hyak processing produces
    """
    #First calcualte the pol rate and base timeline
    timebofbaseadition, timeline, sequence = calculate_pol_rate(rundictionary)
    #Rescale all of the time vectors in the rundictionary to align with timelin
    rundictionary = rescale_time_vectors(rundictionary, timeline)
    #Do a final polishing step to make sure that there are not duplicate time points
    #we do this because we want all of the structural shifts to be discrete
    rundictionary = change_duplicate_time_points(rundictionary)
    return rundictionary, timebofbaseadition, sequence

def change_duplicate_time_points(dictionaryofruns):
    """This will scan through all of the runs and add a small
    amount of time to timepoints which are reported as being
    identical"""
    for runnumber in dictionaryofruns:
        timelist = dictionaryofruns[runnumber]['time']
        previoustime = 0
        for index, time in enumerate(timelist):
            currenttime = time
            if previoustime == currenttime:
                #need to shift current time up
                timelist[index] += 0.001
            previoustime = currenttime
        #ictionaryofruns[runnumber]['time'] = timelist
    return dictionaryofruns

def calculate_pol_rate(dictionaryofruns):
    """ function scans through run dictionary to calcuate pol rate if this
    value wasn't previously provided
    """
    timeofbaseaddition = []
        #Polymerization rate
    for runnumber in dictionaryofruns:
        tempbaseadditionlist = []
        dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
        timelist = dictionaryofruns[runnumber]['time']
        if dictionaryofruns[runnumber]['sequence']:
            sequence = dictionaryofruns[runnumber]['sequence']
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
    #Make the timeline
    timeline = []
    for basecount in range(len(sequence)):
        timeline.append(basecount * timeofbaseaddition)
    return timeofbaseaddition, timeline, sequence

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
    totalsize = len(timelist)
    presize = 0
    for count, dotbracket in enumerate(dotbracketlist):
        cursize = len(dotbracket)
        #Second condition accounts for first pass
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
        elif cursize > totalsize:
            timelist[count] += timeline[-1]
    return timelist

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
def compress_run_dictionary(rundictionary, baseadditiontime, completesequence):
    dotdict = OrderedDict()
    energydict = {}
    timelist = calculate_time_list(rundictionary, baseadditiontime)
    #Initialize the dictionary
    for time in timelist:
        dotdict[time] = Counter()
    maxtime = max(timelist)
    for runnumber in rundictionary:
        #print runnumber
        temptimelist = rundictionary[runnumber]['time']
        totalnumber = len(temptimelist)
        #add unstructured sequences structures where appropriate
        mintime = temptimelist[0]
        numberofsinglebasestoadd = int(mintime/baseadditiontime)
        previousbasecountindex = 0
        for basecount in range(numberofsinglebasestoadd):
            temptime = basecount * baseadditiontime
            structure = '.'*(basecount + 1)
            dotdict, previousbasecountindex = add_structure_timedictionary(dotdict, structure,
                        [temptime, temptime + baseadditiontime], previousbasecountindex)
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
    return outdict

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

def calculate_time_list(dictofruns, baseaddition, minstep=25):
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
    timepointstoadd = np.arange(0, mintime, baseaddition)
    outlist = list(timepointstoadd)
    previousstep = 0
    for time in templist:
        # if time > previousstep + minstep:
        #     #print 'this happened'
        #     while time > previousstep + minstep:
        #         previousstep += minstep
        #         outlist.append(previousstep)
        #     outlist.append(time)
        #     previousstep = time
        if time > previousstep:
            #outlist.append(time - 0.005)
            outlist.append(time)
            previousstep = time
    return np.array(outlist)

###############################################################################
###############################################################################
# COMPRESS RUN DICTIONARY BLOCK stop
###############################################################################
###############################################################################
def rate_of_folding_func(t, maxvalue, rate):
    return maxvalue*(1 - np.exp(-rate*(t)))

def caluate_time_vector_for_structures(compresseddict, timewindow,
                   firstexposure = False, windowstartstop = None):
    mintime, maxtime = timewindow
    if windowstartstop:
        lengthofsequence = windowstartstop[1]
    startindex = None
    stopindex = None
    timelist = compresseddict['dotbracket'].keys()
    for index, timepoint in enumerate(timelist):
        if ((not startindex) and timepoint >= mintime):
            if not firstexposure:
                startindex = index
            elif len(compresseddict['dotbracket'][timepoint].most_common(1)[0][0]) == lengthofsequence:
                startindex = index
        elif ((not stopindex) and timepoint >= maxtime):
            stopindex = index
            break
    if not stopindex:
        stopindex = len(timelist) - 1
    return timelist[startindex : stopindex + 1]


def calculate_indexs_from_time(dictionaryofruns, timewindow):
    """ This serves to find the last possible base of the window
    sequence """
    output = [1, 0]
    maxtime = timewindow[1]
    #Cycle through the time points
    for timepoint in dictionaryofruns['dotbracket']:
        maxlength = len(
                dictionaryofruns['dotbracket'][timepoint].most_common[1][0][0])
        if maxtime <= timepoint:
            break
    output[1] = maxlength
    return output

import pandas
def structure_evaluation(dictionaryofruns, structurewindow, timewindow=None,
                         cutoff=0., rescale=False, firstexposure=False):
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
    if timewindow:
        mintime, maxtime = timewindow
    else:
        mintime = 0
        #Grab the max time of the run
        maxtime = dictionaryofruns['dotbracket'].keys()[-1]
    timevector = caluate_time_vector_for_structures(dictionaryofruns, [mintime, maxtime], firstexposure, structurewindow)
    #Now create structure dicitonary and populate it
    structuredataframe = pandas.DataFrame(timevector, columns=['time'])
    structuredataframe = structuredataframe.set_index('time')
    timesize = len(timevector)
    for timepoint in timevector:
        if timepoint >= mintime:
            structuredataframe = add_structure_data(dictionaryofruns['dotbracket'][timepoint], structurewindow, structuredataframe, timepoint, timesize)
        elif timepoint <= maxtime:
            structuredataframe = add_structure_data(dictionaryofruns['dotbracket'][timepoint], structurewindow, structuredataframe, timepoint, timesize)
            break
    return structuredataframe.fillna(0)

def add_structure_data(counterofstructures, structurewindow, structuredataframe, timepoint, sizeoftimevector):
    tempstrudict = Counter()
    for structure in counterofstructures:
        if len(structure) < structurewindow[1]:
            return structuredataframe
        tempstru = structure[structurewindow[0]-1: structurewindow[1]]
        try:
            tempstrudict[tempstru] += counterofstructures[structure]
        except KeyError:
            tempstrudict[tempstru] = counterofstructures[structure]
    for structure in tempstrudict:
        if structure not in structuredataframe:
            structuredataframe[structure] = np.zeros(sizeoftimevector)
        structuredataframe.ix[timepoint, structure] = tempstrudict[structure]
    # print structuredataframe
    # raw_input()
    return structuredataframe

def dna_to_rna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq



# def _add_energy_data(energydict, structure, energy):
#     if structure not in energydict:
#         energydict[structure] = energy
#     return energydict

# def _add_structure_timedictionary(ordereddict, structure, times):
#     for time in ordereddict:
#         if (time >= times[0] and time < times[1]):
#             if structure in ordereddict[time]:
#                 ordereddict[time][structure] += 1
#             else:
#                 ordereddict[time][structure] = 1
#     return ordereddict



# def _calculate_time_list(dictofruns, timewindow=None, rescale=False,
#         minstep=25):
#     templist = [0]
#     for run in dictofruns:
#         templist.extend(dictofruns[run]['time'])
#         #print dictofruns[run]['time']
#     templist = set(templist)
#     templist = list(templist)
#     templist.sort()
#     #print templist
#     if timewindow:
#         maxtime = timewindow[1]
#         mintime = timewindow[0]
#     else:
#         maxtime = max(templist)
#         mintime = min(templist)
#     if rescale and maxtime > max(templist):
#         #asking for data outside of timewindow
#         entry = max(templist)
#         while entry < maxtime:
#             entry += minstep
#             templist.append(entry)
#     outlist = [mintime]
#     previousstep = mintime
#     for time in templist:
#         if time >= maxtime:
#             #outlist.append(maxtime - 0.005)
#             outlist.append(maxtime)
#             break
#         elif time > previousstep:
#             #outlist.append(time - 0.005)
#             outlist.append(time)
#             previousstep = time
#         if time > previousstep + minstep:
#             while time > previousstep + minstep:
#                 previousstep += minstep
#                 outlist.append(previousstep)
#             outlist.append(time)
#             previousstep = time
#     return np.array(outlist)

# def calculate_stats(timecoursedictionary):
#     structs = {}
#     stats = {}
#     for structure in timecoursedictionary:
#         timecountarray = timecoursedictionary[structure]
#         #Find the first instance of it's max
#         maxindex = np.argmax(timecountarray[:, 1])
#         maxfeq = timecountarray[maxindex, 1]
#         maxtime = timecountarray[maxindex, 0]
#         structs[structure] = [maxtime, maxfeq]
#         #finding the first appearance
#         for counter, countofstru in enumerate(timecountarray[:, 1]):
#             if countofstru > 0:
#                 timeoffirst = timecountarray[counter, 0]
#                 break
#         #Finding all of the extrema
#         indexofmaxes = local_maxima_index(timecountarray[:, 1])
#         maxcountvalues = timecountarray[:, 1][indexofmaxes]
#         timeofvalues = timecountarray[:, 0][indexofmaxes]
#         stats[structure] = [timeoffirst, maxcountvalues, timeofvalues]
#         #sorttime dict by time of first appearance
#         structkeys = stats.keys()
#         sortlist = zip([stats[struct][0] for struct in structkeys], structkeys)
#         sortlist.sort()
#         outtimedict = OrderedDict()
#         for time, struct in sortlist:
#             outtimedict[struct] = timecoursedictionary[struct]

#     return structs, stats, outtimedict

# def reduce_size_of_dict(timecoursedictionary):
#     """
#     This function will go through the time and frequency data and remove points
#     Which provide no insight into how the structure changes with time.
#     Basically it will only have time and structure values if there is a dynamic
#     change.
#     """
#     outdict = {}
#     minimumtimestep = 100
#     previousstep = 0
#     for structure in timecoursedictionary:
#         outdict[structure] = []
#         firstentry = False
#         for counter, timefreq in enumerate(timecoursedictionary[structure]):
#             if not firstentry:
#                 firstentry = True
#                 outdict[structure].append(timefreq)
#                 previousstep = timefreq[0]
#                 previousfreq = timefreq[1]
#             #If there is a change in the value we need to add them to the list
#             elif timefreq[1] != previousfreq:
#                 #append the previous value and the current value to appreciate
#                 #the curve
#                 outdict[structure].append(timecoursedictionary[structure][counter - 1])
#                 outdict[structure].append(timefreq)
#                 previousstep = timefreq[0]
#                 previousfreq = timefreq[1]
#             #If we are past the minimum time step we need to add it to the list
#             elif timefreq[0] > previousstep + minimumtimestep:
#                 outdict[structure].append(timefreq)
#                 previousstep = timefreq[0]
#                 previousfreq = timefreq[1]
#     return outdict

# def local_maxima_index(array):
#     gradients = np.diff(array)
#     leftgradeindex = 0
#     listofmaxes = []
#     foundincrease = True
#     for count, grade in enumerate(gradients):
#         if grade == 0:
#             pass
#         if grade < 0 and foundincrease:
#             #Time to restart
#             listofmaxes.append(leftgradeindex)
#             foundincrease = False
#         if grade > 0:
#             leftgradeindex = count + 1
#             foundincrease = True
#     if foundincrease:
#         listofmaxes.append(len(gradients))
#     return listofmaxes


# def add_structures_to_timecountlist(dictionary, structure, timepoints):
#     """
#     This will scan through the timepoints of the dictionary and add the
#     structure and adjust the counts moving forward
#     """
#     sizeofdictionary = len(dictionary[structure])
#     for counter, timecount in enumerate(dictionary[structure]):
#         if timepoints[0] < timecount[0] < timepoints[1]:
#             dictionary[structure][counter][1] += 1
#         elif timepoints[0] == timecount[0]:
#             dictionary[structure][counter][1] += 1
#         elif timecount[0] > timepoints[1]:
#             break
#     return dictionary





###############################################################################
#                               DEPRECIATED
###############################################################################
    # def generate_time_to_structure_dict(self):
    #     outdict = {}
    #     for structure in self.timedata:
    #         for timefreq in self.timedata[structure]:
    #             if timefreq[0] in outdict:
    #                 outdict[timefreq[0]][structure] = timefreq[1]
    #             else:
    #                 outdict[timefreq[0]] = {}
    #                 outdict[timefreq[0]][structure] = timefreq[1]
    #     self.timetostructures = outdict
    #     #return outdict
# def _compensate_for_intial_translation(self, dictionaryofruns):
#     """Kinefold takes timepoint 0 to be the time of the first helix
#     formation. This function accounts for the intial polymerization that
#     has to occur beforehand
#      """
#     #We need to identify the polymerization rate
#     newdict = {}
#     timeofbaseaddition = []
#     #Polymerization rate
#     for runnumber in dictionaryofruns:
#         tempbaseadditionlist = []
#         dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
#         timelist = dictionaryofruns[runnumber]['time']
#         for counter, dotbracket in enumerate(dotbracketlist):
#             currentlength = len(dotbracket)
#             shift = 1
#             while counter-shift >= 0:
#                 if len(dotbracketlist[counter - shift]) == currentlength:
#                     shift += 1
#                 else:
#                     tempbaseadditionlist.append(timelist[counter] - timelist[counter-shift])
#                     break
#             if len(timeofbaseaddition) == 8:
#                 timeofbaseaddition.append(tempbaseadditionlist)
#                 break
#     timeofbaseaddition = np.round(np.mean(tempbaseadditionlist), 2)
#     for runnumber in dictionaryofruns:
#         dotbracketlist = deepcopy(dictionaryofruns[runnumber]['dotbracket'])
#         timelist = deepcopy(dictionaryofruns[runnumber]['time'])
#         if not dotbracketlist:
#             continue
#         elif len(dotbracketlist[0]) == 1:
#             sequence = dictionaryofruns[runnumber]['sequence']
#             continue
#         newdict[runnumber] = {}
#         #Now have to add time and dotbrackets to this system
#         newtimelist = []
#         newdotbracketlist = []
#         newenegylist = []
#         for basenumber in range(len(dotbracketlist[0]) - 1):
#             newtimelist.append(timeofbaseaddition * basenumber)
#             newdotbracketlist.append('.' + '.'*basenumber)
#             newenegylist.append(0)
#         basetime = max(newtimelist) + timeofbaseaddition
#         for counter, dotbracket in enumerate(dotbracketlist):
#             newdotbracketlist.append(dotbracket)
#             newtimelist.append(basetime + timelist[counter])
#         newenegylist.extend(dictionaryofruns[runnumber]['energy'])
#         newdict[runnumber]['dotbracket'] = deepcopy(newdotbracketlist)
#         newdict[runnumber]['time'] = \
#                                 np.round(np.array(deepcopy(newtimelist)), 3)
#         newdict[runnumber]['energy'] = np.array(deepcopy(newenegylist))
#         newdict[runnumber]['sequence'] = dictionaryofruns
#         sequence = dictionaryofruns[runnumber]['sequence']
#     #print newdict
#     if not newdict:
#         newdict = dictionaryofruns
#     return newdict, timeofbaseaddition, sequence

# def structure_evaluation(dictionaryofruns, structurewindow,
#                                     timewindow=None, cutoff=0., rescale=False):
#     """
#     :param dictionaryofruns: This is the standard dictionary that is output by
#         Kinefold runnumber:['dotbracket', etc, etc,]:list
#     :type dictionaryofruns: Dictionary
#     :param structurewindow: The start and stop positions of the dotbrackets to
#         mine over time (1 indexed)
#     :type inddexestomine: list of lists of ints
#         NOTE: This is currently broken and only supports a single index
#     :param timewindow: If used this outlines the timepoints of interest
#     :type timewindow: list of lists of floats
#     :param cutoff: This is a way to remove not common structures - anystructure
#         which doesn't appear more than the cutoff is removed from the output
#     :type cutoff: float
#     """
#     #Probably best to have this set up as a dictionary output?
#     tempdict = {} #keys: structures values: [count, time]
#     maxtime = 0
#     indexs = structurewindow
#     for runnumber in dictionaryofruns:
#         if not dictionaryofruns[runnumber]['time'].any():
#             continue
#         tempdict[runnumber] = {}
#         tempdict[runnumber]['struct'] = []
#         tempdict[runnumber]['time'] = []
#         timelist = dictionaryofruns[runnumber]['time']
#         if maxtime < timelist[-1]:
#             maxtime = timelist[-1]
#         dotbracketlist = dictionaryofruns[runnumber]['dotbracket']
#         #sequencelist = dictionaryofruns[runnumber]['sequence']
#         #sequence = self.sequence[indexstomine[0] - 1: indexstomine[1]]
#         for counter, dotbracket in enumerate(dotbracketlist):
#             if len(dotbracket) >= indexs[1]:
#                 #start recording structure
#                 tempdict[runnumber]['struct'].append(
#                                       dotbracket[indexs[0] - 1: indexs[1]])
#                 tempdict[runnumber]['time'].append(timelist[counter])
#     timelist = _calculate_time_list(tempdict, timewindow, rescale)
#     print timelist
#     maxtime = max(timelist)
#     mintime = min(timelist)
#     #Find the total parts
#     totalparts = len(tempdict)
#     outdict = {}
#     for runnumber in tempdict:
#         # previousstructures = None
#         for counter, structure in enumerate(tempdict[runnumber]['struct']):
#             if counter == len(tempdict[runnumber]['time']) - 1:
#                 timepoints = [tempdict[runnumber]['time'][counter], maxtime]
#             else:
#                 timepoints = tempdict[runnumber]['time'][counter:counter + 2]
#             if timepoints[0] > maxtime:
#                 #print 'breaking'
#                 break
#             elif timepoints[0] == timepoints[1]:
#                 #we don't want to double count these
#                 continue
#             elif timepoints[0] < mintime:
#                 #We don't want to consider earlier timepoints
#                 #Check to see if this structure exists into the actual
#                 #Window
#                 if timepoints[1] > mintime:
#                     timepoints = [mintime, timepoints[1]]
#                 else:
#                     continue
#             if structure in outdict:
#                 outdict = add_structures_to_timecountlist(
#                                                 outdict, structure, timepoints)
#             else:
#                 outdict[structure] = [[ls, 0] for ls in timelist]
#                 # print structure
#                 # print timepoints
#                 outdict = add_structures_to_timecountlist(
#                                                 outdict, structure, timepoints)
#     #Reduce the size of the dictionary
#     #outdict = reduce_size_of_dict(outdict)
#     structstoouput = []
#     for structure in outdict:
#         maxfreq = 0
#         for count, timecount in enumerate(outdict[structure]):
#             tempfreq = timecount[1]/float(totalparts)
#             outdict[structure][count][1] = tempfreq
#             if tempfreq > maxfreq:
#                 maxfreq = tempfreq
#         if maxfreq >= cutoff:
#             structstoouput.append(structure)
#     timecoursedatadict = {x:np.array(y) for (x, y) in
#                             outdict.items() if x in structstoouput}
#     # print timecoursedatadict
#     structures, stats, timecoursedatadict = calculate_stats(timecoursedatadict)
#     return timecoursedatadict, structures, stats
