"""
This module exists to do basic processing of timecourse data that is output
from kinefld simulations


Possible things TODO:
add a wildcardstructure thing like ****** when structures not yet exposed

"""

import math
import numpy as np
import pandas
# from scipy.signal import argrelextrema
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
        from ..hyak import process as hyakp

        temprundict, baseadditiontime, completesequence, simulation_type = \
                            hyakp.consolidate_run_dictionary(dictionaryofruns)
        #This step is time intensive and needs to be optimized
        dictionary = hyakp.compress_run_dictionary(temprundict,
                                  baseadditiontime, completesequence,
                                  simulation_type)
        #A list of structures that have been found thus far and their max
        #frequency
        output = cls(dictionary, structurewindow, timewindow, rescale,
                                        firstexposure, maxlength, cutoff)
        return output


    def generate_data(self, structurewindow, timewindow=None, cutoff=0.0,
                      maxlength=None, firstexposure=True):
        """ This function will condense the run dictionary information into a
        pandas dataframe that is windowed by the specified bounds.
        :param structurewindow: The start and stop location of structure (1
            indexed)
        :type structurewindow: (int, int)
        :param timewindow: The timewindow of the simulation to consider when
            generating the data
        :type timewindow: (float, float)
        :param cutoff: The minimum frequency a structure must have to be
            considered
        :type cutoff: float
        :param maxlength:
        :type maxlength:
        :param firstexposure: Should the dataframe start after the first base
            has been exposed?
        :type firstexposure: bool
        :returns: Populates the structure dataframe of the object
        :rtype: None
        """
        self.structurewindow = structurewindow
        self.timewindow = timewindow
        if maxlength:
            self.structurewindow = calculate_indexs_from_time(
                                        self.dictionary, self.timewindow)
        self.structuredataframe = \
                structure_evaluation(self.dictionary, self.structurewindow,
                               self.timewindow, cutoff, firstexposure)
        self.sequence = self.completesequence[
                      self.structurewindow[0] - 1: self.structurewindow[1]]


    def find_sequence_index(self, sequence):
        """this function will look through the parts and identify the start
        and stop location of the part if it exists in the part
        :param sequence: the sequence of the part that you are looking for
        :type sequence: str
        """
        # First convert the string to RNA
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
        popt = [maxvalue, rate]
        pcov an array of the variance of these values
        """
        # First consolidate all of the structures
        timearray = np.array(self.structuredataframe.index.tolist())
        timearray = timearray - min(timearray)
        structurearray = np.zeros(len(timearray))
        for structure in structurelist:
            structurearray += self.structuredataframe[structure].values
        # timearray = timearray[:1000]
        # structurearray = structurearray[:1000]
        # shift = self.baseadditiontime * self.structurewindow[1]
        popt, pcov = curve_fit(rate_of_folding_func, timearray,
                               structurearray, p0=[0.5, 0.0001])
        # popt, pcov = curve_fit(
        #                        lambda t, maxvalue, rate: rate_of_folding_func(
        #                                              t, shift, maxvalue, rate),
        # temptime, freq)
        return popt, pcov

    def maxium_frequency_of_structures(self, structurelist):
        """
        This function searches the populated structuredataframe for all of the
        structures that are requested. Right now this requires that you have
        generated a dataframe with the appropriate window
        :param structurelist: list of all of the structures to be considered
        :type structurelist: list of str
        :return: Maximum frequency of a combined set of structures appearing
        :rtype: float
        """
        structurearray = np.zeros(self.structuredataframe.shape[0])
        for structure in structurelist:
            try:
                structurearray += self.structuredataframe[structure].values
            except KeyError:
                pass
        return structurearray.max()

    def final_frequency_of_structures(self, structurelist):
        """
        This function searches the populated structuredataframe for all of the
        structures that are requested and then returns the final frequency
        of the collection of structures
        :param structurelist: list of all of the structures to be considered
        :type structurelist: list of str
        :return: Maximum frequency of a combined set of structures appearing
        :rtype: float
        """
        return_frequency = 0.
        for structure, frequency in self.structuredataframe.iloc[-1].iteritems():
            if structure in structurelist:
                return_frequency += frequency
        return return_frequency

    def average_frequency_of_structures(self, structurelist):
        """
        This function will process an exisiting structuredataframe. It searches
        for every strcuture that is asked for, creates a vector of their total
        folding frequency, and then integrates under this data and normalizes
        to the toal time that is elapsed
        :param structurelist: list of all of the structures to be considered
        :type structurelist: list of str
        :return: Average frequency structures appeared
        :rtype: float
        """
        structurearray = np.zeros(self.structuredataframe.shape[0])
        for structure in structurelist:
            try:
                structurearray += self.structuredataframe[structure].values
            except KeyError:
                pass
        total_time = self.structuredataframe.index[-1] - \
                     self.structuredataframe.index[0]
        return np.trapz(x=self.structuredataframe.index, y=structurearray) / \
               total_time

    def final_structures_seen(self, structurewindow, cutoff=0.0):
        """Generates data for the specific window and then returns the windowed
        structures for the devices
        :param structurewindow:
        :type structurewindow:
        :param cutoff: The maximum frequency of structure to be considered
        :type cutoff: float
        :returns: List of final dotbracket structures seen
        :rtype: list
        """
        return_list = []
        self.generate_data(structurewindow,  cutoff=cutoff,
                           firstexposure=True)
        for structure, frequency in self.structuredataframe.iloc[-1].iteritems():
            if frequency > cutoff:
                return_list.append(structure)
        return return_list

    def __repr__(self):
        # Crude Sorting system
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

def rate_of_folding_func(t, maxvalue, rate):
    return maxvalue*(1 - np.exp(-rate*(t)))

def calcuate_time_vector_for_structures(compresseddict, timewindow,
                                        firstexposure=False,
                                        windowstartstop=None):
    mintime, maxtime = timewindow
    if windowstartstop:
        lengthofsequence = windowstartstop[1]
    startindex = None
    stopindex = None
    timelist = compresseddict['dotbracket'].keys()
    for index, timepoint in enumerate(timelist):
        if firstexposure and startindex is None:
            size = len(compresseddict['dotbracket'][timepoint].most_common(1)[0][0])
            if size == lengthofsequence:
                startindex = index
        elif startindex is None and timepoint >= mintime:
            startindex = index
        elif stopindex is None and timepoint >= maxtime:
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


def structure_evaluation(dictionaryofruns, structurewindow, timewindow=None,
                         cutoff=0., firstexposure=False):
    """
    :param dictionaryofruns: This is the standard dictionary that is output by
        Kinefold runnumber:['dotbracket', etc, etc,]:list
    :type dictionaryofruns: Dictionary
    :param structurewindow: The start and stop positions of the dotbrackets to
        mine over time (1 indexed)
    :param timewindow: If used this outlines the timepoints of interest
    :type timewindow: list of lists of floats
    :param cutoff: This is a way to remove not common structures - anystructure
        which doesn't appear more than the cutoff is removed from the output
    :type cutoff: float
    """
    # Probably best to have this set up as a dictionary output?
    if timewindow:
        mintime, maxtime = timewindow
    else:
        mintime = 0
        # Grab the max time of the run
        maxtime = dictionaryofruns['dotbracket'].keys()[-1]
    timevector = calcuate_time_vector_for_structures(dictionaryofruns, [mintime, maxtime], firstexposure, structurewindow)
    # Now create structure dicitonary and populate it
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

def add_structure_data(counterofstructures, structurewindow,
                       structuredataframe, timepoint, sizeoftimevector):
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
