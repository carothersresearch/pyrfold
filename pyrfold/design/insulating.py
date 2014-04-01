"""
2014-01-14 12:34 WEV
This module should contain a lot of the core elements needed for analyzing
an insulating section

A good deal of work will have to be done on classes for this system to make
operations smooth and easy to use
"""

########################################################################
# All of the functions that are used to process sequence strings
########################################################################

###
# Modules
###
from random import choice
from .. data import RNAINTERACTION
from copy import copy
from collections import Counter

def recount_mer_interactions(designclassobj):
    """
    This only considers mer in a personal context of mers
    This takes a counter and generates a new counter
    """
    counterofmers = designclassobj.mers
    counterreference = copy(designclassobj.mers)
    for fwmer in counterofmers:
        #Turn find all complement mers
        tempcomps = designclassobj.complexing_mers[fwmer]['complement']
        tempcounter = Counter()
        for compmer in tempcomps:
            if compmer in counterofmers:
                #Add to the the total
                tempcounter[compmer] = counterofmers[compmer]
        #Find the scaling factor
        if sum(tempcounter.values()) == 0:
            continue
        tempscale = float(counterofmers[fwmer])/sum(tempcounter.values())
        for revmer in tempcounter:
            counterreference[revmer] -= tempscale
    for mer in counterreference:
        if counterreference[mer] < 0:
            counterreference[mer] = 0
    return counterreference

def find_scaled_complexing_mers(context, part, exposedparts=False):
    """
    must be designsequence class
    """
    mersofconcern = Counter()
    rescaledcontext = recount_mer_interactions(context)
    #For every partmer find every possible context mer
    for fwmer in rescaledcontext:
        tempcounter = Counter()
        revmers = RNAINTERACTION[len(fwmer)][fwmer]
        for revmer in revmers:
            if revmer in part.mers:
                #The number of partmers that can complex
                tempcounter[revmer] = part.mers[revmer]
        mersofconcern[fwmer] = rescaledcontext[fwmer]*sum(tempcounter.values())
    remlist = []
    for mer in mersofconcern:
        if mersofconcern[mer] == 0.:
            remlist.append(mer)
    for rem in remlist:
        del mersofconcern[rem]
    return mersofconcern

    # for fwmer in part:
    #     tempcounter = Counter()
    #     for revmer in part[fwmer]['complement']:
    #         if revmer in rescaledcontext:
    #             tempcounter[revmer] = rescaledcontext[revmer]




def find_complexing_mers(context, part, exposedparts=False):
    """
    These must be design sequence classes
    """
    mersofconcern = []
    tempcontext = context.complexing_mers
    if exposedparts:
        temppart = part.complexing_exp_mers
    else:
        temppart = part.complexing_mers
    for fwdmercon in tempcontext:
        for fwdmerpart in temppart:
            if fwdmercon in temppart[fwdmerpart]['complement']:
                mersofconcern.extend(tempcontext[fwdmercon]['complement'])
    return set(mersofconcern)

def fraction_of_similarity(denominator, numerator):
    """takes two sets returns the union/denominator
    """
    intersection = set(numerator).intersection(set(denominator))
    return len(intersection)/float(len(denominator))

def sequence_generator(n, nucleictype = 'DNA'):
    """Generates random sequences of nucleic acids
    user defined if DNA or RNA is desired
    nuclictype == 'D' #DNA
    nucleictype == 'R' #RNA
    """
    returnstring = ""
    if nucleictype == 'D' or 'd' or 'DNA':
        possiblebases = ['A', 'C', 'T', 'G']
    if nucleictype == 'R' or 'r' or 'RNA':
        possiblebases = ['A', 'C', 'U', 'G']
    for i in range(n):
        returnstring += choice(possiblebases)
    return returnstring

def dna_to_rna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq

def rna_to_dna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("U","T")
    return seq

def reverse_complement(seq, nucleictype = 'RNA'):
    if nucleictype == 'RNA':
        if isinstance(seq, list) or isinstance(seq, set):
            outlist = []
            for sequence in seq:
                sequence = sequence.upper()[::-1]
                tempstring = ''
                for base in sequence:
                    if base == 'A':
                        tempstring += 'U'
                    elif base == 'U':
                        tempstring += 'A'
                    elif base == 'G':
                        tempstring += 'C'
                    elif base == 'C':
                        tempstring += 'G'
                outlist.append(tempstring)
            return set(outlist)
        if isinstance(seq, str):
            tempstring = ''
            seq = seq.upper()[::-1]
            for base in seq:
                if base == 'A':
                    tempstring += 'U'
                elif base == 'U':
                    tempstring += 'A'
                elif base == 'G':
                    tempstring += 'C'
                elif base == 'C':
                    tempstring += 'G'
            return tempstring
    else:
        if isinstance(seq, list) or isinstance(seq, set):
            outlist = []
            for sequence in seq:
                sequence = sequence.upper()[::-1]
                tempstring = ''
                for base in sequence:
                    if base == 'A':
                        tempstring += 'T'
                    elif base == 'T':
                        tempstring += 'A'
                    elif base == 'G':
                        tempstring += 'C'
                    elif base == 'C':
                        tempstring += 'G'
                outlist.append(tempstring)
            return set(outlist)
        if isinstance(seq, str):
            tempstring = ''
            seq = seq.upper()[::-1]
            for base in seq:
                if base == 'A':
                    tempstring += 'T'
                elif base == 'T':
                    tempstring += 'A'
                elif base == 'G':
                    tempstring += 'C'
                elif base == 'C':
                    tempstring += 'G'
            return tempstring

"""
Functions to process dotbracket data to find exposed sequences
"""
# Core functions
def extract_exposed_seqeunces(sequence, distlist, cutofffreq):
    """(str,[float, str], float) ->[str,str...]
    2013-12-02 18:51 WEV this should index the dotbrackets
    and then pull the sequences from the components
    """
    minimumnumberofbases = 1
    #print sequence
    #print distlist
    listofstartstop = []
    for dist in distlist:
        #dist = [freq, dotbracket]
        if dist[1] > cutofffreq:
            tempdot = dist[1]
            foundstart = 0
            tempstartstop = [0, 0]
            for index, dotorbracket in enumerate(tempdot):
                if (dotorbracket == '.' and not foundstart):
                    #print "FoundSTART"
                    #print index
                    foundstart = 1
                    tempstartstop[0] = index
                elif (dotorbracket != '.' and foundstart):
                    #print "FoundSTOP"
                    #print index
                    foundstart = 0
                    tempstartstop[1] = index - 1
                    #print tempstartstop
                    if (tempstartstop[1] - tempstartstop[0]) > \
                                                        minimumnumberofbases:
                        listofstartstop.append(tempstartstop[:])
            if tempstartstop[0] > tempstartstop[1]:
                tempstartstop[1] = len(sequence)
                if (tempstartstop[1] - tempstartstop[0]) > \
                                                        minimumnumberofbases:
                    listofstartstop.append(tempstartstop[:])
    outputsequences = []
    #print listofstartstop
    for startstop in listofstartstop:
        outputsequences.append(sequence[startstop[0]: startstop[1] + 1])
    #print outputsequences
    outputsequences = set(outputsequences)
    #print outputsequences
    #raw_input("press enter")
    return outputsequences

def unique_mers(sequencelist, sizeofmer):
    """ WEV
    This will create a list of all of the unique Nmers of the object
    The goal is to have this be flexible enough to take a list or string
    """
    output = []
    if isinstance(sequencelist, list) or isinstance(sequencelist, set):
        for sequence in sequencelist:
            sizeofsequence = len(sequence)
            for index in range(sizeofmer, sizeofsequence + 1):
                output.append(sequence[(index - sizeofmer): index])
            #print output
            #raw_input()
        return set(output)
    elif isinstance(sequencelist, str):
        sizeofsequence = len(sequencelist)
        for index in range(sizeofmer, sizeofsequence + 1):
            output.append(sequencelist[(index - sizeofmer): index])
        return set(output)

class PartData:
    def __init__(self, sequence, listofdist, cutofffreq = 0):
        """2013-12-10 11:56 WEV
        sequence is the string of bases
        listofdist is a [float, string] which contains all .()...
        cutoff is the frequency cutoff for exposed parts
        """
        self.sequence = rna_to_dna(sequence.upper())
        self.allthreemers = unique_mers(self.sequence, 3)
        self.exposedsequence = \
            extract_exposed_seqeunces(self.sequence, listofdist, cutofffreq)
        self.exposedthreemers = unique_mers(self.exposedsequence, 3)

    def mers(self, fworrev, size):
        """ collects the unique mers of a specific size in either forward
        or reverse orentation
        """
        temp = unique_mers(self.sequence, size)
        if fworrev == 'rev':
            return reverse_complement(temp, 'DNA')
        else:
            return temp
    def expmers(self, fworrev, size):
        """ collects all of the unique mers of a given size in either forward
        or reverse orientation that are exposed """
        temp = unique_mers(self.exposedsequence, size)
        if fworrev  == 'rev':
            return reverse_complement(temp, 'DNA')
        else:
            return temp

class FlankingData:
    def __init__(self, sequence):
        """2013-12-10 12:04 WEV
        This will process the flanking sequence i.e. 5' and 3'
        to pull out all of unique mers
        """
        self.sequence = rna_to_dna(sequence.upper())
        self.allthreemers = unique_mers(self.sequence, 3)
    def mers(self, fworrev, size):
        temp = unique_mers(self.sequence, size)
        if fworrev == 'rev':
            return reverse_complement(temp, 'DNA')
        else:
            return temp

"""Function to extract all of the similar mers between data sets"""
def similar_mers(desiredcomplementset, revcomplementset):
    """2013-12-04 14:21 WEV
    This should return a list of Nmer sequences which match up
    DesiredCompliment <- should be sequences which are oriented the way
    That the parts are oriented
    RevCompliment <- This is the revComplement - We are looking at binding spot
    """
    # First chance to use set functions!
    return desiredcomplementset.intersection(revcomplementset)

def sequence_similarity_fraction(sequence, listofsequences, tolerance,
    aboveorbelow):
    """WEV this will count the number of sequences from the list which appear
    within the test insulating sequence"""
    totalnumberofsequences = len(listofsequences)
    numberofhits = 0
    for seq in listofsequences:
        if seq in sequence:
            numberofhits += 1
            if float(numberofhits)/totalnumberofsequences >= tolerance:
                if aboveorbelow == 'above':
                    return True
                if aboveorbelow == 'below':
                    return False
    if aboveorbelow == 'below':
        return True
    if aboveorbelow == 'above':
        return False
