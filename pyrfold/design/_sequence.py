import collections
from ..data import RNAINTERACTION

class Designsequence:
    def __init__(self, sequence, listofdist=None, cutofffreq=0,
            initialize_mers=None):
        """
        :param sequence: RNA sequence of a part
        :type sequence: str
        :param listofdist: list of structures with distributions
        :type listofdist: list
        :param cutofffreq: distributions will only be considered if freq > this
        :type cutofffreq: float
        """
        self.sequence = dna_to_rna(sequence.upper())
        self.exp_sequence = None
        self.exp_mers = None
        self.complexing_exp_mers = None
        self.mers = None
        self.complexing_mers = None
        if listofdist:
            #need to build in this functionality
            self.exp_sequence = extract_exposed_seqeunces(self.sequence,
                                                        listofdist, cutofffreq)
            if initialize_mers:
                self.exp_mers = self.mers_count(initialize_mers,
                    complements=True, expmers=True)
        if initialize_mers:
            self.mers = self.mers_count(initialize_mers, complements=True)

    def mers_count(self, sizes, fworrev='fw', complements=False, expmers=False):
        """ collects the unique mers of a specific size in either forward
        or reverse orentation
        :param size: list of mer sizes to return
        :type size: list of ints
        :param fworrev: orientation of the desired mers
        :type fworrev: str ['fw' or 'rev']
        """
        temp = collections.Counter()
        if expmers:
            for size in sizes:
                temp += mer_count(self.exp_sequence, size)
            if complements:
                self.complexing_exp_mers = generate_complexing_mers(temp)
            if fworrev == 'rev':
                return reverse_complement(temp, 'DNA')
            else:
                return temp
        for size in sizes:
            temp += mer_count([self.sequence], size)
        if complements:
            self.complexing_mers = generate_complexing_mers(temp)
        if fworrev == 'rev':
            return reverse_complement(temp, 'DNA')
        else:
            return temp

    def expmers(self, fworrev, size):
        """ collects all of the unique mers of a given size in either forward
        or reverse orientation that are exposed """
        temp = mer_count(self.exp_sequence, size)
        if fworrev  == 'rev':
            return reverse_complement(temp, 'DNA')
        else:
            return temp

    def complexing_mers_to_list(self, expmers=False):
        output = []
        if expmers:
            for mer in self.complexing_exp_mers:
                output.extend(self.complexing_exp_mers[mer]['complement'])
        else:
            for mer in self.complexing_mers:
                output.extend(self.complexing_mers[mer]['complement'])
        return set(output)



def generate_complexing_mers(mers):
    """
    Takes a Counter which contains mers and in turn generates a dictionary
    which relates the initial seqeunce to all possible interactions
    :mers: some kind of counter

    """
    outdict = {}
    for mer in mers:
        outdict[mer] = {}
        outdict[mer]['complement'] = RNAINTERACTION[len(mer)][mer]
        outdict[mer]['count'] = mers[mer]
    return outdict



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

def mer_count(sequences, size):
    """
    :param sequences: sequences to be analyzed
    :type sequences: list
    :param size: size of the mers to consider
    :type size: int
    """
    output = []
    for sequence in sequences:
        sizeofsequence = len(sequence)
        if sizeofsequence < size:
            continue
        for index in range(size, sizeofsequence + 1):
            output.append(sequence[(index - size): index])
    return collections.Counter(output)

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
            output.append(seq)


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
