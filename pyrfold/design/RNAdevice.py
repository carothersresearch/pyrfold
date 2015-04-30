"""
This module should aide in the design of RNAdevices for the design of
RNAparts and devices.
"""
import random

class RNAdevice(object):
    """Basic class for the design of these parts"""
    def __init__(self, partlist, sequencelist):
        self.partlist = partlist
        self.sequencelist = sequencelist
        self.parttoposition = {}
        self.parttosequence = {}
        self.sequence = ''
        #Make dictionary of partname to sequence and partname to sequence
        self.update_part_index_and_sequence_dict()
        #self.update_sequence()

    def __str__(self):
        return convert_to_RNA(self.sequence)

    def change_part_sequence(self, partname, sequence):
        partindex = self.partlist.index(partname)
        self.sequencelist[partindex] = sequence
        self.update_part_index_and_sequence_dict()

    def add_part(self, part, sequence, requestedposition):
        """Will add a part and a sequence to the partlist and sequencelists"""
        self.partlist.insert(requestedposition, part)
        self.sequencelist.insert(requestedposition, sequence)
        self.update_part_index_and_sequence_dict()

    def remove_part(self, part):
        """finds the part that is requested and removes it from the partlist"""
        partindex = self.partlist.index(part)
        self.partlist.pop(partindex)
        self.sequencelist.pop(partindex)
        self.update_part_index_and_sequence_dict()

    def part(self, part):
        return self.parttosequence[part]

    def parts(self, partlist):
        out = ''
        for part in partlist:
            out += self.parttosequence[part]
        return out

    def update_part_index_and_sequence_dict(self):
        listofindexs = []
        leftindex = 0
        for index, part in enumerate(self.partlist):
            rightindex = leftindex + len(self.sequencelist[index])
            listofindexs.append([leftindex, rightindex])
            leftindex = rightindex
        self.parttoposition = dict(zip(self.partlist, listofindexs))
        self.parttosequence = dict(zip(self.partlist, [seq.upper() for seq in self.sequencelist]))
        self.sequence = ''.join(self.sequencelist)


class RNApart(object):
    """Basic class for the design of these parts"""
    def __init__(self, partlist, sequencelist):
        self.partlist = partlist
        self.sequencelist = sequencelist
        self.parttoposition = {}
        self.parttosequence = {}
        self.sequence = ''
        #Make dictionary of partname to sequence and partname to sequence
        self.update_part_index_and_sequence_dict()
        #self.update_sequence()

    def __str__(self):
        return convert_to_RNA(self.sequence)

    def change_part_sequence(self, partname, sequence):
        partindex = self.partlist.index(partname)
        self.sequencelist[partindex] = sequence
        self.update_part_index_and_sequence_dict()

    def add_part(self, part, sequence, requestedposition):
        """Will add a part and a sequence to the partlist and sequencelists"""
        self.partlist.insert(requestedposition, part)
        self.sequencelist.insert(requestedposition, sequence)
        self.update_part_index_and_sequence_dict()

    def remove_part(self, part):
        """finds the part that is requested and removes it from the partlist"""
        partindex = self.partlist.index(part)
        self.partlist.pop(partindex)
        self.sequencelist.pop(partindex)
        self.update_part_index_and_sequence_dict()

    def part(self, part):
        return self.parttosequence[part]

    def parts(self, partlist):
        out = ''
        for part in partlist:
            out += self.parttosequence[part]
        return out

    def update_part_index_and_sequence_dict(self):
        listofindexs = []
        leftindex = 0
        for index, part in enumerate(self.partlist):
            rightindex = leftindex + len(self.sequencelist[index])
            listofindexs.append([leftindex, rightindex])
            leftindex = rightindex
        self.parttoposition = dict(zip(self.partlist, listofindexs))
        self.parttosequence = dict(zip(self.partlist, [seq.upper() for seq in self.sequencelist]))
        self.sequence = ''.join(self.sequencelist)

class RNAsequence(object):
    def __init__(self, sequence):
        self.sequence = convert_to_RNA(sequence)

    def reverse_complement(self):
        tempseq = self.sequence.replace('A', 'x')
        tempseq = tempseq.replace('U', 'A')
        tempseq = tempseq.replace('x', 'U')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        self.sequence = tempseq[::-1]

    def __str__(self):
        return self.sequence

class Helix(object):
    """This class will simply keep track of two sides of helix and
    autoamtically generate one half based on the other half
    """
    def __init__(self, helix0 = '', helix1 = ''):
        """helix0 and helix1 are the left and right half of a helix
        """
        self.helixes = [helix0, helix1]

    def generate_helix(self, basedon = 0, randomly_substitue_u_for_c = False):
        """ this will synthesize the other half of the helix based on
        the squence based on the other half """
        templatehelix = RNAsequence(self.helixes[basedon])
        templatehelix.reverse_complement()
        rev_comp = templatehelix.sequence
        #print rev_comp
        if randomly_substitue_u_for_c:
            sequence = []
            for base in rev_comp:
                if base == 'C':
                    base = random.choice(['C', 'C', 'U'])

                sequence.append(base)
            rev_comp = ''.join(sequence)
        if basedon == 0:
            self.helixes[1] = rev_comp
        else:
            self.helixes[0] = rev_comp

    def randomize_helix(self, sizerange, randomly_substitue_u_for_c=False):
        """This will generate a random helix """
        self.size = random.choice(range(sizerange[0],
                                        sizerange[1] + 1))
        self.helixes[0] = random_sequence(self.size)
        self.generate_helix(self, randomly_substitue_u_for_c=False)

class Unpaired(object):
    """This class will have the basic functionality of making unpaired
    sequences specifically with random generation in mind
    :param sizerange: The lower and upper bound of size of a given sequence
    :type sizerange: A list of two ints.
    :param GC_range: The lower and upper bounds of GC fraction allowed
    :type GC_range: A list of two ints
    """
    def __init__(self, sizerange, sequence='', GC_range=None):
        self.sequence = sequence
        self.size = None
        self.sizerange = sizerange
        self.gcrange = GC_range

    def __str__(self):
        return self.sequence

    def generate_random_sequence(self):
        """
        Randomly creates a sequence based on requirements
        """
        self.size = random.choice(range(self.sizerange[0],
                                        self.sizerange[1] + 1))
        self.sequence = random_sequence(self.size, self.gcrange)

def random_sequence(size, GC_range=None):
    """
    Simple Random RNA generating random sequence generator
    """
    for i in range(10000):
        out = []
        for placecounter in range(size):
            out.append(random.choice(['A', 'U', 'G', 'C']))
        if GC_range:
            if GC_range[0] <= GC_content(out) <= GC_range[1]:
                return ''.join(out)
        else:
            return ''.join(out)

def generate_docking_site(site_size = 20):
    for i in range(10000): #long counter to make sure we find a sequence that works
        out = []
        #First base biased to be a G or C
        first_base = random.choice(['G', 'C', 'G', 'C', 'G', 'C', 'A', 'U'])
        out.append(first_base)
        for count in range(site_size - 1):
            out.append(random.choice(['A', 'U', 'G', 'C']))
        out = ''.join(out)
        #Make sure that the sequence content makes sense
        if 0.35 <= GC_content(out) <= 0.5:
            return ''.join(out)

def GC_content(sequence):
    GC_count = sequence.count('G')
    GC_count += sequence.count('C')
    return float(GC_count)/len(sequence)

def convert_to_RNA(sequence):
    sequence = sequence.upper()
    return sequence.replace('T', 'U')

def convert_to_DNA(sequence):
    sequence = str(sequence).upper()
    return sequence.replace('U', 'T')









