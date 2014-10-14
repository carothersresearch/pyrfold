"""
This module should aide in the design of RNAdevices for the design of
RNAparts and devices.
"""

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
        return self.sequence

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


class RNAsequence():
    def __init__(self, sequence):
        self.sequence = self.convert_to_RNA(sequence)

    def reverse_complement(self):
        tempseq = self.sequence.replace('A', 'x')
        tempseq = tempseq.replace('U', 'A')
        tempseq = tempseq.replace('x', 'U')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        self.sequence = tempseq[::-1]

    def convert_to_RNA(self, sequence):
        sequence = sequence.upper()
        return sequence.replace('T', 'U')

    def __str__(self):
        return self.sequence


import random

class helix():
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
