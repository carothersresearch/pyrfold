"""
This module should aide in the design of RNAdevices for the design of
RNAparts and devices.
"""
import random
import copy
from ..utilities import convert_to_RNA, random_sequence, \
                        RNA_sequences_complementary, \
                        randomly_substitue_u_for_c, RNAsequence

# UNIT TESTED CODE

class Device(object):
    """Basic class for the design of these parts"""
    def __init__(self, partnamelist, sequencelist):
        """
        Initializes the object with parnames and corresponding sequences.
        Note that the order of the sequence_list and partlist required to be
        in the same order.

        :param partnamelist: List of partnames
        :type partnamelist: list of str
        :param sequencelist: List of sequences corresponding to names
        :type sequencelist:  lst of str
        """
        self.partnamelist = copy.copy(partnamelist)
        self.sequencelist = copy.copy(sequencelist)
        self.parttoposition = {}
        self.parttosequence = {}
        self.sequence = ''
        # Make dictionary of partname to sequence and partname to sequence
        self.update_part_index_and_sequence_dict()


    def __str__(self):
        return convert_to_RNA(self.sequence)

    def change_part_sequence(self, partname, sequence):
        partindex = self.partnamelist.index(partname)
        self.sequencelist[partindex] = sequence
        self.update_part_index_and_sequence_dict()

    def add_part(self, partname, sequence, partposition):
        """
        Will add a part and a sequence to the partlist and sequencelists.
        :param partname: Name of the part that is being inserted
        :type partname: str
        :param sequence: Sequence of part that is being inserted
        :type sequence: str
        :param partposition: Index of the part to be added (0 indexed)
        :type partposition: int
        """
        self.partnamelist.insert(partposition, partname)
        self.sequencelist.insert(partposition, sequence)
        self.update_part_index_and_sequence_dict()

    def remove_part(self, partname):
        """finds the part that is requested and removes it from the partlist"""
        partindex = self.partnamelist.index(partname)
        self.partnamelist.pop(partindex)
        self.sequencelist.pop(partindex)
        self.update_part_index_and_sequence_dict()

    def part_sequence(self, part):
        return self.parttosequence[part]

    def combined_part_sequences(self, partnamelist):
        out = ''
        for part in partnamelist:
            out += self.parttosequence[part]
        return out

    def update_part_index_and_sequence_dict(self):
        listofindexs = []
        leftindex = 0
        for index, part in enumerate(self.partnamelist):
            rightindex = leftindex + len(self.sequencelist[index])
            listofindexs.append([leftindex, rightindex])
            leftindex = rightindex
        self.parttoposition = dict(zip(self.partnamelist, listofindexs))
        self.parttosequence = dict(zip(self.partnamelist, [seq.upper() for seq in self.sequencelist]))
        self.sequence = ''.join(self.sequencelist)


class Helix(object):
    """This class will simply keep track of two sides of helix and
    autoamtically generate one half based on the other half
    """
    def __init__(self, helix0='', helix1=''):
        """helix0 and helix1 are the left and right half of a helix
        """
        self.helixes = [helix0, helix1]

    def generate_helix(self, based_on_half, random_sub_u_for_c=False):
        """ this will synthesize the other half of the helix based on
        the squence based on the other half """
        templatehelix = RNAsequence(self.helixes[based_on_half])
        templatehelix.reverse_complement()
        rev_comp = templatehelix.sequence

        if random_sub_u_for_c:
            rev_comp = randomly_substitue_u_for_c(rev_comp, probability=0.5)
        if based_on_half == 0:
            self.helixes[1] = rev_comp
        else:
            self.helixes[0] = rev_comp

    def randomize_helix(self, sizerange, random_sub_u_for_c=False):
        """This will generate a random helix """
        self.size = random.choice(range(sizerange[0],
                                        sizerange[1] + 1))
        self.helixes[0] = random_sequence(self.size)
        self.generate_helix(self, random_sub_u_for_c=False)


class Unpaired(object):
    """This class will have the basic functionality of making unpaired
    sequences specifically with random generation in mind
    :param sizerange: The lower and upper bound of size of a given sequence
    :type sizerange: A list of two ints.
    :param GC_range: The lower and upper bounds of GC fraction allowed
    :type GC_range: A list of two ints
    """
    def __init__(self, size_range, sequence='', GC_range=None):
        self.sequence = sequence
        self.size = None
        self.sizerange = size_range
        self.gcrange = GC_range

    def __str__(self):
        return self.sequence

    def generate_random_sequence(self):
        """
        Randomly creates a sequence based on requirements
        """
        self.size = random.choice(range(self.sizerange[0],
                                        self.sizerange[1] + 1))
        self.sequence = random_sequence(self.size, self.gcrange,
                                        strand_type='RNA')

# NOT UNIT TESTED CODE


# OLD code that could be repurposed for future use

# """
# This module should aide in the design of RNAdevices for the design of
# RNAparts and devices.
# """
# from random import choice
# from ..utilities import GC_content, reverse_complement, convert_to_RNA, \
#     convert_to_DNA, RNA_sequences_complementary, randomly_substitue_u_for_c, \
#     random_RNA_sequence


# class Part():
#     """a class which peices together unpaired, helix, and raw_sequence"""
#     def __init__(self, structure_name_ordered_list, partlist):
#         self.structure_name_list = structure_name_ordered_list
#         self.partdictionary = {}
#         for structure_name, structure in zip(structure_name_ordered_list,
#                                              partlist):
#             # Check if it's a helix
#             if 'H' in structure_name:
#                 if 'p' in structure_name:
#                     # We have a postive strand the part will be kept in the
#                     # dictionary on the non postive strand
#                     continue
#                 self.partdictionary[structure_name] = structure
#             else:
#                 self.partdictionary[structure_name] = structure

#     def __str__(self):
#         to_print = []
#         for name in self.structure_name_list:
#             if 'H' in name:
#                 if 'p' in name:
#                     name = name.split('p')[0]
#                     if self.partdictionary[name]:
#                         to_print.append(self.partdictionary[name].helixes[1])
#                 else:
#                     if self.partdictionary[name]:
#                         to_print.append(self.partdictionary[name].helixes[0])
#             else:
#                 to_print.append(str(self.partdictionary[name]))
#         return ''.join(to_print)

#     def sequence_of_component(self, component_name):
#         if 'H' in component_name:
#             if 'p' in component_name:
#                 component_name = component_name.split('p')[0]
#                 if self.partdictionary[component_name]:
#                     return self.partdictionary[component_name].helixes[1]
#             if self.partdictionary[component_name]:
#                 return self.partdictionary[component_name].helixes[0]
#         return str(self.partdictionary[component_name])

#     def reset_sequences(self, list_of_parts):
#         for part in list_of_parts:
#             self.partdictionary[part] = ''


# class Helix():
#     """This class will simply keep track of two sides of helix and
#     autoamtically generate one half based on the other half
#     """
#     def __init__(self, helix0='', helix1='', GC_content=None, sizerange=None):
#         """helix0 and helix1 are the left and right half of a helix
#         """
#         # First Check Complementarity
#         if not RNA_sequences_complementary(helix0, helix1):
#             # Write real error function
#             print 'sequences not complentary'

#         self.helixes = [helix0, helix1]
#         self.GC_content = GC_content
#         self.sizerange = sizerange


# class Unpaired():
#     """This class will have the basic functionality of making unpaired
#     sequences specifically with random generation in mind
#     :param sizerange: The lower and upper bound of size of a given sequence
#     :type sizerange: A list of two ints.
#     :param GC_range: The lower and upper bounds of GC fraction allowed
#     :type GC_range: float
#     """
#     def __init__(self, sequence, sizerange=None, GC_range=None):
#         self.sequence = sequence
#         self.GC_range = GC_range
#         self.sizerange = sizerange

#     def __str__(self):
#         return self.sequence

#     def random_sequence(self):
#         """
#         Randomly creates a unpaired sequence object based on requirements
#         """
#         # Need to write error cases for if no sizerange is included
#         if self.sizerange[0] == self.sizerange[0]:
#             size = self.sizerange[0]
#         else:
#             size = choice(range(self.sizerange[0], self.sizerange[1] + 1))
#         self.sequence = random_RNA_sequence(size, self.GC_range)


# def generate_random_helix(sizerange, randomly_sub_u_for_c=False,
#                           GC_range=None):
#     """This will generate a random helix object
#     :return type: object
#     """
#     size = choice(range(sizerange[0], sizerange[1] + 1))

#     neg_helix = random_RNA_sequence(size, GC_range)
#     pos_helix = reverse_complement(neg_helix)

#     if randomly_substitue_u_for_c:
#         pos_helix = randomly_substitue_u_for_c(pos_helix)
#     # This last random choice is to make sure there is no biase with G_C sub
#     if choice([True, False]):
#         return Helix(neg_helix, pos_helix, GC_range, sizerange)
#     else:
#         return Helix(pos_helix, neg_helix, GC_range, sizerange)


# def generate_half_helix(sequence, strand=0, randomly_sub_u_for_c=False):
#     """ this will synthesize the other half of the helix based on
#     the squence based on the other half """
#     rev_comp = reverse_complement(sequence)
#     # print rev_comp
#     if randomly_sub_u_for_c:
#         rev_comp = randomly_substitue_u_for_c(rev_comp)
#     if strand == 0:
#         return Helix(sequence, rev_comp)
#     elif strand == 1:
#         return Helix(rev_comp, sequence)


