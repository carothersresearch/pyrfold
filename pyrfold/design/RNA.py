"""
This module should aide in the design of RNAdevices for the design of
RNAparts and devices.
"""
from random import choice
from .utilities import GC_content, reverse_complement, convert_to_RNA, \
    convert_to_DNA, RNA_sequences_complementary, randomly_substitue_u_for_c, \
    random_RNA_sequence


class Part():
    """a class which peices together unpaired, helix, and raw_sequence"""
    def __init__(self, structure_name_ordered_list, partlist):
        self.structure_name_list = structure_name_ordered_list
        self.partdictionary = {}
        for structure_name, structure in zip(structure_name_ordered_list,
                                             partlist):
            # Check if it's a helix
            if 'H' in structure_name:
                if 'p' in structure_name:
                    # We have a postive strand the part will be kept in the
                    # dictionary on the non postive strand
                    continue
                self.partdictionary[structure_name] = structure
            else:
                self.partdictionary[structure_name] = structure

    def __str__(self):
        to_print = []
        for name in self.structure_name_list:
            if 'H' in name:
                if 'p' in name:
                    name = name.split('p')[0]
                    if self.partdictionary[name]:
                        to_print.append(self.partdictionary[name].helixes[1])
                else:
                    if self.partdictionary[name]:
                        to_print.append(self.partdictionary[name].helixes[0])
            else:
                to_print.append(str(self.partdictionary[name]))
        return ''.join(to_print)

    def sequence_of_component(self, component_name):
        if 'H' in component_name:
            if 'p' in component_name:
                component_name = component_name.split('p')[0]
                if self.partdictionary[component_name]:
                    return self.partdictionary[component_name].helixes[1]
            if self.partdictionary[component_name]:
                return self.partdictionary[component_name].helixes[0]
        return str(self.partdictionary[component_name])

    def reset_sequences(self, list_of_parts):
        for part in list_of_parts:
            self.partdictionary[part] = ''


class Helix():
    """This class will simply keep track of two sides of helix and
    autoamtically generate one half based on the other half
    """
    def __init__(self, helix0='', helix1='', GC_content=None, sizerange=None):
        """helix0 and helix1 are the left and right half of a helix
        """
        # First Check Complementarity
        if not RNA_sequences_complementary(helix0, helix1):
            # Write real error function
            print 'sequences not complentary'

        self.helixes = [helix0, helix1]
        self.GC_content = GC_content
        self.sizerange = sizerange


class Unpaired():
    """This class will have the basic functionality of making unpaired
    sequences specifically with random generation in mind
    :param sizerange: The lower and upper bound of size of a given sequence
    :type sizerange: A list of two ints.
    :param GC_range: The lower and upper bounds of GC fraction allowed
    :type GC_range: float
    """
    def __init__(self, sequence, sizerange=None, GC_range=None):
        self.sequence = sequence
        self.GC_range = GC_range
        self.sizerange = sizerange

    def __str__(self):
        return self.sequence

    def random_sequence(self):
        """
        Randomly creates a unpaired sequence object based on requirements
        """
        # Need to write error cases for if no sizerange is included
        if self.sizerange[0] == self.sizerange[0]:
            size = self.sizerange[0]
        else:
            size = choice(range(self.sizerange[0], self.sizerange[1] + 1))
        self.sequence = random_RNA_sequence(size, self.GC_range)


def generate_random_helix(sizerange, randomly_sub_u_for_c=False,
                          GC_range=None):
    """This will generate a random helix object
    :return type: object
    """
    size = choice(range(sizerange[0], sizerange[1] + 1))

    neg_helix = random_RNA_sequence(size, GC_range)
    pos_helix = reverse_complement(neg_helix)

    if randomly_substitue_u_for_c:
        pos_helix = randomly_substitue_u_for_c(pos_helix)
    # This last random choice is to make sure there is no biase with G_C sub
    if choice([True, False]):
        return Helix(neg_helix, pos_helix, GC_range, sizerange)
    else:
        return Helix(pos_helix, neg_helix, GC_range, sizerange)


def generate_half_helix(sequence, strand=0, randomly_sub_u_for_c=False):
    """ this will synthesize the other half of the helix based on
    the squence based on the other half """
    rev_comp = reverse_complement(sequence)
    # print rev_comp
    if randomly_sub_u_for_c:
        rev_comp = randomly_substitue_u_for_c(rev_comp)
    if strand == 0:
        return Helix(sequence, rev_comp)
    elif strand == 1:
        return Helix(rev_comp, sequence)





def generate_docking_site(site_size=20):
    for i in range(10000):  # long counter tfind a sequence that works
        out = []
        # First base biased to be a G or C
        first_base = choice(['G', 'C', 'G', 'C', 'G', 'C', 'A', 'U'])
        out.append(first_base)
        for count in range(site_size - 1):
            out.append(choice(['A', 'U', 'G', 'C']))
        out = ''.join(out)
        # Make sure that the sequence content makes sense
        if 0.35 <= GC_content(out) <= 0.5:
            return ''.join(out)
