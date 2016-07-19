"""
This module should aide in the design of RNAdevices for the design of
RNAparts and devices.

Notes about expanding the functionality of the device. I think it could be
useful to provide a specificpartname list to the device. So that it is possible
to name the specific version of a part that is being used. For example

the aptazyme part could be one of 4 specific rna devices...

===============
Need to make the table of devices that is stored be more descriptive

======
Need to add error conditions for trying to randomize parts that don't
exist within a device

========
Need to write test code for randomizing parts

=========

FIX window part for when the reference part 3 shift is 0

"""
import random
import copy
from ..utilities import convert_to_RNA, random_sequence, \
                        RNA_sequences_complementary, \
                        randomly_substitue_u_for_c, RNAsequence, \
                        reverse_complement
from ..foldingsub import FoldingSubData

# UNIT TESTED CODE


class Device(object):
    """Basic class for the design of these parts"""
    def __init__(self, partnamelist, sequencelist):
        """
        Initializes the object with parnames and corresponding sequences.
        Note that the order of the sequence_list and partnamelist required to be
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

        # should also have Unpaired and helix list
        self.helix_list = []
        self.unpaired_list = []

        # Make dictionary of partname to sequence and partname to sequence
        self.update_part_index_and_sequence_dict()

    def __str__(self):
        return convert_to_RNA(self.sequence)

    def add_helix_part(self, helix, partnames):
        """
        Define the sequence of two parts of the RNA Device with a helix
        object. This will be useful in randomizing sequences of the device.

        :param helix: Helix object with the specified characteristics of
            the parts you want.
        :type helix: Helix object
        :param partnames: This is list of TWO partnames the helix will object
            will be mapped to.
        :type partnames: list of str
        """
        to_add = copy.copy(helix)
        self.helix_list.append((partnames, to_add))
        self.change_part_sequence(partnames[0], to_add.helix0)
        self.change_part_sequence(partnames[1], to_add.helix1)

    def add_unpaired_part(self, unpaired, partname):
        """
        Define the sequence of one part of the RNA Device with an unpaired
        object. This will be useful in randomizing sequences of the device.

        :param helix: Helix object with the specified characteristics of
            the parts you want.
        :type helix: Helix object
        :param partnames: This is list of TWO partnames the helix will object
            will be mapped to.
        :type partnames: list of str
        """
        # to_add = copy.copy(unpaired)
        to_add = unpaired
        self.unpaired_list.append((partname, to_add))
        self.change_part_sequence(partname, to_add)

    def randomize_parts(self, list_of_parts='all'):
        """
        This will randomize all of the RNA Helix and Unpaired parts that
        compose this device
        """
        if list_of_parts == 'all':
            for partnames, helix in self.helix_list:
                helix.randomize()
                self.change_part_sequence(partnames[0], helix.helix0)
                self.change_part_sequence(partnames[1], helix.helix1)
            for partname, unpaired in self.unpaired_list:
                unpaired.randomize()
                # I think this might be unnecessary
                self.change_part_sequence(partname, unpaired)

        else:
            for partnames, helix in self.helix_list:
                for partname in partnames:
                    if partname in list_of_parts:
                        helix.randomize()
                        self.change_part_sequence(partnames[0], helix.helix0)
                        self.change_part_sequence(partnames[1], helix.helix1)

            for partname, unpaired in self.unpaired_list:
                if partname in list_of_parts:
                    unpaired.randomize()
                    # I think this might be unnecessary
                    self.change_part_sequence(partname, unpaired)

        self.update_part_index_and_sequence_dict()

    def change_part_sequence(self, partname, sequence):
        """
        Perminently changes the sequence of a static part

        """
        partindex = self.partnamelist.index(partname)
        self.sequencelist[partindex] = sequence
        self.update_part_index_and_sequence_dict()

    def add_part(self, partname, sequence, partposition):
        """
        Will add a part and a sequence to the partnamelist and sequencelists.
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
        """
        finds the part that is requested and removes it from the partnamelist
        """
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
        self.parttosequence = dict(zip(self.partnamelist, [str(seq).upper() for seq in self.sequencelist]))
        self.sequence = ''
        for part_sequence in self.sequencelist:
            self.sequence += str(part_sequence)

    def return_parts_sequences(self):
        """
        Returns a list of the parts and sequences used to define a sequence
        """
        return ([str(name) for name in self.partnamelist],
                [str(sequence) for sequence in self.sequencelist])

    def return_windowed_sequence(self, fiveprimeshift=None,
                                 fiveprimerefpart=None,
                                 threeprimeshift=None,
                                 threeprimerefpart=None):
        """
        A function which will return a sequence that has been truncated
        based on specificaitons provided.

        Need to have a call for the part sequence thta can be used too
        """

        windowstart_0 = 0
        windowstop_0 = len(self.sequence) + 1

        if fiveprimeshift:
            if fiveprimerefpart:
                # need an index
                start, stop = self.parttoposition[fiveprimerefpart]
            else:
                # if partcontexttofold == 'all':
                #     # Just grab the outside part
                #     start, stop = self.parttoposition[self.partnamelist[0]]
                # else:
                start, stop = self.parttoposition[self.partnamelist[0]]
            # We have the start stop of the device that we'll make the decision
            # on so now we have to shift everything
            # We are going to shift based on the fiveprime side of this
            rel_shift_index = start + fiveprimeshift

            windowstart_0 += rel_shift_index

        if threeprimeshift:
            if threeprimerefpart:
                # need an index
                start, stop = self.parttoposition[threeprimerefpart]
            else:
                # if partcontexttofold == 'all':
                #     # Just grab the outside part
                #     start, stop = self.parttoposition[self.partnamelist[-1]]
                # else:
                start, stop = self.parttoposition[self.partnamelist[-1]]
            # We have the start stop of the device that we'll make the decision
            # on so now we have to shift everything
            # We are going to shift based on the fiveprime side of this
            rel_shift_index = stop + threeprimeshift

            windowstop_0 = rel_shift_index

        return self.sequence[windowstart_0:windowstop_0]

    def create_kinefold_submission_object(self, device_name,
                                          partcontexttofold='all',
                                          partstofold=None,

                                          fiveprimeshift=False,
                                          fiveprimerefpart=None,
                                          threeprimeshift=False,
                                          threeprimerefpart=None,

                                          polrate=30, foldtimeafter=1,
                                          experimenttype=2,
                                          pseudoknots=0, entanglements=0,
                                          numberofsimulations=10,
                                          helix_min_free_engery=6.3460741):
        """
        Creates a kinefold submission object that when added to a list can
        be used to create a submission object for kinefold
        folding. This has all of the functionality of FoldingSubData object -
        that additionally will name simulation objects in a way that will be
        easily deconstructed into a pandas dataframe.

        :param device_name: A base name for the device that is being folded
        :type device_name: str
        :param partcontexttofold: List of all parts housed within the object
            that are to be folded.
        :type partcontexttofold: list of str
        :param partstofold: List of parts that you intend to calcuate folding
            data on.
        :type partstofold: list of str
        :param fiveprimeshift: Shifts the folding window of the seqeunce from
            the five prime end of the sequence. This shifting is always
            relative to the fiveprime portion of a part.
        :type fiveprimeshift: int
        :param fiveprimerefpart: The part that this shifting should be done
            relative to. NOTE shifting is always relative to the five prime
            portion of the part.
        :type fiveprimerefpart: str
        :param threeprimeshift: Shifts the folding window of the seqeunce from
            the three prime end of the sequence. This shifting is always
            relative to the three prim portion of a part.
        :type threeprimeshift: str
        :param threeprimerefpart: The part that this shifting should be done
            relative to. NOTE shifting is always relative to the three prime
            portion of the part.
        :type threeprimerefpart: str
        :param polrate: The rate of polymerase for the reaction [nt/s].
        :type polrate: float
        :param foldtimeafter: Time the sequence is allowed to fold after total
            elongation (for co-trans folding). The total time of folding (melt
            and anneal folding) [s].
        :type foldtimeafter: float
        :param experimenttype: 1 for melt and anneal folding, 2 for co-trans
            folding
        :type experimenttype: int
        :param pseudoknots: 0 for no pseudoknot considerations, 1 for
            pseudoknots to be considered.
        :type pseudoknots: int
        :param entanglements: 0 for entanglements to be ignored. 1 for
            entanglements to be considered.
        :type entanglements: int
        :param numberofsimulations: Number of simulations to be completed for
            a given sequence. Recommended sqrt(length of sequence)
        :type numberofsimulations: int
        :param helix_min_free_engery: The minimum free energy a helix is
            required to contribute for it to be considered [kcal/mol]
        :type helix_min_free_engery: float
        """
        # First create the name of the simulation object
        top_character = '---'
        bottom_character = '#'

        def additional_element_to_name(class_name, value,
                                       top_character=top_character,
                                       bottom_character=bottom_character):
            return top_character + class_name + bottom_character + str(value)

        out_name = 'dev' + bottom_character + device_name

        if partcontexttofold == 'all':
            partcontexttofold = self.partnamelist

        if experimenttype == 2:
            out_name += additional_element_to_name('pol', polrate)
            out_name += additional_element_to_name('fold_time_after',
                                                   foldtimeafter)
        else:
            out_name += additional_element_to_name('anneal_time',
                                                   foldtimeafter)

        if fiveprimeshift:
            out_name += additional_element_to_name('five_prime_shift',
                                                   fiveprimeshift)
            if fiveprimerefpart:
                out_name += additional_element_to_name('rel_five_prime_part',
                                                       fiveprimerefpart)
            else:
                out_name += additional_element_to_name('rel_five_prime_part',
                                                       partcontexttofold[0])
            # calculate five prime shift

        if threeprimeshift:
            out_name += additional_element_to_name('three_prime_shift',
                                                   threeprimeshift)
            if threeprimerefpart:
                out_name += additional_element_to_name('rel_three_prime_part',
                                                       threeprimerefpart)
            else:
                out_name += additional_element_to_name('rel_three_prime_part',
                                                       partcontexttofold[-1])
            # calculate five prime shift
        if partstofold:
            for part in partstofold:
                out_name += additional_element_to_name('part_to_fold',
                                                       part)

        if helix_min_free_engery != 6.3460741:
            out_name += additional_element_to_name('helix_min_energy',
                                                   helix_min_free_engery)

        if pseudoknots:
            out_name += additional_element_to_name('pseudoknots', 'True')
        if entanglements:
            out_name += additional_element_to_name('entanglements', 'True')

        # NOTES FOR determining the indexes to simulate
        # The self.parttoposition is in 0 index and the start-stop windows
        # Are given in 1 indexed

        # I might not actually need to have a relative window stop
        sequence_to_fold = self.combined_part_sequences(partcontexttofold)
        rel_windowstart_0 = self.parttoposition[partcontexttofold[0]][0]
        rel_windowstop_0 = self.parttoposition[partcontexttofold[-1]][1]
        windowstart_0 = 0
        windowstop_0 = self.parttoposition[partcontexttofold[-1]][1] - \
            rel_windowstart_0
        max_window_size = windowstop_0

        if fiveprimeshift:
            if fiveprimerefpart:
                # need an index
                start, stop = self.parttoposition[fiveprimerefpart]
            else:
                # if partcontexttofold == 'all':
                #     # Just grab the outside part
                #     start, stop = self.parttoposition[self.partnamelist[0]]
                # else:
                start, stop = self.parttoposition[partcontexttofold[0]]
            # We have the start stop of the device that we'll make the decision
            # on so now we have to shift everything
            # We are going to shift based on the fiveprime side of this
            rel_shift_index = start + fiveprimeshift

            windowstart_0 += rel_shift_index - rel_windowstart_0

        if threeprimeshift:
            if threeprimerefpart:
                # need an index
                start, stop = self.parttoposition[threeprimerefpart]
            else:
                # if partcontexttofold == 'all':
                #     # Just grab the outside part
                #     start, stop = self.parttoposition[self.partnamelist[-1]]
                # else:
                start, stop = self.parttoposition[partcontexttofold[-1]]
            # We have the start stop of the device that we'll make the decision
            # on so now we have to shift everything
            # We are going to shift based on the fiveprime side of this
            rel_shift_index = stop + threeprimeshift
            windowstop_0 += rel_shift_index - rel_windowstop_0

        if (windowstart_0 < 0) or (windowstart_0 > max_window_size):
            raise IndexError
            print "Shifting window outside of sequence"

        if (windowstop_0 > max_window_size) or (windowstop_0 < 0):
            raise IndexError
            print "Shifting window outside of sequence"

        if windowstart_0 >= windowstop_0:
            raise IndexError
            print "Requested windowstart is greater than or equal to" +\
                " windowstop"

        # Populate the folding sub object
        # Creating the start_stop_list
        part_list = []
        part_start_stop_list = []
        if partstofold:
            for part in partstofold:
                part_list.append(part)
                start, stop = self.parttoposition[part]
                rel_start = start - rel_windowstart_0
                rel_stop = stop - rel_windowstart_0
                # The plus 1 is to shift to 1 indexed
                part_start_stop_list.append([rel_start+1, rel_stop])

        return FoldingSubData(name=out_name, sequence=sequence_to_fold,
                              windowstart=windowstart_0+1,
                              windowstop=windowstop_0,
                              partstartstoplist=part_start_stop_list,
                              partnamelist=part_list,
                              referencepart=None,
                              forcedhelixes=[],
                              polrate=polrate,
                              foldtimeafter=foldtimeafter,
                              experimenttype=experimenttype,
                              pseudoknots=pseudoknots,
                              entanglements=entanglements,
                              numberofsimulations=numberofsimulations,
                              helix_min_free_eng=helix_min_free_engery)


class Helix(object):
    """This class will simply keep track of two sides of helix and
    autoamtically generate one half based on the other half
    """
    def __init__(self, helix0='', helix1='', size_range=None,
                 GC_range=None):
        """
        helix0 and helix1 are the left and right half of a helix.
        :param helix0: RNA sequence of the 5' strand of helix
        :type helix0: str
        :param helix1: RNA sequence of the 3' strand of helix
        :type helix1: str
        """
        self.helix0 = RNAsequence(helix0)
        self.helix1 = RNAsequence(helix1)
        self.sizerange = size_range
        self.gcrange = GC_range

    def generate_helix(self, from_helix_half, random_sub_u_for_c=False):
        """
        This creates a new helix half based on the sequence of a
        helix that was previously entered.
        :param from_helix_half: Either 0 or 1 for the 5' or 3' helix half
            to base the other half on.
        :type from_helix_half: int
        :param random_sub_u_for_c: Should 'U' bases be randomly submitted
            in place of 'C' for the helix compliment. If not False, a
            frequency of substitution should be included.
        :type random_sub_u_for_c: bool or float
        """
        if from_helix_half == 0:
            templatehelix = self.helix0
        elif from_helix_half == 1:
            templatehelix = self.helix1
        else:
            raise IndexError
            print "Requested helix half not 0 or 1"

        rev_comp = reverse_complement(templatehelix)

        if random_sub_u_for_c:
            rev_comp = randomly_substitue_u_for_c(rev_comp,
                                                  random_sub_u_for_c)

        if from_helix_half == 0:
            self.helix1 = rev_comp
        else:
            self.helix0 = rev_comp

    def randomize(self, random_sub_u_for_c=False):
        """
        This will generate a random helix (both the 5' and 3' half).

        :param random_sub_u_for_c: Should 'U' bases be randomly submitted
            in place of 'C' for the helix compliment. If not False, a
            frequency of substitution should be included.
        :type random_sub_u_for_c: bool or float
        """
        self.size = random.choice(range(self.sizerange[0],
                                        self.sizerange[1] + 1))
        self.helix0 = random_sequence(self.size, GC_range=self.gcrange,
                                      strand_type='RNA')
        self.generate_helix(from_helix_half=0,
                            random_sub_u_for_c=random_sub_u_for_c)


class Unpaired(object):
    """This class will have the basic functionality of making unpaired
    sequences specifically with random generation in mind
    """
    def __init__(self, size_range, sequence='', GC_range=None):
        """
        :param sizerange: The lower and upper bound of size of a given sequence
        :type sizerange: A list of two ints.
        :param GC_range: The lower and upper bounds of GC fraction allowed
        :type GC_range: A list of two ints
        """
        self.sequence = sequence
        self.size = None
        self.sizerange = size_range
        self.gcrange = GC_range

    def __repr__(self):
        return 'unpa.' + self.sequence

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def randomize(self):
        """
        Randomly creates a sequence based on requirements assigned
        at object initilization.
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
#     def __init__(self, structure_name_ordered_list, partnamelist):
#         self.structure_name_list = structure_name_ordered_list
#         self.partdictionary = {}
#         for structure_name, structure in zip(structure_name_ordered_list,
#                                              partnamelist):
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
