"""collection of data and utility functions for general RNA and
DNA work
"""
from random import random, choice
import pandas as pd
import os
import warnings

# Dictionaries of complement sequences
RNA_complement_dict_gu = {'A': 'U',
                          'U': ['A', 'G'],
                          'G': ['C', 'U'],
                          'C': 'G'}

RNA_complement_dict = {'A': 'U',
                       'U': 'A',
                       'G': 'C',
                       'C': 'G'}


DNA_complement_dict = {'A': 'T',
                       'T': 'A',
                       'G': 'C',
                       'C': 'G'}

# UNITTESTED CODE


def index_of_sequence(sequence, sub_sequence):
    """this function will look through the parts and identify the start
    and stop location of the part if it exists in the part.
    NOTE: This is 1 indeded
    :param sequence: The full sequence to search
    :type sequence: str
    :param sub_sequence: The sub sequence to search for
    :type sub_sequence: str
    :returns: Start and stop index (1 indexed).
    :rtype: list
    """
    # First convert the string to RNA
    sequence = convert_to_RNA(sequence)
    sub_sequence = convert_to_RNA(sub_sequence)
    try:
        lowindex = sequence.index(sub_sequence) + 1
        highindex = lowindex - 1 + len(sub_sequence)
    except ValueError:
        return 'Target sequence is not in the complete sequence'
    return [lowindex, highindex]


def random_sequence(size, GC_range=None, strand_type='RNA'):
    """
    Simple random nucleotide sequence generator
    Nucleobases accepts RNA or DNA as input

    :type size: int
    :param size: Length of sequence to be returned
    :type GC_range: list of floats
    :param GC_range: Low and high bounds of GC content to be accepted
    :type nucleobase: str
    :param nucleobase: RNA or DNA as input for the type of sequence to return
    """
    if strand_type == 'RNA':
        base_set = ['A', 'U', 'G', 'C']
    elif strand_type == 'DNA':
        base_set = ['A', 'T', 'G', 'C']
    else:
        raise ValueError('Strand_type accepts RNA or DNA')

    # Check GC bounds
    if GC_range:
        if (GC_range[0] == 0) or (GC_range[1] == 1):
            pass
        elif (GC_range[0]*size > (size - 1) and (GC_range[1]*size < 1)):
            warnings.warn("Impossible GC_range requested, removed constraint")
            GC_range = None
        elif GC_range[0] > GC_range[1]:
            raise ValueError('GR_range must go from low to high')

        # Randomly decide how many Gs to add
        ## pick random GC content
        # lower_bound = int(math.ceil(GC_range[0]*size))
        # upper_bound = int(GC_range[1]*size)
        # number_of_gs = np.random.randint(lower_bound, upper_bound + 1)

    for i in xrange(100000):

        out = []
        for placecounter in range(size):
            out.append(choice(base_set))
        if GC_range:
            if ((GC_content(''.join(out)) <= GC_range[1]) &
               (GC_content(''.join(out)) >= GC_range[0])):
                return ''.join(out)
        else:
            return ''.join(out)


def GC_content(sequence):
    sequence = str(sequence)
    GC_count = sequence.count('G')
    GC_count += sequence.count('C')
    return float(GC_count)/len(sequence)


def reverse_complement(sequence, strand_type='RNA'):
    """ This block simply creates the reverse complement of a strand of either
    DNA or RNA
    :param sequence: The sequence of DNA/RNA to be reversed
    :type sequence: str
    :param strand_type: The specific strand type that is being used
    :type strand_type: str
    """
    if strand_type == 'RNA':
        sequence = convert_to_RNA(sequence)
        tempseq = sequence.replace('A', 'x')
        tempseq = tempseq.replace('U', 'A')
        tempseq = tempseq.replace('x', 'U')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        sequence = tempseq[::-1]
    if strand_type == 'DNA':
        sequence = convert_to_DNA(sequence)
        tempseq = sequence.replace('A', 'x')
        tempseq = tempseq.replace('T', 'A')
        tempseq = tempseq.replace('x', 'T')
        tempseq = tempseq.replace('G', 'x')
        tempseq = tempseq.replace('C', 'G')
        tempseq = tempseq.replace('x', 'C')
        sequence = tempseq[::-1]
    return sequence


def RNA_sequences_complementary(sequence1, sequence2, gu_wobble=True):
    """
    Function to determine if two RNA sequences are complementary

    # TODO

    """
    # Make sure the sequences are RNA
    if gu_wobble:
        complement_dict = RNA_complement_dict_gu
    else:
        complement_dict = RNA_complement_dict
    sequence1 = convert_to_RNA(str(sequence1))
    sequence2 = convert_to_RNA(str(sequence2))
    if len(sequence1) != len(sequence2):
        # write error script
        raise ValueError
        print "Supplied sequences are different lengths"
    for position, base1 in enumerate(sequence1):
        if sequence2[-(position+1)] in complement_dict[base1]:
            continue
        else:
            return False
    return True


# NOT UNITTESTED CODE


def create_dataframe_of_submission(list_of_submission_objs):
    """
    This function will create a pandas dataframe which contains all of the
    sequence information for a device.
    :param list_of_submission_objs: A list containing all of the submission
        objects that are to be recroded.
    :type list_of_submission_objs: list
    """
    sequence = []
    index = []
    for sub_obj in list_of_submission_objs:
        sequence.append(sub_obj.sequence)
        dict_of_name = file_name_splitter(sub_obj.name)
        index.append(dict_of_name['dev'][0])

    return pd.Series(sequence, index)


def create_dataframe_of_part_information(part_information_list, name_list):
    """
    Returns a pandas datafame which has all of the detailed part information is
    stored in columns.
    :param part_information_list: List of part information. Specifically
        containing orderd lists of the parts and corresponding sequences.
    :type part_information_list: list of tuples
    :param name_list: Ordered list of part names corresponding to the
        part_information_list data.
    :type name_list: list of str
    :returns: Dataframe of all part information
    :rtype: Pandas.DataFrame
    """
    output_dataframe = pd.DataFrame()
    for (partnames, sequences), name in zip(part_information_list, name_list):
        output_dataframe = output_dataframe.append(pd.DataFrame(
                                                   data=[sequences],
                                                   index=[name],
                                                   columns=partnames))
    return output_dataframe


def file_name_splitter(input_file_name, top_character='---',
                       bottom_character='#'):
    """
    This function will split the filename into it's part and build a
    a dictionary of its function. The filename is expecting the
    parameters to be be divided by '#' and the value to be devided
    by '&'
    :input_file_name type: string
    :input_file_name param: Filename to be divided
    """
    out_dictionary = {}
    # First removing the filepath
    file_name = os.path.splitext(os.path.basename(input_file_name))[0]
    variables = file_name.split(top_character)
    for variable in variables:
        variable_name, value = variable.split(bottom_character)
        try:
            out_dictionary[variable_name].append(value)
        except:
            out_dictionary[variable_name] = [value]
    return out_dictionary


class RNAsequence(object):
    def __init__(self, sequence):
        self.sequence = convert_to_RNA(sequence)

    def reverse_complement(self):
        self.sequence = reverse_complement(self.sequence, strand_type='RNA')

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)


def randomly_substitue_u_for_c(sequence, probability=0.5):
    """
    Function for randomly_substituing Us for Cs for a given sequence
    """
    sequence = str(sequence)
    outsequence = []
    for base in sequence:
        if base == 'C':
            if random() <= probability:
                base = 'U'
        outsequence.append(base)
    return ''.join(outsequence)


def convert_to_RNA(sequence):
    """
    Concerts DNA sequence to RNA Sequence
    """
    sequence = str(sequence)
    sequence = sequence.upper()
    return sequence.replace('T', 'U')


def convert_to_DNA(sequence):
    """
    Converts RNA to DNA
    """
    sequence = str(sequence)
    sequence = sequence.upper()
    return sequence.replace('U', 'T')


def random_RNA_sequence(size, GC_range=None):
    """
    Simple Random RNA generating random sequence generator
    """
    return random_sequence(size, GC_range, strand_type='RNA')


def random_DNA_sequence(size, GC_range=None):
    """
    Simple Random DNA generating random sequence generator
    """
    return random_sequence(size, GC_range, strand_type='DNA')
