"""
This script contains functions which are useful for RNA folding devices
"""
from .vienna import Vienna

def structure_seen_in_mfe(sequence, mfe_structures, temp=37):
    """Will fold a sequence using vienna package and look through the list
    of MFE_structures are given to see if those structures exist

    NOTE: This does not account for specific regions of structure

    :param sequence: Sequence of RNA to be sampled
    :type sequence: str
    :param mfe_sturctures:
    :type mfe_structures:
    :returns: Boolean answer if the structure was seen or not
    :rtype: bool
    """
    temp_vienna = Vienna([sequence])
    structure = temp_vienna.mfe(temp=temp, returnstructure=True)[0][1]
    for mfe_sturcture in mfe_structures:
        if mfe_sturcture in structure:
            return True
    return False


def freq_of_structures_appearing(sequence, list_of_target_structures,
                                 num_structures_to_generate, temp=37):
    """Function will calculate the frequency that structures appear with a
    certain number of folding simulations

    :param sequence: Sequence of RNA to be sampled
    :type sequence: str
    :param list_of_target_structures: All target structures that are considered
    :type list_of_target_structures: list
    :param num_structures_to_generate: The number of randomly sampled
        structures to generate
    :type num_structures_to_generate: int
    :param temp: Temperature to carry out the simulation
    :type temp: float
    :returns: Frequency of the collection of structures appearing
    :rtype: float
    """
    num = num_structures_to_generate
    tofold = Vienna([sequence])
    number_present = 0
    for dotbracket in tofold.pobabalistic_structures(temp, num):
        for target_structure in list_of_target_structures:
            if target_structure in dotbracket:
                number_present += 1
                break
    return float(number_present)/num_structures_to_generate


def dic_of_sampled_structures(sequence, num_structures_to_generate, temp=37):
    """Function uses vienna sampling to build a large dictionary of randomly
    sampled structures. These structures may be used to approximate a
    partition function of structures

    :param sequence: sequence of RNA to be sampled from
    :type sequence: str
    :param num_structures_to_generate: Number of structures to sample.
    :type num_structures_to_generate: int
    :returns: Dictionary of structures with count
    :rtype: Dict
    """
    returndict = {}
    num = num_structures_to_generate
    tofold = Vienna([sequence])
    for dotbracket in tofold.pobabalistic_structures(temp, num):
        try:
            returndict[dotbracket] += 1
        except KeyError:
            returndict[dotbracket] = 1
    return returndict
