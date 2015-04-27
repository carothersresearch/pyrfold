"""
This script contains functions which are useful for RNA folding devices
"""
from .vienna import Vienna

def structure_seen_in_mfe(sequence, mfe_structures, region=None, temp=37):
    """Will fold a sequence using vienna package and look through the list
    of MFE_structures are given to see if those structures exist

    NOTE: This does not account for specific regions of structure

    :param sequence: Sequence of RNA to be sampled
    :type sequence: str
    :param mfe_sturctures: List of structures to look for
    :type mfe_structures: list
    :param region: The index of a subpart to look for (1 indexed)
    :type region: list
    :returns: Boolean answer if the structure was seen or not
    :rtype: bool
    """
    temp_vienna = Vienna([sequence])
    structure = temp_vienna.mfe(temp=temp, returnstructure=True)[0][1]
    if region != None:
        structure = structure[region[0]-1: region[1]]
        size = len(structure)
    for mfe_structure in mfe_structures:
        if len(mfe_structure) > size:
            raise Exception('MFE_structure larger than region')
        if mfe_structure in structure:
            return True
    return False


def freq_of_structures_appearing(sequence, list_of_target_structures,
                                 num_structures_to_generate, region=None,
                                 temp=37):
    """Function will calculate the frequency that structures appear with a
    certain number of folding simulations

    :param sequence: Sequence of RNA to be sampled
    :type sequence: str
    :param list_of_target_structures: All target structures that are considered
    :type list_of_target_structures: list
    :param num_structures_to_generate: The number of randomly sampled
        structures to generate
    :type num_structures_to_generate: int
    :param region: The index of a subpart to look for (1 indexed)
    :type region: list
    :param temp: Temperature to carry out the simulation
    :type temp: float
    :returns: Frequency of the collection of structures appearing
    :rtype: float
    """
    num = num_structures_to_generate
    tofold = Vienna([sequence])
    number_present = 0

    for dotbracket in tofold.pobabalistic_structures(temp, num):
        if region != None:
            dotbracket = dotbracket[region[0]-1: region[1]]
            size = len(dotbracket)
        for target_structure in list_of_target_structures:
            if len(target_structure) > size:
                raise Exception('target_structure larger than region')
            if target_structure in dotbracket:
                number_present += 1
                break
    return float(number_present)/num_structures_to_generate


def dic_of_sampled_structures(sequence, num_structures_to_generate,
                              region=None, temp=37, returnfreq=False):
    """Function uses vienna sampling to build a large dictionary of randomly
    sampled structures. These structures may be used to approximate a
    partition function of structures

    :param sequence: sequence of RNA to be sampled from
    :type sequence: str
    :param num_structures_to_generate: Number of structures to sample.
    :param region: The index of a subpart to look for (1 indexed)
    :type region: list
    :type num_structures_to_generate: int
    :param returnfreq: Allows for the dict to be in terms of frequency
    :type returnfreq: dict
    :returns: Dictionary of structures with count
    :rtype: Dict
    """
    returndict = {}
    num = num_structures_to_generate
    tofold = Vienna([sequence])
    for dotbracket in tofold.pobabalistic_structures(temp, num):
        if region != None:
            dotbracket = dotbracket[region[0]-1: region[1]]
            size = len(dotbracket)
        try:
            returndict[dotbracket] += 1
        except KeyError:
            returndict[dotbracket] = 1
    if returnfreq:
        total = float(sum(returndict.values()))
        for key in returndict:
            returndict[key]= returndict[key]/total
    return returndict

def list_of_dominant_sampled_structures(sequence, num_structures_to_generate,
                                        threshold=0.01, region=None, temp=37):
    """ Function uses vienna RNA to sample many structures and then return
    The structures that appear at a higher frequency than the threshold
    :param sequence: sequence of RNA to be sampled from
    :type sequence: str
    :param num_structures_to_generate: Number of structures to sample.
    :param region: The index of a subpart to look for (1 indexed)
    :param threshold: The threshold frequency to be included in the list
    :type threshold: float
    :type region: list
    :type num_structures_to_generate: int
    :param returnfreq: Allows for the dict to be in terms of frequency
    :type returnfreq: dict
    :returns: List of all of the structures to over the threshold
    :rtype: list
    """
    dict_of_structures = dic_of_sampled_structures(sequence,
                                                   num_structures_to_generate,
                                                   region=None, temp=37,
                                                   returnfreq=True)
    retrun_list = []
    for key in dict_of_structures:
        if dict_of_structures[key] > 0.01:
            retrun_list.append(key)
    return return_list
