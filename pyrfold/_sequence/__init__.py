class RNA(Sequence):
    def __init__(self, rna, features):


class Feature():
    def __init__(self, name, start, stop, structurelist, featuretype):
        """
        :param name: name describing the feature
        :type name: str
        :param start: start position of the feature (in 1 indexing)
        :type start: int
        :param stop: stop position of the feature (in 1 indexing)
        :type stop: int
        :param structurelist: a list of the structures with distributin
        :type structurelist: list
        :subtype structurelist: [frequency, sequence, structure, free energy]

        """

class RNAStructure():
    def __init__(self, sequence, dotbracket, freeenergy,
                    frequency=None, calculate=False):
        """
        :param sequence: sequence of the RNA
        :type sequence: str
        :param dotbracket: dot bracket notation of structure
        :type dotbracket: str
        :param freeenergy: free energy of the structure [kcal/mol]
        :type freeenergy: float
        :param calculate: Uses vienna RNA to calculate free energy of struct
        :type calculate: Boolean

        """

