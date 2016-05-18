"""
Testing package specifically for all utility functions
"""
import unittest
from pyrfold import utilities


class Random_sequenceTestFunction(unittest.TestCase):

    def setUp(self):
        # Need to load in a test case
        self.lengths_to_test = range(5, 20)
        self.GC_ranges_to_test = [None, (0.1, 0.4), (0.4, 0.7)]
        self.strand_types = ['DNA', 'RNA']
        self.strand_dict = {'RNA': ['A', 'U', 'G', 'C'], 'DNA': ['A', 'T', 'G',
                                                                 'C']}

    def test_random_sequences(self):
        for length in self.lengths_to_test:
            for gc_range in self.GC_ranges_to_test:
                for strand_type in self.strand_types:
                    output = utilities.random_sequence(size=length,
                                                       GC_range=gc_range,
                                                       strand_type=strand_type)

                    # Test Strand type
                    if strand_type == 'DNA':
                        self.assertFalse('U' in output)
                    else:
                        self.assertFalse('T' in output)

                    # GR_Range
                    if gc_range is None:
                        continue
                    else:
                        self.assertTrue((utilities.GC_content(output) <=
                                         gc_range[1]) &
                                        (utilities.GC_content(output) >=
                                        gc_range[0]))
                    # legnth
                    self.assertEqual(len(output), length)


class Reverese_ComplementTestFunction(unittest.TestCase):

    def setUp(self):

        self.testSequencesDNA = [('ACTTGAGGACCACAGAGAGACAGA',
                                  'TCTGTCTCTCTGTGGTCCTCAAGT'),
                                 ('ACTTGAGGACCACAGAGAGACAGAAGACgAGGAgGAGGACGGA',
                                  'TCCGTCCTCCTCCTCGTCTTCTGTCTCTCTGTGGTCCTCAAGT')]

        self.testSequencesRNA = [('AUAGCAUAGUAUCAUUGAUCUAGU',
                                  'ACUAGAUCAAUGAUACUAUGCUAU'),
                                 ('AUUUAUUAUCAGAGCGAUACG',
                                  'CGUAUCGCUCUGAUAAUAAAU')]

    def test_reverse_complement(self):
        # DNA
        for input_seq, output_seq in self.testSequencesDNA:
            self.assertEqual(utilities.reverse_complement(input_seq,
                                                          strand_type='DNA'),
                             output_seq)

        for input_seq, output_seq in self.testSequencesRNA:
            self.assertEqual(utilities.reverse_complement(input_seq,
                                                          strand_type='RNA'),
                             output_seq)


class Index_of_sequenceTestFunction(unittest.TestCase):

    def setUp(self):
        self.sequence = 'ACGACUAGCAGAGCACGGCAGGGACUUCAGUCAGGCAUCAGGCGCGCUCGAGCGA'
        self.part = 'GCAGAGCACGGCA'
        self.index = [8, 20]

    def test_index_of_sequence(self):
        self.assertEqual(self.index, utilities.index_of_sequence(self.sequence,
                                                                 self.part))


class GC_content_TestFunction(unittest.TestCase):

    def setUp(self):
        self.sequence = 'ACCCGUAGUACUGU'

    def test_gc_content(self):
        self.assertEqual(0.5, utilities.GC_content(self.sequence))


if __name__ == '__main__':
    unittest.main()
