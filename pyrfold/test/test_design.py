import unittest
from pyrfold.design import RNA
from pyrfold.utilities import RNA_sequences_complementary, GC_content

class Device_testcase(unittest.TestCase):

    def setUp(self):
        # Need to load in a test case
        self.partnames = ['test1', 'test2', 'test3', 'test4']
        self.partsequences = ['ACAGGU', 'ACCGGGAGA', 'GGGG', 'ACAUAC']
        self.init_obj = RNA.Device(partnamelist=self.partnames,
                                   sequencelist=self.partsequences)

    def test_init(self):
        combined_string = ''
        for sequence in self.partsequences:
            combined_string += sequence

        self.assertEqual(combined_string, str(self.init_obj))

    def test_change_part(self):
        part_to_change = 'test2'
        part_sequence = 'AGAGA'
        new_string = self.partsequences[0] + part_sequence + \
                     self.partsequences[2] + self.partsequences[3]

        self.init_obj.change_part_sequence(partname=part_to_change,
                                           sequence=part_sequence)
        self.assertEqual(new_string, str(self.init_obj))

    def test_add_part(self):
        new_part = 'test5'
        new_sequence = 'AGGG'
        position = 5

        combined_string = ''
        for sequence in self.partsequences:
            combined_string += sequence

        combined_string += new_sequence

        self.init_obj.add_part(partname=new_part, sequence=new_sequence,
                               partposition=position)

        self.assertEqual(combined_string, str(self.init_obj))

    def test_part_sequence(self):
        for partname, partsequence in zip(self.partnames, self.partsequences):
            self.assertEqual(self.init_obj.part_sequence(partname),
                             partsequence)

    def combined_part_sequences(self):
        # Testing individual parts
        for partname, partsequence in zip(self.partnames, self.partsequences):
            self.assertEqual(self.init_obj.combined_part_sequences([partname]),
                             partsequence)
        # Testing multiple parts
        for high_index in range(1, 4):
            partnamelist = self.partnames[0:high_index]
            self.assertEqual(self.init_obj.combined_part_sequences(partnamelist) ,''.join(self.partsequences[0:high_index]))

    def test_remove_part(self):
        self.init_obj.remove_part(partname='test1')

        self.assertEqual(str(self.init_obj), ''.join(self.partsequences[1:5]))


class Helix_testcase(unittest.TestCase):

    def setUp(self):
        self.list_of_helix_helix0 = ['ACGACUGU', 'ACGUACGUAU', 'ACGAUCGGACG']
        self.list_of_helix_helix1 = ['ACAGUCGU', 'AUACGUACGU', 'CGUCCGAUCGU']

    def test_generate_helix(self):
        for helix0, helix1 in zip(self.list_of_helix_helix0,
                                  self.list_of_helix_helix1):

            test_obj = RNA.Helix(helix0=helix0)
            test_obj.generate_helix(based_on_half=0)
            self.assertEqual(test_obj.helixes[1], helix1)
            test_obj.generate_helix(based_on_half=0, random_sub_u_for_c=True)
            self.assertTrue(RNA_sequences_complementary(test_obj.helixes[0],
                                                        test_obj.helixes[1]))

            test_obj = RNA.Helix(helix1=helix1)
            test_obj.generate_helix(based_on_half=1)
            self.assertEqual(test_obj.helixes[0], helix0)
            test_obj.generate_helix(based_on_half=1, random_sub_u_for_c=True)
            self.assertTrue(RNA_sequences_complementary(test_obj.helixes[0],
                                                        test_obj.helixes[1]))


class Unpaired_testcase(unittest.TestCase):

    def setUp(self):
        self.size_ranges = [(1, 10), (2, 4), (7, 8)]
        self.GC_ranges_to_test = [None, (0.1, 0.4), (0.4, 0.7)]

    def test_generate_helix(self):
        for size_range in self.size_ranges:
            for gc_range in self.GC_ranges_to_test:
                test_obj = RNA.Unpaired(size_range=size_range, GC_range=gc_range)
                test_obj.generate_random_sequence()

                self.assertTrue((size_range[0] <= len(str(test_obj))) &
                                (len(str(test_obj)) <= size_range[1]))

                if gc_range:

                    self.assertTrue((gc_range[0] <= GC_content(str(test_obj))) &
                                    (GC_content(str(test_obj)) <= gc_range[1]))


































