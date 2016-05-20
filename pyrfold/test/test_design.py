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

    def test_combined_part_sequences(self):
        # Testing individual parts
        for partname, partsequence in zip(self.partnames, self.partsequences):
            self.assertEqual(self.init_obj.combined_part_sequences([partname]),
                             partsequence)
        # Testing multiple parts
        for high_index in range(1, 4):
            partnamelist = self.partnames[0:high_index]
            self.assertEqual(self.init_obj.combined_part_sequences(partnamelist) ,''.join(self.partsequences[0:high_index]))

    def test_create_kinefold_submisison_object(self):
        # First testing the 'all' case
        device_name = 'test_device'
        partcontexttofold = 'all'
        partstofold = None

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold)

        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'

        self.assertEqual(expected_name, test_obj.name)
        self.assertEqual(self.init_obj.sequence, test_obj.sequence)

        device_name = 'test_device'
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold)
        expected_sequence = self.init_obj.combined_part_sequences(
            partcontexttofold)
        self.assertEqual(expected_sequence, test_obj.sequence)
        self.assertEqual(1, test_obj.windowstart)
        self.assertEqual(len(expected_sequence), test_obj.windowstop)

        # TESTING error cases for devices
        fiveprimeshift = -20
        fiveprimerefpart = None
        threeprimeshift = False
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        self.assertRaises(IndexError,
                          self.init_obj.create_kinefold_submission_object,
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        fiveprimeshift = 50
        fiveprimerefpart = None
        threeprimeshift = False
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        self.assertRaises(IndexError,
                          self.init_obj.create_kinefold_submission_object,
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        fiveprimeshift = -3
        fiveprimerefpart = None
        threeprimeshift = 30
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        self.assertRaises(IndexError,
                          self.init_obj.create_kinefold_submission_object,
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        fiveprimeshift = -3
        fiveprimerefpart = None
        threeprimeshift = -50
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        self.assertRaises(IndexError,
                          self.init_obj.create_kinefold_submission_object,
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        # Testing naming conventions

        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None
        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'+\
            '#five_prime_shift&2#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)
        self.assertEqual(expected_name, test_obj.name)

        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = 'all'
        partstofold = None

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        # Testing that the start and stop are correct
        # Requested five 2 and three -3
        # ACAGGU-ACCGGGAGA-GGGG-ACAUAC
        # --+---------------------+---
        expected_start = 3
        expected_stop = 22
        self.assertEqual(expected_start, test_obj.windowstart)
        self.assertEqual(expected_stop, test_obj.windowstop)

        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        # Testing that the start and stop are correct
        # Requested five 2 and three -3
        # ACCGGGAGA-GGGG-ACAUAC
        # --+------&----&--+---
        expected_start = 3
        expected_stop = 16
        self.assertEqual(expected_start, test_obj.windowstart)
        self.assertEqual(expected_stop, test_obj.windowstop)

        fiveprimeshift = -4
        fiveprimerefpart = 'test2'
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = 'all'
        partstofold = None

        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'+\
            '#five_prime_shift&-4#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        self.assertEqual(expected_name, test_obj.name)

        # Testing that the start and stop are correct
        # Requested five 2 and three -3
        # ACAGGU-ACCGGGAGA-GGGG-ACAUAC
        # --+---------------------+---
        expected_start = 3
        expected_stop = 22
        self.assertEqual(expected_start, test_obj.windowstart)
        self.assertEqual(expected_stop, test_obj.windowstop)

        fiveprimeshift = -4
        fiveprimerefpart = 'test3'
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None

        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'+\
            '#five_prime_shift&-4#rel_five_prime_part&test3' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)

        self.assertEqual(expected_name, test_obj.name)

        # Testing that the start and stop are correct
        # Requested five -4 and three -3
        # ACCGGGAGA-GGGG-ACAUAC
        # -----+-----------+---
        expected_start = 6
        expected_stop = 16
        self.assertEqual(expected_start, test_obj.windowstart)
        self.assertEqual(expected_stop, test_obj.windowstop)

        # Testing folding parts
        fiveprimeshift = -4
        fiveprimerefpart = 'test2'
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = 'all'
        partstofold = ['test3']

        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'+\
            '#five_prime_shift&-4#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4' + \
            '#part_to_fold&test3'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)
        self.assertEqual(expected_name, test_obj.name)
        # ACAGGU-ACCGGGAGA-GGGG-ACAUAC
        # --+---------------------+---
        partstart = 16
        partstop = 16+3
        self.assertEqual(test_obj.partstartstoplist[0], [partstart, partstop])

        fiveprimeshift = -4
        fiveprimerefpart = 'test2'
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = 'all'
        partstofold = ['test3', 'test4']

        expected_name = 'dev&'+device_name+'#pol&30#fold_time_after&1'+\
            '#five_prime_shift&-4#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4' + \
            '#part_to_fold&test3#part_to_fold&test4'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart)
        self.assertEqual(expected_name, test_obj.name)
        # ACAGGU-ACCGGGAGA-GGGG-ACAUAC
        # --+---------------------+---
        partstart = 16
        partstop = 16+3
        self.assertEqual(test_obj.partstartstoplist[0], [partstart, partstop])
        partstart = 20
        partstop = 20+5
        self.assertEqual(test_obj.partstartstoplist[1], [partstart, partstop])

        # Testing pseudoknots and entaglemtns
        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None
        polrate = 30
        foldtimeafter = 10
        experimenttype = 1
        pseudoknots = 0
        entanglements = 0

        expected_name = 'dev&'+device_name+'#anneal_time&10' + \
            '#five_prime_shift&2#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart,
                                                    polrate, foldtimeafter,
                                                    experimenttype,
                                                    pseudoknots,
                                                    entanglements)

        self.assertEqual(expected_name, test_obj.name)
        self.assertEqual(experimenttype, test_obj.experimenttype)
        self.assertEqual(pseudoknots, test_obj.pseudoknots)
        self.assertEqual(entanglements, test_obj.entanglements)

        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None
        polrate = 30
        foldtimeafter = 10
        experimenttype = 1
        pseudoknots = 1
        entanglements = 0

        expected_name = 'dev&'+device_name+'#anneal_time&10' + \
            '#five_prime_shift&2#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4' + \
            '#pseudoknots&True'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart,
                                                    polrate, foldtimeafter,
                                                    experimenttype,
                                                    pseudoknots,
                                                    entanglements)

        self.assertEqual(expected_name, test_obj.name)
        self.assertEqual(experimenttype, test_obj.experimenttype)
        self.assertEqual(pseudoknots, test_obj.pseudoknots)
        self.assertEqual(entanglements, test_obj.entanglements)

        fiveprimeshift = 2
        fiveprimerefpart = None
        threeprimeshift = -3
        threeprimerefpart = None
        partcontexttofold = ['test2', 'test3', 'test4']
        partstofold = None
        polrate = 30
        foldtimeafter = 10
        experimenttype = 1
        pseudoknots = 1
        entanglements = 1

        expected_name = 'dev&'+device_name+'#anneal_time&10' + \
            '#five_prime_shift&2#rel_five_prime_part&test2' + \
            '#three_prime_shift&-3#rel_three_prime_part&test4' + \
            '#pseudoknots&True' + '#entanglements&True'

        test_obj = self.init_obj.create_kinefold_submission_object(
                                                    device_name,
                                                    partcontexttofold,
                                                    partstofold,
                                                    fiveprimeshift,
                                                    fiveprimerefpart,
                                                    threeprimeshift,
                                                    threeprimerefpart,
                                                    polrate, foldtimeafter,
                                                    experimenttype,
                                                    pseudoknots,
                                                    entanglements)

        self.assertEqual(expected_name, test_obj.name)
        self.assertEqual(experimenttype, test_obj.experimenttype)
        self.assertEqual(pseudoknots, test_obj.pseudoknots)
        self.assertEqual(entanglements, test_obj.entanglements)

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
        self.size_ranges = [(5, 10), (12, 20), (7, 8)]
        self.GC_ranges_to_test = [None, (0.1, 0.4), (0.4, 0.7)]

    def test_generate_random_sequence(self):
        for size_range in self.size_ranges:
            for gc_range in self.GC_ranges_to_test:

                test_obj = RNA.Unpaired(size_range=size_range,
                                        GC_range=gc_range)
                test_obj.generate_random_sequence()


                self.assertTrue((size_range[0] <= len(str(test_obj))) &
                                (len(str(test_obj)) <= size_range[1]))

                if gc_range:

                    self.assertTrue((gc_range[0] <= GC_content(str(test_obj))) &
                                    (GC_content(str(test_obj)) <= gc_range[1]))





























