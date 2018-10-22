#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 16:58:25 2018

@author: suliat16
"""
import os
import aminoCons as am
import biskit.test
import constool

class test_amino_conservation(biskit.test.BiskitTest):

    """
    Test suite testing the behaviour of the aminoCons module
    """

    TAGS = [biskit.test.EXE, biskit.test.LONG]

    @classmethod
    def setUpClass(cls):
        cls.filepath = os.getcwd() + os.sep + 'example_data'

    def test_clean_alignment(self):
        """Tests to see that alignment file and the folder created are deleted
        after calling clean argument"""
        extra_aln = am.build_alignment(self.filepath +os.sep + 'multiFasta.fasta')
        self.assertTrue(os.path.exists(extra_aln))
        am.clean_alignment(extra_aln, cache=False)
        self.assertFalse(os.path.exists(extra_aln))

    def test_get_num_good(self):
        """Tests that get num can retrieve a single float from a string"""
        digit = constool.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        """Tests that get_num can retrieve multiple floats and integers from a string """
        digit = constool.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        """Tests that a string containing no digits returns no digits"""
        digit = constool.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])

    def test_get_num_negative(self):
        """Tests that get_num can retrieve negative and positive numbers from within a string"""
        digit = constool.get_num("Hello, how are all 1234.56 and -42 of you today")
        self.assertEqual(digit, [1234.56, -42])

    def test_get_alpha_g(self):
        """Tests that get alpha can retrieve the alpha parameter from the file"""
        digit = am.get_alpha(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(digit, 2.83688)
    
    def test_get_alpha_nofile(self):
        """Test that an error is raised if a non-existent file is given to get_alpha"""
        with self.assertRaises(FileNotFoundError):
            am.get_alpha(self.filepath + os.sep + 'TheTree.txt')

    def test_getalpha_badfile(self):
        """Tests that an error is raised if get_alpha is called on an incompatible file"""
        with self.assertRaises(am.Rate4SiteError) as cm:
            am.get_alpha(self.filepath + os.sep + 'Fak2Human.fasta')
        err = cm.exception
        self.assertEqual(str(err), 'File format is not supported')

    def test_r2mat_g(self):
        """Tests that r2mat outputs the correct defaults for the dictionary"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta')
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0:('A', 0.6979),
                    1:('A', 0.6979),
                    2:('C', 0.7026),
                    3:('C', 0.7026),
                    4:('G', 0.1769),
                    5:('G', 0.1769),
                    6:('T', -1.577),
                    7:('T', -1.577)}
        self.assertDictEqual(dictionary, output)

    def test_r2mat_multi(self):
        """Tests that the correct dictionary is output when some parameters are set"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', score=True, qqint=True, gapped=True)
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0:('A', 0.6979, (-1.946, 2.735), '3/3'),
                           1:('A', 0.6979, (-1.946, 2.735), '3/3'),
                           2:('C', 0.7026, (-1.946, 2.735), '3/3'),
                           3:('C', 0.7026, (-1.946, 2.735), '3/3'),
                           4:('G', 0.1769, (-2.332, 1.853), '3/3'),
                           5:('G', 0.1769, (-2.332, 1.853), '3/3'),
                           6:('T', -1.577, (-3.889, -0.7852), '3/3'),
                           7:('T', -1.577, (-3.889, -0.7852), '3/3')}
        self.assertDictEqual(output, dictionary)

    def test_r2mat_score(self):
        """Tests that correct dictionary is output when other parameters are set"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', score=False, qqint=True, gapped=True, std=True)
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0:('A', (-1.946, 2.735), 2.836, '3/3'),
                           1:('A', (-1.946, 2.735), 2.836, '3/3'),
                           2:('C', (-1.946, 2.735), 2.836, '3/3'),
                           3:('C', (-1.946, 2.735), 2.836, '3/3'),
                           4:('G', (-2.332, 1.853), 2.725, '3/3'),
                           5:('G', (-2.332, 1.853), 2.725, '3/3'),
                           6:('T', (-3.889, -0.7852), 2.309, '3/3'),
                           7:('T', (-3.889, -0.7852), 2.309, '3/3')}
        self.assertDictEqual(output, dictionary)

    def test_r2prof(self):
        """tests that rate2profile correctly creates profile collections with the correct output"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', qqint=True, std=True, gapped=True)
        output = test.rate2profile(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(5, len(output))
        self.assertTrue('Amino Acid' in output)
        self.assertTrue('Conservation Score' in output)
        self.assertTrue('Standard Deviation' in output)
        self.assertTrue('QQ interval' in output)
        self.assertTrue('Gapped' in output)
        self.assertTrue(isinstance(output, biskit.ProfileCollection))

    def test_r2prof_default(self):
        """Tests that the default rate2site instantiation produces the right profile"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta')
        output = test.rate2profile(self.filepath + os.sep + 'multiFasta.res')
        self.assertTrue('Amino Acid' in output)
        self.assertTrue('Conservation Score' in output)
        self.assertEqual(2, len(output))
        self.assertTrue(isinstance(output, biskit.ProfileCollection))

    def test_r4s_close(self):
        """Tests to see that close deletes the correct files"""
        r4sobject = am.Rate4Site(self.filepath + os.sep + 'multiFasta.aln')
        r4sobject.run()
        self.assertTrue(os.path.isfile(os.getcwd() + os.sep + 'multiFasta.res'))
        r4sobject.close()
        self.assertFalse(os.path.isfile(os.getcwd() + os.sep + 'multiFasta.res'))
        self.assertFalse(os.path.isdir(os.getcwd() + os.sep + 'multiFasta'))

    @classmethod
    def tearDownClass(cls):
        am.clean_alignment(os.getcwd() + os.sep + 'multiFasta.aln', cache=False)

if __name__ == '__main__':
    biskit.test.localTest()
