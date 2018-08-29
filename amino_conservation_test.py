#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 16:58:25 2018

@author: suliat16
"""
import os
import unittest
import aminoCons as am
import numpy as np

class test_amino_conservation(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(cls):
        cls.filepath = os.getcwd() + os.sep + 'example_data'
        cls.r4sobject = am.Rate4Site(cls.filepath + os.sep + 'multiFasta.aln')
        cls.rt2mat = cls.r4sobject.run()
        cls.test = am.build_alignment(cls.filepath + os.sep + 'multiFasta.fasta')

    def test_tcof_output(self):
        self.assertTrue(os.path.isdir(self.test))
        with open(self.test + os.sep + 'multiFasta.aln', 'r') as f:
            contents = f.read()
            self.assertIn('AACCGGTT', contents)

    def test_get_num_good(self):
        digit = am.Rate4Site.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        digit = am.Rate4Site.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        digit = am.Rate4Site.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])

    def test_get_num_negative(self):
        digit = am.Rate4Site.get_num("Hello, how are all 1234.56 and -42 of you today")
        self.assertEqual(digit, [1234.56, -42])

    def test_get_alpha_g(self):
        digit = self.r4sobject.get_alpha(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(digit, 2.83688)
    
    def test_get_alpha_nofile(self):
        with self.assertRaises(FileNotFoundError) as cm:
            self.r4sobject.get_alpha(self.filepath + os.sep + 'TheTree.txt')

    def test_getalpha_badfile(self):
        with self.assertRaises(am.Rate4SiteError) as cm:
            self.r4sobject.get_alpha(self.filepath + os.sep + 'Fak2Human.fasta')
        err = cm.exception
        self.assertEqual(str(err), 'File format is not supported')

    def test_r2mat_g(self):
        output = self.r4sobject.read2matrix(self.filepath + os.sep + 'multiFasta.res')
        matrix = np.array([['A', '0.6979'],
                          ['A', '0.6979'],
                          ['C', '0.7026'],
                          ['C', '0.7026'],
                          ['G', '0.1769'],
                          ['G', '0.1769'],
                          ['T', '-1.577'],
                          ['T', '-1.577']])
        np.testing.assert_array_equal(matrix, output)

    def test_r2mat_multi(self):
        output = self.r4sobject.read2matrix(self.filepath + os.sep + 'multiFasta.res', score=True, qqint=True, msa=True)
        matrix = np.array([['A', '0.6979', '[-1.946, 2.735]', '3/3'],
                             ['A', '0.6979', '[-1.946, 2.735]', '3/3'],
                             ['C', '0.7026', '[-1.946, 2.735]', '3/3'],
                             ['C', '0.7026', '[-1.946, 2.735]', '3/3'],
                             ['G', '0.1769', '[-2.332, 1.853]', '3/3'],
                             ['G', '0.1769', '[-2.332, 1.853]', '3/3'],
                             ['T', '-1.577', '[-3.889, -0.7852]', '3/3'],
                             ['T', '-1.577', '[-3.889, -0.7852]', '3/3']])
        np.testing.assert_array_equal(output, matrix)

    def test_r2mat_score(self):
        output = self.r4sobject.read2matrix(self.filepath + os.sep + 'multiFasta.res', score=False, qqint=True, msa=True, std=True)
        matrix = np.array([['A', '[-1.946, 2.735]', '2.836', '3/3'],
                           ['A', '[-1.946, 2.735]', '2.836', '3/3'],
                           ['C', '[-1.946, 2.735]', '2.836', '3/3'],
                           ['C', '[-1.946, 2.735]', '2.836', '3/3'],
                           ['G', '[-2.332, 1.853]', '2.725', '3/3'],
                           ['G', '[-2.332, 1.853]', '2.725','3/3'],
                           ['T', '[-3.889, -0.7852]', '2.309', '3/3'],
                           ['T', '[-3.889, -0.7852]', '2.309', '3/3']])
        np.testing.assert_array_equal(output, matrix)

    def test_r4s_close(self):
        self.r4sobject.close()
        self.assertFalse(os.path.isfile(os.getcwd() + os.sep + 'multiFasta.res'))
        self.assertFalse(os.path.isdir(os.getcwd() + os.sep + 'multiFasta'))

        #TODO: Figure out a way to close the example folder created in testing

    def test_del_garbage(self):
        """ Tests to see if the files are deleted after the pointer to them is gone
        """
        pointer = am.Rate4Site(os.getcwd() + os.sep + 'multiFasta.aln')
        pointer.run()
        pointer = 64
        self.assertFalse(os.path.exists(os.getcwd() + os.sep + 'multiFasta'))

    def test_clean_argument(self):
        """Tests to see that alignment file and the folder created are deleted
        after calling clean argument"""
      #  am.clean_alignment(self.test)
       # self.assertFalse(os.path.isdir(self.test))

    @classmethod
    def tearDownClass(cls):
        cls.r4sobject.__del__()
        am.clean_alignment(cls.test)

if __name__ == '__main__':
    unittest.main()
