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

    def setUp(self):
        self.filepath = os.getcwd() + os.sep + 'example_data' 
        self.testscore = am.Rate4Site(self.filepath +os.sep + 'multiFasta.aln')
        self.rt2mat = self.testscore.run()

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
    
    def test_tcof_output(self):
        test = am.build_alignment(self.filepath + os.sep + 'multiFasta.fasta')
        self.assertTrue(os.path.isdir(test))
        with open(test + os.sep + 'multiFasta.fasta', 'r') as f:
            contents = f.read()
            self.assertIn('AACCGGTT', contents)

    def test_get_alpha_g(self):
        digit = self.testscore.get_alpha(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(digit, 2.83688)
    
    def test_get_alpha_nofile(self):
        with self.assertRaises(FileNotFoundError) as cm:
            self.testscore.get_alpha(self.filepath + os.sep + 'TheTree.txt')

    def test_getalpha_badfile(self):
        with self.assertRaises(am.Rate4SiteError) as cm:
            self.testscore.get_alpha(self.filepath + os.sep + 'Fak2Human.fasta')
        err = cm.exception
        self.assertEqual(str(err), 'File format is not supported')

    def test_r2mat_g(self):
        output = self.testscore.read2matrix(self.filepath + os.sep + 'multiFasta.res')
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
        output = self.testscore.read2matrix(self.filepath + os.sep + 'multiFasta.res', score=True, qqint=True, msa=True)
        matrix = np.array([['A', '0.6979', '[-1.946, 2.836]', '3/3'],
                             ['A', '0.6979', '[-1.946, 2.836]', '3/3'],
                             ['C', '0.7026', '[-1.946, 2.836]', '3/3'],
                             ['C', '0.7026', '[-1.946, 2.836]', '3/3'],
                             ['G', '0.1769', '[-2.332, 2.725]', '3/3'],
                             ['G', '0.1769', '[-2.332, 2.725]', '3/3'],
                             ['T', '-1.577', '[-3.889,-0.7852]', '3/3'],
                             ['T', '-1.577', '[-3.889,-0.7852]', '3/3']])
        print(output)
        np.testing.assert_array_equal(output, matrix)

    def test_r2mat_score(self):
        output = self.testscore.read2matrix(self.filepath + os.sep + 'multiFasta.res', score=True, qqint=True, msa=True)
        print(output)


    #TODO: Test r2mat, please

if __name__ == '__main__':
    unittest.main()
