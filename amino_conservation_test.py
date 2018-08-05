#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 16:58:25 2018

@author: suliat16
"""
import os
import unittest
import aminoCons as am

class test_amino_conservation(unittest.TestCase):

    """
    """

    def setUp(self):
        self.filepath = os.getcwd() + os.sep + 'example_data' 
        self.testscore = am.Rate4Site(self.filepath +os.sep + 'multiFasta.aln')
        self.rt2mat = self.testscore.run()

    def test_ortho_chk(self):
        output = am.ortholog_checker("Hello World")
        self.assertEqual(output, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        with self.assertRaises(am.SequenceError) as cm:
            am.ortholog_checker("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        with self.assertRaises(am.SequenceError) as sm:
            am.ortholog_checker("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a FASTA sequence. Please try again")

    def  test_ortho_already_fasta(self):
        output = am.ortholog_checker(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(output, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")

    def test_get_num_good(self):
        digit = am.Rate4Site.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        digit = am.Rate4Site.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        digit = am.Rate4Site.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])
    
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



if __name__ == '__main__':
    unittest.main()
