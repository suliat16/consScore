#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 16:58:25 2018

@author: suliat16
"""

import unittest
import aminoCons as am

class test_amino_conservation(unittest.TestCase):

    """
    """

    def setUp(self):
        self.tester = am.AminoConservation()

    def test_ortho_chk(self):
        self.tester.ortholog_checker("Hello World")
        self.assertEqual(self.tester.sequences, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        with self.assertRaises(am.SequenceError) as cm:
            self.tester.ortholog_checker("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        with self.assertRaises(am.SequenceError) as sm:
            self.tester.ortholog_checker("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a FASTA sequence. Please try again")

    def  test_ortho_already_fasta(self):
        self.tester.ortholog_checker(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(self.tester.sequences, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")

    def test_get_num_good(self):
        digit = am.rate4site.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        digit = am.rate4site.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        digit = am.rate4site.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])
        
    #TODO: Test get_alpha

if __name__ == '__main__':
    unittest.main()
