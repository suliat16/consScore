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


if __name__ == '__main__':
    unittest.main()