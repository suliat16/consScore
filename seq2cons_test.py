#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs Sept 13 2018

@author: suliat16
"""

import unittest
import os
import seq2conservation
import cons
pst


class Test(unittest.TestCase):
    """
    """
    def setUp(self):
        pass

    @classmethod
    def setUpClass(cls):

        with open((os.getcwd()+os.sep+'example_data'+os.sep+'CDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = cons.OrthologFinder(arabidopsisCDC48A)

    @patch('cons.requests.get')
    def test_call_orthologs(self, requests_mock):
        requests_mock().status_code = 200
        requests_mock().text = 1 #Todo: some response



if __name__ == '__main__':
    unittest.main()
