#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs Sept 13 2018

@author: suliat16
"""

import biskit.test
import os
from consScore import aminoCons
from consScore import seq2conservation as sq
from requests import exceptions
from unittest.mock import patch


class Test(biskit.test.BiskitTest):
    """
    Test suite for the seq2conservation pipeline
    """

    TAGS = [biskit.test.NORMAL, biskit.test.LONG]

    @classmethod
    def setUpClass(cls):
        cls.cwd = os.getcwd()

        with open((os.getcwd()+os.sep+'example_data'+os.sep+'CDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = sq.ConservationPipe(arabidopsisCDC48A, cache=False)

        with open((os.getcwd()+os.sep+'example_data'+os.sep+'multiFasta.fasta'), 'r') as file:
            cls.ex_seq = file.read()

    @patch('seq2conservation.oma.OrthologFinder.get_HOGs')
    def test_call_orthologs(self, HOG_mock):
        """tests that call_orthologs correctly calls methods to make an output file of 'orthologs'"""
        HOG_mock.return_value = self.ex_seq
        tester = self.CDC48A.call_orthologs()
        self.assertTrue(os.path.isfile(tester))
        self.assertTrue('Protein_Sequence.orth' in tester)

    @patch('consLogo.oma.OrthologFinder.get_HOGs')
    @patch('consLogo.oma.OrthologFinder.get_orthologs')
    def test_call_orthologs_except(self, orth_mock, HOG_mock):
        """Tests that call_orthologs properly calls oma.get_orthologs if get_HOGs fails"""
        HOG_mock.side_effect = exceptions.RequestException('There was an issue querying the database. Status code 401')
        orth_mock.return_value = self.ex_seq
        test = self.CDC48A.call_orthologs()
        self.assertTrue(os.path.isfile(test))
        self.assertTrue(HOG_mock.called)
        self.assertTrue(orth_mock.called)
        self.assertTrue('Protein_Sequence.orth' in test)

    @patch('seq2conservation.aminoCons.build_alignment')
    def test_call_alignment(self, mock_aln):
        """tests that call_alignment calls the correct methods and generates the correct output"""
        mock_aln.return_value = os.getcwd()+os.sep+'example_data'+os.sep + 'multiFasta.aln'
        tester = self.CDC48A.call_alignment(os.getcwd()+os.sep+'example_data'+os.sep+'multiFasta.fasta')
        self.assertTrue(mock_aln.called)
        self.assertTrue(os.path.isfile(tester))
        self.assertTrue('multiFasta.aln' in tester)
        aminoCons.clean_alignment(tester, cache=False)

    @patch('seq2conservation.aminoCons.get_alpha')
    def test_call_rate4site(self, mock_score):
        """tests that call_rate4site calls the correct methods and generates the correct output"""
        mock_score.return_value = 2.83688
        tester= self.CDC48A.call_rate4site(os.getcwd()+os.sep+'example_data'+os.sep + 'multiFasta.aln')
        self.assertTrue(mock_score.called)
        self.assertEqual(2.83688, tester)

    @patch('seq2conservation.oma.OrthologFinder.get_HOGs')
    @patch('seq2conservation.aminoCons.build_alignment')
    @patch('seq2conservation.aminoCons.Rate4Site.run')
    def test_pipe(self, mock_run, mock_aln, mock_hog):
        """tests that the pipe calls the correct methods and generates the correct output"""
        mock_hog.return_value = self.ex_seq
        mock_aln.return_value = os.getcwd()+os.sep+'example_data'+os.sep + 'multiFasta.aln'
        mock_run.return_value = {0:('A', 0.6979),
                    1:('A', 0.6979),
                    2:('C', 0.7026),
                    3:('C', 0.7026),
                    4:('G', 0.1769),
                    5:('G', 0.1769),
                    6:('T', -1.577),
                    7:('T', -1.577)}
        tester = self.CDC48A.pipe()
        self.assertTrue(mock_run.called)
        self.assertTrue(mock_aln.called)
        self.assertTrue(mock_hog.called)
        self.assertTrue(type(tester), dict)

    @classmethod
    def tearDownClass(cls):
        os.remove(os.getcwd()+ os.sep + 'Protein_Sequence.orth')

if __name__ == '__main__':
    biskit.test.localTest()
