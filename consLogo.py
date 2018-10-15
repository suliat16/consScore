#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import oma
import aminoCons
import os
import constool
import biskit.test
from unittest.mock import patch
from requests import RequestException
from biskit.exe import Executor

class OrthoLogo():

    """

    """
    def __init__(self, protein, motif):
        if os.path.isfile(protein):
            filename = os.path.basename(protein)
            self.name = filename.split('.')[0]
        else:
            self.name = "Protein_Sequence"
        self.input = protein
        self.motif = motif
        self.sequence = ""
        self.cwd = os.getcwd()
        self.start = 0

    def get_sequence(self):
        """
        Extracts the protein sequence from the input. Modifies the object by taking the sequence input,
        checking whether the input is a string or a file, and extracting the sequence accordingly
        """
        if os.path.isfile(self.input):
            with open(self.input, "r") as file:
                self.sequence = file.read()
        else:
            self.sequence = self.input

    def call_orthologs(self):
        """
        Retrieves the HOG or the orthologs of the entered sequence by querying the OMA database
        """
        ortholog_call = oma.OrthologFinder(self.sequence)
        try:
            self.orthologs = ortholog_call.get_HOGs()
        except RequestException:
            self.orthologs = ortholog_call.get_orthologs()
        with open("%s.orth" %(self.name), 'w') as o_file:
            o_file.write(self.orthologs)
        return self.cwd + os.sep + '%s.orth'%(self.name)

    def call_alignment(self, orthologs):
        """
        Calls Mafft to create a multiple sequence alignment of the orthologs of the input sequence
        Args:
            Orthologs(str): The path to the file containing the orthologs of the input sequence
        Returns:
            The path to the file containing the alignment of the orthologs
        """
        alignment = aminoCons.build_alignment(orthologs)
        self.alignment = alignment
        return alignment

    def find_motif(self, msa, motif):
        """
        Searches a multiple sequence alignment for the given motif- if found, it returns the index at which the motif
        was found. Note: The indexing refers to the consensus sequence, or to each individual protein, as they are all
        the same length
        Args:
            msa (str): the path to the multiple sequence alignment
            motif (str): the motif to be located in the msa
        Returns:
            The index at which the motif can be found
        """
        with open(msa, "r") as file:
            msa = file.read()
        msa_list = constool.indv_block(msa)
        index = -1
        for prot in msa_list:
            prot = constool.seqnwl_strip(prot)
            prot = constool.get_fasta_sequence(prot)
            index = prot.find(motif)
            if not index == -1:
                break
        return index

    def run_seq2logo(self, msa):
        """Use Executor to call and run Seq2Logo on the selected section of the msa
        Args:
            msa(str): The path to the multiple sequence alignment
            start(int): The starting position of section of interest
            end (int) : The ending position of the section of interest
        Returns:
            The path to the folder containing the output files of Seq2Logo
            """
        aln_file = os.path.basename(msa)
        self.dir_name = aln_file.split('.')[0]
        self.start = self.find_motif(msa, self.motif)
        self.end = self.start + len(self.motif)
        logo = Executor(name='Seq2Logo.py', args='-f %s -o %s -I 2 -C 0 -T 0 -b 0 -c %d-%d'%(msa, self.dir_name,
                                                                                             self.start, self.end), strict=0)

        directory = self.cwd + os.sep + 'Seq2Logo'
        if not os.path.isdir(directory):
           os.makedirs(directory)
        os.chdir(directory)
        logo.run()
        os.chdir(self.cwd)
        return directory

    def get_logo(self):
        """
        Queries the OMA database, MAFFT and Seq2logo in sequence
        Returns:
            The path to the folder containing the seq2logo output files.
        """
        self.get_sequence()
        orth = self.call_orthologs()
        msa = self.call_alignment(orth)
        return self.run_seq2logo(msa)

class ConsLogoTest(biskit.test.BiskitTest):

    """Test suite for the Conservation score logo generator"""

    @classmethod
    def setUpClass(cls):
        cls.path = os.getcwd() + os.sep + 'example_data' + os.sep
        #This path is subject to change, based on where example data is stored for this file
        cls.test_logo = OrthoLogo(cls.path + "atn1seq.fa", "PHHHQHSHIHSHLHLHQ")

    def test_get_sequence_seqinput(self):
        """Tests that given a pure sequence input, that get_sequence saves just the sequence information"""
        seq = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
        test = OrthoLogo(seq, "UUU")
        test.get_sequence()
        self.assertEqual(seq, test.sequence)
        self.assertEqual('Protein_Sequence', test.name)

    def test_get_sequence_fileinput(self):
        """Tests that given a file as input, get_sequence retrieves the sequence inside the file-"""
        test = OrthoLogo(self.path + 'atn1seq.fa', "UUU")
        test.get_sequence()
        self.assertTrue("PASSSAPAPPMRFPYSSSSSSSAAASSSSSSSSSSASPFPASQALPSYPHSFPPPTSLSV" in test.sequence)
        self.assertEqual('atn1seq', test.name)

    @patch('consLogo.oma.OrthologFinder.get_HOGs')
    def test_call_orthologs(self, HOG_mock):
        HOG_mock.return_value = self.path + 'atn1seq.orth'
        test = self.test_logo.call_orthologs()
        self.assertTrue(os.path.isfile(test))
        self.assertTrue(HOG_mock.called)

    @patch('consLogo.aminoCons.build_alignment')
    def test_call_alignment(self, mock_aln):
        """Tests that build_alignment is properly called"""
        mock_aln.return_value = self.path + 'atn1seq.aln'
        test = self.test_logo.call_alignment( self.path + 'atn1seq.orth')
        self.assertTrue(os.path.isfile(test))
        self.assertTrue(mock_aln.called)

    def test_find_motif_present(self):
        """tests that find motif returns the correct index if the motif is present in the alignment"""
        index = self.test_logo.find_motif(self.path + 'atn1seq.aln', "PHHHQHSHIHSHLHLHQ")
        self.assertEqual(1254, index)

    def test_find_motif_absent(self):
        """tests that find_motif returns -1 if the motif cannot be found in the alignment"""
        index = self.test_logo.find_motif(self.path + 'atn1seq.aln', "THEMOTIFISINTHEMOTIVE")
        self.assertEqual(-1, index)

    # @patch('consLogo.OrthoLogo.find_motif')
    # def test_run_seq2logo(self, mock_motif):
    #     mock_motif.return_value = 1254
    #     test = self.test_logo.run_seq2logo(self.path + "atn1seq.aln")
    #     self.assertTrue(os.path.isdir(test))
    #     self.assertTrue(os.path.isfile(os.getcwd() + os.sep + 'Seq2Logo' + os.sep + 'atn1seq_freq.mat'))
    #     self.assertTrue(mock_motif.called)

    # @patch('consLogo.OrthoLogo.get_sequence')
    # @patch('consLogo.OrthoLogo.call_orthologs')
    # @patch('consLogo.OrthoLogo.call_alignment')
    # @patch('consLogo.OrthoLogo.run_seq2logo')
    # def test_get_logo(self, mock_run, mock_aln, mock_orth, mock_seq):
    #     """Tests that get_logo calls the right methods in sequence"""
    #     mock_seq.return_value = self.path + 'atn1seq.fa'
    #     mock_orth.return_value = self.path + 'atn1seq.orth'
    #     mock_aln.return_value = self.path + 'atn1seq.aln'
    #     mock_run.return_value = os.getcwd() + os.sep + 'Seq2Logo'
    #     test = OrthoLogo(self.path + 'atn1seq.fa', 'PHHHQHSHIHSHLHLHQ')
    #     test.get_logo()
    #     self.assertTrue(mock_run.called)
    #     self.assertTrue(mock_aln.called)
    #     self.assertTrue(mock_orth.called)
    #     self.assertTrue(mock_run.called)
    #     self.assertTrue(os.path.isdir(test))

if __name__== '__main__':
    biskit.test.localTest()
