#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:45:57 2018

@author: suliat16
"""

from Bio.Align.Applications import TCoffeeCommandline
from biskit.exe import Executor, RunError
from biskit.errors import BiskitError

class SequenceError(BiskitError):
    pass

class AminoConservation:

    """
    """

    def __init__(self):
        """
        """
        self.sequences = ""
        


    def ortholog_checker(self, sequences):
        """
        Checks to see if input string has a fasta identification line- if not, 
        supplies a default identification line 
        
        Args:
            sequences (str): The sequence(s), in fasta format
        Returns:
            The string, with default identification line, ">Input Sequence\n"
            added if needed
        """
        ##Insert mutated sequence before aligning- make it the reference sequence :)
        if sequences.startswith(">"):
            self.sequences = sequences
        elif not sequences: raise SequenceError("Empty Sequence entered.")
        elif sequences[0].isalpha():
            self.sequences = ">Input Sequence\n"+sequences
        else:
            raise TypeError("Not a FASTA sequence. Please try again")
            
    #TODO: Convert from fasta string with newlines to fasta file with proper formatting

    def build_alignment(self, file):
        """
        """
        #TODO: Either change the implementation, or call it directly
        tcoffee_cline = None
        tcoffee_cline = TCoffeeCommandline(infile= file,
                                           output='fasta_seq',
                                           outfile='aligned.aln')
        print(tcoffee_cline)
        tcoffee_cline()
        
class Rate4SiteWrapper(Executor):
    
    """
    """
    
    def __init__(self, msa, **params):
        """
        """
        # what is a template? Do I need one to call input?
        super().__init__(name='rate4site', tempdir=1, args=msa, f_in=msa, **params, cwd='/tmp')
        
        