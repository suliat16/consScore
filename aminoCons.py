#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:45:57 2018

@author: suliat16
"""

from Bio.Align.Applications import TCoffeeCommandline
from biskit.exe import Executor
from biskit.errors import BiskitError
import biskit.tools as t

class SequenceError(BiskitError):
    pass

class Rate4SiteError(Exception):
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
            raise SequenceError("Not a FASTA sequence. Please try again")
            
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
        super().__init__(name='rate4site', tempdir=True, args=msa, f_in=msa, **params, cwd='/tmp')
    
    def prepare(self):
        """
        """
        super().prepare()
        self.dir_out = self.tempdir
        ##Insert mutated sequence before aligning- make it the reference sequence :)
        ## If I am feeding the sequence directly into Rate4Site wrapper, create
        # an alignment
        # I suppose write a test input file, I GUESS
        
    def finish(self):
        """
        """
        super().finish()
        self.read2matrix()

    def read2matrix(self):
        """
        Take the output from rate4site and convert it into a numpy array, mapping
        each conservation score onto its corresponding amino acid
        """
        #TODO: Complete
        
        
    def isfailed(self):
        """
        Return True if the external program has finished successfully, False 
        otherwise
        """
        if self.returncode == 0: return False
        else: return True
        
    def fail(self):
        """
        Called if external program has failed
        """
        s = 'Rate4Site failed. Please check the program output in the '+\
            'field `output` of this Rate4Site instance, (eg. `print x.output`)!'
        self.log.add(s)
        raise Rate4SiteError(s)
        
    def cleanup(self):
        super().cleanup()
        ## t.tryRemove(self.any_defined_variables)
        
        
        
        
        
        
        
        