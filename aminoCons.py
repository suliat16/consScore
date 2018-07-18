#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:45:57 2018

@author: suliat16
"""

import re
from Bio.Align.Applications import TCoffeeCommandline
from biskit.exe import Executor
from biskit.errors import BiskitError
import numpy as np


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
        elif not sequences:
            raise SequenceError("Empty Sequence entered.")
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
        tcoffee_cline = TCoffeeCommandline(infile=file,
                                           output='fasta_seq',
                                           outfile='aligned.aln')
        print(tcoffee_cline)
        tcoffee_cline()

class rate4site(Executor):

    """
    """

    def __init__(self, msa, *args, **kw):
        """
        """
        #Note: This is super inefficient if an alignment is done multiple times,
        #becuase every time rate4site is called, it wipes the last score analysis.
        #TODO: implement an r4s name and use it to name different output files.
        super().__init__(name='rate4site', args='-s %s -o aligned.res'% (msa), tempdir=True,
                             catch_out=1, **kw, verbose=1, cwd= '/tmp')
        
        self.alpha = 0
        self.score_output = self.cwd+'/aligned.res'

    def finish(self):
        """
        """
        super().finish()     
        self.alpha= self.get_alpha(self.score_output)
        self.result= self.read2matrix(self.score_output)

    def isfailed(self):
        """
        Return True if the external program has finished successfully, False
        otherwise
        """
        return 4
#        if self.returncode == 0:
#            return False
#        else:
#            return True

    def fail(self):
        """
        Called if external program has failed
        """
        s = 'Rate4Site failed. Please check the program output in the '+\
            'field `output` of this Rate4Site instance, (eg. `print x.output`)!'
        self.log.add(s)
        raise Rate4SiteError(s)

    def cleanup(self):
      #  super().cleanup()
      return 5
        ## t.tryRemove(self.any_defined_variables)  
        
    def get_alpha(self, r4s):
        """
        Get the alpha parameter of the conservation scores
        Note: This method is especially susceptible to changes in the format of
        the output file
        """
        #TODO: when the reckoning happens, open the r4s file once, call get alpha
        # and read2matrix- none of this multiple opening nonsense
        with open(r4s, 'r') as f:
            contents = f.read()
            splitted = contents.split('\n')
            parameter = ''
            for s in splitted:
                if re.search('alpha parameter', s):
                    parameter = s
                    break
                else: continue
            parameter = rate4site.get_num(parameter)
            return parameter[0]

    @staticmethod
    def get_num(string):
        """
        Takes a string and returns the first number in that string, decimal included
        Args:
            string(str): The string containing a number
        Returns:
            A list of all the numbers contained in that string, as floats
        """
        digits = r"[0-9]*\.?[0-9]+"
        parameter = re.findall(digits, string)
        return list(map(float, parameter))

    def read2matrix(self, file, identity=True, score=True, qqint=False, std=False,
                    msa=False):
        """
        Take the output from rate4site and convert it into a numpy array, mapping
        each conservation score onto its corresponding amino acid
        Args:
            file: The output from the Rate4Site program, version 2.01
        Returns:
            An array, where the entry at each index contains information about 
            the amino acid at that position. Note that the 
        """
        #TODO: Complete- add booleans for parameters in array
        with open(file, 'r') as f:
            contents = f.read()
            residues = rate4site.extract_resi(contents)
            num = self.count_trues(identity, score, qqint, std, msa)
            r2mat = np.empty([1,num])
            for r in residues:
                resi = np.array([])
                #TODO: This is a LOT of repeated code- ask how to make it more 
                #abstract
                if identity == True:
                    amino = rate4site.extract(r,1) 
                    resi=np.append(resi, amino)
                if score == True:
                    conse = rate4site.extract(r, 2)
                    resi=np.append(resi, conse)
                if qqint == True:
                    intqq = rate4site.extract(r, 3)
                    resi=np.append(resi, intqq)
                if std == True:
                    stdev = rate4site.extract(r, 4)
                    resi=np.append(resi, stdev)
                if msa == True:
                    align = rate4site.extract(r, 5)
                    resi=np.append(resi, align)
                resi = resi.reshape((1,num))
                r2mat = np.concatenate((r2mat, resi), axis=0)
            r2mat = np.delete(r2mat, 0, axis=0)
            return r2mat
        
    def count_trues(*args):
        """
        Counts the number of arguments that are true
        """
        i = 0
        for a in args:
            if a is True:  
                i = i+1
        return i

    @staticmethod
    def extract(string, parameter):
        """
        Pull the specified word from a string. Words are defined as 
        Args:
            string (str): The string from which words will be extracted
            parameter (int): The index of the desired word in the split string
        Returns:
            The word at the specified index. Note that the string is split along
            spaces, so any group of characters surrounded with spaces on both sides
            is considered a word. 
        """
        splitted = string.split()
        return splitted[parameter]

    @staticmethod
    def extract_resi(string):
        """
        """
        splitted = string.split('\n')
        residues = []
        for s in splitted:
            if not s.startswith('#'):
                residues.append(s)
        residues = list(filter(lambda x: x is not '', residues))
        return residues 
