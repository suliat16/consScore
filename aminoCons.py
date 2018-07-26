#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:45:57 2018

@author: suliat16
"""

import os
import re
import numpy as np
import biskit.tools as t
from Bio.Align.Applications import TCoffeeCommandline
from biskit.exe import Executor
from biskit.errors import BiskitError

class SequenceError(BiskitError):
    pass

class Rate4SiteError(Exception):
    pass


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
        self.sequences = ">Input Sequence"+ os.linesep +sequences
    else:
        raise SequenceError("Not a FASTA sequence. Please try again")
#TODO: Convert from fasta string with newlines to fasta file with proper formatting

def build_alignment(self, file):
    """
    Calls the TCoffee program to build an alignment of protein sequences
    Args:
        file: The absolute file path to the collection of protein sequences
    Returns: 
        
    """
    #TODO: Either change the implementation, or call it directly
    tcoffee_cline = None
    tcoffee_cline = TCoffeeCommandline(infile=file,
                                       output='fasta_seq',
                                       outfile='aligned.aln')
    tcoffee_cline()

class rate4site(Executor):

    """
    Wraps the Rate4Site program. Calling run() executes the program, which creates
    a folder containing the rate4site output information and the tree, and returns
    an array that maps the data onto each amino acid, which by default is the
    conservation score. 
    
    Rate4Site is used here for academic purposes. Citation: 
    Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of 
    site-specific rate-inference methods: Bayesian methods are superior.Mol Biol 
    Evol 21: 1781-1791.
    """

    def __init__(self, msa, cwdir=None, *args, **kw):
        
        aln_file = os.path.basename(msa)
        self.dir_name = (aln_file.split('.'))[0]
        super().__init__(name='rate4site', args='-s %s -o %s.res'% (msa, self.dir_name),
                         catch_out=1, **kw, tempdir=self.dir_name)
        self.alpha = 0
        if not cwdir:
            self.cwd = os.getcwd() + os.sep + self.dir_name
        self.score_output = self.cwd + os.sep + '%s.res'%(self.dir_name)
        self.keep_tempdir = True
        self.has_run = False


    def run(self):
        """
        Calls the executor run method if it is a first run, otherwise just calls 
        the post execution methods on the cached files
        """
        if self.has_run:
            self.finish()
            return self.result
        else:
            return super().run()


    def finish(self):
        """
        Overwrites Executor method. Called when the program is done executing.
        """
        super().finish()
        self.alpha = self.get_alpha(self.score_output)
        self.result = self.read2matrix(self.score_output)
        self.keep_tempdir = True
        self.has_run = True

    def isfailed(self):
        """
        Overwrites Executor method. Return True if the external program has finished 
        successfully, False otherwise
        """
        ret = None
        if self.returncode == 0:
            ret = False
        else:
            ret = True
        return ret

    def fail(self):
        """
        Overwrites Executor method. Called if external program has failed
        """
        s = 'Rate4Site failed. Please check the program output in the '+\
            'field `output` of this Rate4Site instance, (eg. `print x.output`)!'
        self.log.add(s)
        raise Rate4SiteError(s)

    def cleanup(self):
        """
        Overwrites Executor method. Cleans up files created during program execution.
        """
        #t.tryRemove(self.cwd + os.sep + 'TheTree.txt')
        t.tryRemove(self.cwd + os.sep + 'r4s.res')
        t.tryRemove(self.cwd + os.sep + 'r4sOrig.res')
        super().cleanup()

    def close(self):
        """
        Deletes the output files of rate4site- the alignment tree and the score
        sheet. 
        """
        t.tryRemove(self.cwd + os.sep + 'TheTree.txt')
        t.tryRemove(self.tempdir, tree=True)
        self.has_run = False

    def __del__(self):
        """
        Deletes output files for rate4stie. When called by garbage collector it also
        deletes the rate4site instance
        """
        t.tryRemove(self.cwd + os.sep + 'TheTree.txt')
        t.tryRemove(self.tempdir, tree=True)

    def get_alpha(self, r4s):
        """
        Get the alpha parameter of the conservation scores
        Note: This method is especially susceptible to changes in the format of
        the output file
        """
        with open(r4s, 'r') as f:
            contents = f.read()
            splitted = contents.split(os.linesep)
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
        each conservation score onto its corresponding amino acid.

        Args:
            file: The output from the Rate4Site program, version 2.01
        Returns:
            An array, where the entry at each index contains information about
            the amino acid at that position.

        The document is parsed by pulling splitting the text from the output file
        using newlines as delimiters, and grabbing only the lines that do not start
        with a # symbol. Whats left are the rows of the table, where each row contains
        information about an amino acid. The rows are then split up depending on
        what information it carries.
        """
        with open(file, 'r') as f:
            contents = f.read()
            residues = rate4site.extract_resi(contents)
            num = self.count_trues(identity, score, qqint, std, msa)
            r2mat = np.empty([1, num])
            for r in residues:
                resi = np.array([])
                #TODO: This is a LOT of repeated code- ask how to make it more
                #abstract
                if identity == True:
                    amino = rate4site.extract(r, 1)
                    resi = np.append(resi, amino)
                if score == True:
                    conse = rate4site.extract(r, 2)
                    resi = np.append(resi, conse)
                if qqint == True:
                    intqq = rate4site.extract(r, 3)
                    resi = np.append(resi, intqq)
                if std == True:
                    stdev = rate4site.extract(r, 4)
                    resi = np.append(resi, stdev)
                if msa == True:
                    align = rate4site.extract(r, 5)
                    resi = np.append(resi, align)
                resi = resi.reshape((1, num))
                r2mat = np.concatenate((r2mat, resi), axis=0)
            r2mat = np.delete(r2mat, 0, axis=0)
            return r2mat

    def count_trues(self, *args):
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
        Pull the specified word from a string.
        """
        splitted = string.split()
        return splitted[parameter]

    @staticmethod
    def extract_resi(string):
        """
        Grabs the lines of the table that correspond to amino acid data, and puts
        them in a list.
        Args:
            string(str): The contents of the rate4site file
        Returns:
            The rows of the amino acid table as a list of strings.
        """
        splitted = string.split(os.linesep)
        residues = []
        for s in splitted:
            if not s.startswith('#'):
                residues.append(s)
        residues = list(filter(lambda x: x is not '', residues))
        return residues
