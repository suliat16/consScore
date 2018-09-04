#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 09:45:57 2018

@author: suliat16
"""

import os
import re
import warnings
import numpy as np
import biskit.tools as t
from Bio.Align.Applications import TCoffeeCommandline
from biskit.exe import Executor
from biskit.errors import BiskitError

class SequenceError(BiskitError):
    pass

class Rate4SiteError(BiskitError):
    pass

def build_alignment(file):
    """
    Calls the TCoffee program to build an alignment of protein sequences
    Args:
        file: The absolute file path to the collection of protein sequences
    Returns:
        A string detailing the path to the folder containing the alignment file. The
        alignment is output in the current working directory.

    """
    filename = os.path.basename(file)
    filename = filename.split('.')[0]
    tcoffee_cline = TCoffeeCommandline(infile=file,
                                       output='clustalw',
                                       outfile='%s.aln' %(filename))

    old_path = os.getcwd()
    directory = old_path + os.sep + str(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chdir(directory)
        tcoffee_cline()
        os.chdir(old_path)
        return directory
    else:
        os.chdir(directory)
        tcoffee_cline()
        os.chdir(old_path)
        return directory

def clean_alignment(directory):
    """
    """
    filename = os.path.basename(directory)
    filename = filename.split('.')[0]

    if len(os.listdir(directory)) == 0:
        t.tryRemove(directory)
    else:
        t.tryRemove(os.getcwd() + os.sep + '%s.aln' %(filename))
        t.tryRemove(os.getcwd() + os.sep + '%s.dnd' %(filename))

class Rate4Site(Executor):

    """
    Wraps the Rate4Site program. Calling run() executes the program, which creates
    a folder containing the rate4site output information and the tree, and returns
    an array that maps the data onto each amino acid, which by default is the
    conservation score.

    Rate4Site is used for academic purposes. Citation:
    Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of
    site-specific rate-inference methods: Bayesian methods are superior.Mol Biol
    Evol 21: 1781-1791.
    """

    def __init__(self, msa, cwdir=None, cache=True, identity=True, score=True, qqint=False, std=False,
                    gapped=False, **kw):

        aln_file = os.path.basename(msa)
        self.dir_name = aln_file.split('.')[0]
        super().__init__(name='rate4site', args='-s %s -o %s.res'% (msa, self.dir_name),
                         catch_out=1, **kw, tempdir=self.dir_name)
        self.alpha = 0
        if not cwdir:
            self.cwd = os.getcwd() + os.sep + self.dir_name
        else:
            self.cwd = cwdir
        self.score_output = self.cwd + os.sep + '%s.res'% self.dir_name
        self.has_run = False
        self.cache = cache

        self.identity = identity
        self.score= score
        self.qqint = qqint
        self.std = std
        self.gapped = gapped

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
        self.result = self.rate2dict(self.score_output, identity=self.identity, score=self.score,
                                     qqint=self.qqint, std= self.std, gapped = self.gapped)
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
        t.tryRemove(self.cwd + os.sep + 'r4s.res')
        t.tryRemove(self.cwd + os.sep + 'r4sOrig.res')
        super().cleanup()

    def close(self):
        """
        Deletes the output files of rate4site- the alignment tree and the score
        sheet.
        """
        t.tryRemove(self.cwd + os.sep + 'TheTree.txt')
        t.tryRemove(self.cwd + os.sep + '%s.res' %(self.dir_name))
        t.tryRemove(self.tempdir, tree=True)
        self.has_run = False

    def __del__(self):
        """
        Deletes output files for rate4site. When called by garbage collector it also
        deletes the rate4site instance
        """
        if not self.cache:
            self.close()
            t.tryRemove(self.cwd + os.sep + self.dir_name, tree=True)

    @classmethod
    def get_alpha(cls, r4s):
        """
        Get the alpha parameter of the conservation scores
        Args:
            r4s (file): The absolute file path to the alignment file
        Returns:
            The alpha parameter of the conservation score.
        """
        try:
            if os.path.isfile(r4s):
                with open(r4s, 'r') as f:
                    contents = f.read()
                splitted = contents.split(os.linesep)
                for s in splitted:
                    if re.search('alpha parameter', s):
                        parameter = Rate4Site.get_num(s)
                        return parameter[0]
                raise Rate4SiteError('File format is not supported')
            else: 
                raise FileNotFoundError
        except IndexError:
            raise Rate4SiteError('File format is not supported')
        finally:
            warnings.warn("This method is especially susceptible to changes in the format of \
            the output file", Warning)

    @classmethod
    def get_num(cls, string):
        """
        Takes a string and returns the first number in that string, decimal included
        Args:
            string(str): The string containing a number
        Returns:
            A list of all the numbers contained in that string, as floats
        """
        digits = r"-?[0-9]*\.?[0-9]+"
        parameter = re.findall(digits, string)
        return list(map(float, parameter))

    @classmethod
    def rate2dict(cls, r4s, identity=True, score=True, qqint=False, std=False,
                  gapped=False):
        """
        Take the output from rate4site and convert it into a numpy array, mapping
        each conservation score onto its corresponding amino acid.

        Args:
            r4s (str): The absolute filepath to the output file from the Rate4Site program, version 2.01
            If the following parameters are true, the resulting array will contain that information, in the
            order of the arguments
                identity (boolean): The identity of the amino acid (Single letter code) at each position
                score(boolean): The conservation scores. lower value = higher conservation.
                qqint(boolean): QQ-INTERVAL, the confidence interval for the rate estimates. The default interval is 25-75 percentiles
                std(boolean): The standard deviation of hte posterior rate distribution
                gapped(boolean): MSA DATA, the number of aligned sequences having an amino acid (non-gapped) from the overall
                    number of sequences at each position
        Returns:
            An array, where the entry at each index contains information about
            the amino acid at that position.

        The document is parsed by pulling splitting the text from the output file
        using newlines as delimiters, and grabbing only the lines that do not start
        with a # symbol. Whats left are the rows of the table, where each row contains
        information about an amino acid. The rows are then split up depending on
        what information it carries.
        """
        try:
            if os.path.isfile(r4s):
                with open(r4s, 'r') as file:
                    contents = file.read()
                residues = Rate4Site.extract_resi(contents)
                r2dict = {}
                i=0
                for r in residues:
                    r2mat = []
                    aa_data = re.split(r'[\]\[\s,]', r)
                    aa_data = list(filter(lambda x: x is not '', aa_data))
                    if identity:
                        amino = aa_data[1]
                        r2mat.append(amino)
                    if score:
                        conse = float(aa_data[2])
                        r2mat.append(conse)
                    if qqint:
                        intqq = (float(aa_data[3]), float(aa_data[4]))
                        r2mat.append(intqq)
                    if std:
                        stdev = float(aa_data[5])
                        r2mat.append(stdev)
                    if gapped:
                        align = aa_data[6]
                        r2mat.append(align)
                    r2mat = tuple(r2mat)
                    r2dict[i] = r2mat
                    i+=1
                return r2dict
            else:
                raise FileNotFoundError
        finally:
            warnings.warn("This method is especially susceptible to changes in the format of the output file", Warning)

    @classmethod
    def read2matrix(cls, r4s, identity=True, score=True, qqint=False, std=False,
                    gapped=False):
        """
        Take the output from rate4site and convert it into a numpy array, mapping
        each conservation score onto its corresponding amino acid.

        Args:
            r4s (str): The absolute filepath to the output file from the Rate4Site program, version 2.01
            If the following parameters are true, the resulting array will contain that information, in the
            order of the arguments
                identity (boolean): The identity of the amino acid (Single letter code) at each position
                score(boolean): The conservation scores. lower value = higher conservation.
                qqint(boolean): QQ-INTERVAL, the confidence interval for the rate estimates. The default interval is 25-75 percentiles
                std(boolean): The standard deviation of hte posterior rate distribution
                gapped(boolean): MSA DATA, the number of aligned sequences having an amino acid (non-gapped) from the overall
                    number of sequences at each position
        Returns:
            An array, where the entry at each index contains information about
            the amino acid at that position.

        The document is parsed by pulling splitting the text from the output file
        using newlines as delimiters, and grabbing only the lines that do not start
        with a # symbol. Whats left are the rows of the table, where each row contains
        information about an amino acid. The rows are then split up depending on
        what information it carries.
        """
        try:
            if os.path.isfile(r4s):
                with open(r4s, 'r') as file:
                    contents = file.read()
                residues = Rate4Site.extract_resi(contents)
                func_args = [identity, score, qqint, std, gapped]
                dimension = func_args.count(True)
                r2mat = np.empty([1, dimension])
                for r in residues:
                    aa_data = re.split(r'[\]\[\s,]', r)
                    aa_data = list(filter(lambda x: x is not '', aa_data))
                    np_residues = np.array([])
                    if identity:
                        amino = aa_data[1]
                        np_residues = np.append(np_residues, amino)
                    if score:
                        conse = aa_data[2]
                        np_residues = np.append(np_residues, conse)
                    if qqint:
                        intqq = '[{0}, {1}]'.format(aa_data[3] , aa_data[4])
                        np_residues = np.append(np_residues, intqq)
                    if std:
                        stdev = aa_data[5]
                        np_residues = np.append(np_residues, stdev)
                    if gapped:
                        align = aa_data[6]
                        np_residues = np.append(np_residues, align)
                    np_residues = np_residues.reshape((1, dimension))
                    r2mat = np.concatenate((r2mat, np_residues), axis=0)
                r2mat = np.delete(r2mat, 0, axis=0)
                return r2mat
            else:
                raise FileNotFoundError
        finally:
            warnings.warn("This method is especially susceptible to changes in the format of the output file", Warning)

    @classmethod
    def extract_resi(cls, string):
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
            #Remove comments and empty lines
            if s is not '' and not s.startswith('#'):
                #Remove comments
                residues.append(s)
        return residues
