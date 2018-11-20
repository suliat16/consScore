#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

from consScore import oma
from consScore import aminoCons
import os
from consScore import constool
import argparse
from requests import exceptions
from biskit.exe import Executor


class OrthoLogo:

    """
    Initializing this class creates an object that stores the parameters of the methods in the pipe
    as fields. Calling the get_logo method runs the pipe, which takes a file in fasta format,
    as well as a motif as input and returns an image of the conserved and depleted acids in the motif.
    """
    def __init__(self, protein, motif, hogs=False, rel_type=None):
        filename = os.path.basename(protein)
        self.name = filename.split('.')[0]
        self.input = protein
        self.motif = motif
        self.hogs = hogs
        self.rel_type = rel_type
        self.sequence = ""
        self.cwd = os.getcwd()
        self.start = 0

    def get_sequence(self):
        """
        Extracts the protein sequence from the input. Modifies the object by taking the sequence input
        and extracting the sequence accordingly
        """
        if os.path.isfile(self.input):
            with open(self.input, "r") as file:
                self.sequence = file.read()
        else:
            raise oma.SequenceError("Cannot open {0}".format(self.input))

    def call_orthologs(self):
        """
        Retrieves the HOG or the orthologs of the entered sequence by querying the OMA database
        """
        ortholog_call = oma.OrthologFinder(self.sequence, rel_type=self.rel_type)
        if self.hogs:
            self.orthologs = ortholog_call.get_HOGs()
        else:
            self.orthologs = ortholog_call.get_orthologs()

        with open("%s.orth" %(self.name), 'w') as o_file:
            o_file.write(self.orthologs)
        return os.getcwd() + os.sep + '%s.orth'%(self.name)

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
        if self.start == -1:
            raise oma.SequenceError("The motif could not be found. Please check the multiple sequence alignment.")
        self.end = self.start + len(self.motif)
        logo = Executor(name='Seq2Logo.py', args='-f %s -o %s -I 2 -C 0 -T 0 -b 0 -c %d-%d'%(msa, self.dir_name,
                                                                                             self.start, self.end))
        logo.run()

    def get_logo(self):
        """
        Queries the OMA database, MAFFT and Seq2logo in sequence
        Returns:
            The path to the folder containing the seq2logo output files.
        """
        self.get_sequence()

        directory = self.cwd + os.sep + 'Seq2Logo'
        if not os.path.isdir(directory):
            os.makedirs(directory)
        os.chdir(directory)

        msa = directory +  os.sep + '%s.aln' % (self.name)
        if os.path.isfile(msa):
            aln = msa
            self.run_seq2logo(aln)
        else:
            orth = self.call_orthologs()
            aln = self.call_alignment(orth)
            self.run_seq2logo(aln)
        os.chdir(self.cwd)
        return directory

if __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("protein", help="Input the path to the file containing the protein")
    parser.add_argument("motif", help="Input the sequence of the motif of interest")
    parser.add_argument("--hogs", action="store_true",
                        help="When specified, the Hierarchical Orthologous Group (HOGs) of the sequence is returned. Defaults to True.")
    parser.add_argument("--reltype", help="Specify the type of orthologs retrieved from OMA, if hogs is false. Defaults to"
                        "None, meaaning all ortholog types are included. Other options are '1:1', 'm:1', 'm:n'")
    args = parser.parse_args()

    logo = OrthoLogo(protein=args.protein, motif=args.motif, hogs=args.hogs, rel_type=args.reltype)
    logo.get_logo()
