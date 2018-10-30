#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline which retrieves the conservation score of the amino acids in sequence
given an input sequence or fasta file. Intakes the sequence or fasta string of a protein using the standard, single
letter alphabet, gets the orthologs from an external (online) database, the OMA browser. It then generates a
multiple sequence alignment (MSA) using Mafft, and calculates the conservation score of each
amino acid at each position using Rate4Site, with respect to the entered sequence.

The steps of the pipeline are as follows:
Sequence --> Orthologs
Orthologs --> Alignment
Alignment --> Conservation scores

Citations
OMA database:
Altenhoff A et al., The OMA orthology database in 2018: retrieving evolutionary relationships among
all domains of life through richer web and programmatic interfaces Nucleic Acids Research, 2018,
46 (D1): D477-D485 (doi:10.1093/nar/gkx1019).

Mafft:
*Katoh, Misawa, Kuma and Miyata (Nucleic Acids Res. 30:3059-3066, 2002) MAFFT: a novel method for rapid
multiple sequence alignment based on fast Fourier transform (describes the FFT-NS-1, FFT-NS-2 and FFT-NS-i
strategies)


Rate4Site:
Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference methods:
Bayesian methods are superior. Mol Biol Evol 21: 1781-1791.
"""

from consScore import oma
from consScore import aminoCons
import os
from biskit.errors import BiskitError
from requests import RequestException


class PipelineError(BiskitError):
    pass


class ConservationPipe:

    """
    Initializing this class creates an object that stores the parameters of the methods in the pipeline as
    fields. Calling the pipe method runs the pipe, which takes a sequence or file as input and returns various
    conservation related scores of each amino acid.
    """

    def __init__(self, sequence, name=None, cache=True, profile=True, identity=True, score=True, qqint=False, std=False,
                gapped=False):
        """
        Args:
            sequence (str): The sequence of the protein of interest, or the filepath of the fasta file containing
            the sequence of interest
            name (str): The name of the output alignment file. If a file is given, the name of the fasta file is used.
            Defaults to Protein Sequence
            cache(boolean): When true, generated MSA is saved in a folder called protein sequences. When false, all
            files are deleted. If the following parameters are true, the output dictionary will contain that information,
            in the order of the arguments
                identity (boolean): The identity of the amino acid (Single letter code) at each position
                score(boolean): The conservation scores. lower value = higher conservation.
                qqint(boolean): QQ-INTERVAL, the confidence interval for the rate estimates. The default interval is 25-75 percentiles
                std(boolean): The standard deviation of hte posterior rate distribution
                gapped(boolean): MSA DATA, the number of aligned sequences having an amino acid (non-gapped) from the overall
                    number of sequences at each position

        Note: If name is more than one word, then the words should not be separated using a space (For example, name="Silly Potatoes")
        becuase this causes errors in the file handling. Instead, workds should be separated using characters such as an
        underscore (for example, name="Silly_Potatoes").
        """
        if name:
            self.name = name
        elif os.path.isfile(sequence):
            filename = os.path.basename(sequence)
            self.name = filename.split('.')[0]
        else:
            self.name = "Protein_Sequence"
        self.input = sequence
        self.cache = cache
        self.profile = profile
        self.identity = identity
        self.score = score
        self.qqint = qqint
        self.gapped = gapped
        self.std = std
        self.orthologs = ""
        self.alignment = ""
        self.scores = None
        self.alpha = None

    def call_orthologs(self):
        """
        Retrieves the HOGS of the input sequence. This is done by querying the OMA online database.
        """
        if os.path.isfile(self.input):
            with open(self.input, "r") as file:
                sequence = file.read()
            ortholog_call = oma.OrthologFinder(sequence)
        else:
            ortholog_call = oma.OrthologFinder(self.input)
        try:
            self.orthologs = ortholog_call.get_HOGs()
        except RequestException:
            self.orthologs = ortholog_call.get_orthologs()
        with open("%s.orth" % (self.name), "w") as o_file:
            o_file.write(self.orthologs)
        return os.getcwd() + os.sep + "%s.orth" % (self.name)

    def call_alignment(self, orthologs):
        """
        Calls Mafft to generate an MSA of the orthologs that have been input.
        Args:
            orthologs(str): The filepath to the file containing the orthologs of the input, in fasta format
        Returns:
            The filepath to the the msa
        """
        alignment = aminoCons.build_alignment(orthologs)
        self.alignment = alignment
        return alignment

    def call_rate4site(self, msa):
        """
        Calls Rate4Site to calculate various statistics of the amino acids in the input sequence
        Args:
            msa(str): The filepath to the file containing the msa
        Returns:
            The alpha parameter of the data
        """
        conservation_score = aminoCons.Rate4Site(msa, profile=self.profile, cache=self.cache, identity=self.identity,
                                                 score=self.score, qqint=self.qqint, gapped=self.gapped, std=self.std)
        self.scores = conservation_score.run()
        self.alpha = conservation_score.alpha
        conservation_score.close()
        return self.alpha

    def pipe(self):
        """
        Queries the OMA database, Mafft and Rate4Site in sequence to get the
        Returns:
            A dictionary containing the various statistical scores mapped to each amino acid, depending
            on which inputs were selected.
        """
        old_dir = os.getcwd()
        directory = old_dir + os.sep + 'Sequence_Alignments'
        if not os.path.isdir(directory):
            os.makedirs(directory)
        os.chdir(directory)

        msa = directory+os.sep+'%s.aln' % (self.name)
        if os.path.isfile(msa):
            aln = msa
            self.call_rate4site(aln)
        else:
            orth = self.call_orthologs()
            aln = self.call_alignment(orth)
            self.call_rate4site(aln)
            os.remove(os.getcwd() + os.sep + "%s.orth" % (self.name))

        aminoCons.clean_alignment(aln, self.cache)
        os.chdir(old_dir)
        if not self.cache and not os.listdir(directory):
            os.rmdir(directory)
        return self.scores
