#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

from Bio.Blast.Applications import NcbiblastpCommandline


class SeqSearch:

    """
    Initializing this class creates an object which takes as an input a protein sequences and produces a protein identity
    as output. This class in particular then converts that protein ID into that of a database of choice. Currently
    supported databases are:
        OMA- The OMA orthology database
    """

    def __int__(self, sequence):
        self.sequence = sequence
        self.blast = None

    def multi_prot(self, seq_list):
        """Concatenates the queries so that blast can search for multiple queries at once.
        Useful if searching multiple proteins, as it only goes over the databases once
        Args:
            seq_list(list): A list of the query sequences to be searched using blast
        Returns:
            (A file or a string??) containing the concatenated queries to be searched using blast
        """


    def identify_protein(self, sequence, evalue, db='refseq_protein'):
        """Use blast to determine tne identity of an amino acid sequence"""
        cline = NcbiblastpCommandline(query=sequence, db=db, evalue= evalue)
        self.blast = cline()

