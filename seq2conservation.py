#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline which retrieves the conservation score of the amino acids in sequence
given an input sequence or fasta file.

The steps of the pipeline are as follows:
Sequence --> Orthologs
Orthologs --> Alignment
Alignment --> Conservation scores
"""

import cons
import aminoCons
import os
from biskit.errors import BiskitError

class PipelineError(BiskitError):
    pass

class ConservationPipe():

    """
    """

    def __init__(self, sequence, name=None, cwdir=None, cache=True, identity=True, score=True, qqint=False, std=False,
                    gapped=False):
        """
        """
        if name:
            self.name = name
        if os.path.isfile(sequence):
            filename = os.path.basename(sequence)
            self.name = filename.split('.')[0]
        else:
            self.name = "Protein Sequence"
        self.input = sequence
        self.old_path = os.getcwd()

        self.cwdir = cwdir
        self.cache = cache
        self.identity = identity
        self.score = score
        self.qqint = qqint
        self.gapped = gapped
        self.std = std


    def call_orthologs(self):
        """
        """
        #TODO: exception handling in case of no hogs- call normal oma
        ortholog_call = cons.OrthologFinder(self.input)
        orthologs = ortholog_call.get_HOGs()
        self.orthologs = orthologs
        with ("%s orthologs" %(self.name), 'w') as o_file:
            o_file.write(orthologs.get_HOGs())
        return os.getcwd() + os.sep + "%s orthologs"%(self.name)


    def call_alignment(self, orthologs):
        """
        """
        #TODO: to be called only if there are no pregenerated msa- they take too long
        alignment = aminoCons.build_alignment(orthologs)
        self.alignment = alignment
        return alignment

    def call_rate4site(self, msa):
        """
        """
        #TODO: Ask about close() method- how to design so that I can call it below
        conservation_score = aminoCons.Rate4Site(msa, cwdir=self.cwdir, cache=self.cache, identity=self.identity,
                                                 score=self.score, qqint=self.qqint, gapped=self.gapped, std= self.std)
        self.cons = conservation_score
        self.scores = conservation_score.run()
        return conservation_score.alpha

    def cleanup(self, msa):
        """
        """
        self.cons.close()
        aminoCons.clean_alignment(msa)
        os.chdir(self.old_path)
