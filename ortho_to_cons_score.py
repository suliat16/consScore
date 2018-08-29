#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intakes the sequence or fasta string of a protein using the standard, single letter alphabet,
gets the orthologs from an external (online) database, the OMA browser. It then generates a
multiple sequence alignment (MSA) using T-Coffee, and calculates the conservation score of each
amino acid at each position using Rate4Site, with respect to the entered sequence.

Citations
OMA database:
Altenhoff A et al., The OMA orthology database in 2018: retrieving evolutionary relationships among
all domains of life through richer web and programmatic interfaces Nucleic Acids Research, 2018,
46 (D1): D477-D485 (doi:10.1093/nar/gkx1019).

T-Coffee:
T-Coffee: A novel method for multiple sequence alignments.
Notredame,Higgins,Heringa,JMB,302(205-217)2000

Rate4Site:
Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference methods:
Bayesian methods are superior. Mol Biol Evol 21: 1781-1791.
"""

import cons
import aminoCons
import argparse
import os

#Convert command line arguments to variables
parser = argparse.ArgumentParser()
parser.add_argument("sequence", help="Input the sequence of the protein of interest")
parser.add_argument("--hogs", action="store_true", help="When specified, the Hierarchical Orthologous Group (HOGs) of the sequence is returned")
parser.add_argument("--name", default="Protein_sequence", help="The name of the protein. All files generated will be based on this name")
args = parser.parse_args()

with open(args.sequence, 'r') as prot_file:
    prot_seq = prot_file.read()

if args.hogs:
    seq2ortho = cons.OrthologFinder(prot_seq)
    orthologs = seq2ortho.get_HOGs()
else:
    seq2ortho = cons.OrthologFinder(prot_seq)
    orthologs = seq2ortho.get_orthologs()

print(type(orthologs))
with open("%s.txt" %(args.name), 'w') as seq_file:
    seq_file.write(orthologs)

print(os.getcwd() + os.sep + "%s.txt"%(args.name))

alignment = aminoCons.build_alignment(os.getcwd() + os.sep + "%s.txt"%(args.name))

cons_matrix = aminoCons.Rate4Site(alignment)
