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

import oma
import aminoCons
import argparse
import os

OLD_DIR = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("sequence", help="Input the sequence of the protein of interest")
parser.add_argument("--hogs", action="store_true", help="When specified, the Hierarchical Orthologous Group (HOGs) of the sequence is returned")
parser.add_argument("--name", default="Protein_sequence", help="The name of the protein. All files generated will be based on this name")
parser.add_argument("--identity", action="store_true", help="The identity of the amino acid (Single letter code) at each position")
parser.add_argument("--score", action="store_true", help="The conservation scores. lower value = higher conservation")
parser.add_argument("--qqint", action="store_true", help="QQ-INTERVAL, the confidence interval for the rate estimates. The default interval is 25-75 percentiles")
parser.add_argument("--std", action="store_true", help="The standard deviation of the posterior rate distribution")
parser.add_argument("--gapped", action="store_true", help="MSA DATA, the number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position")
args = parser.parse_args()

with open(args.sequence, 'r') as prot_file:
    prot_seq = prot_file.read()

if args.hogs:
    seq2ortho = oma.OrthologFinder(prot_seq)
    orthologs = seq2ortho.get_HOGs()
else:
    seq2ortho = oma.OrthologFinder(prot_seq)
    orthologs = seq2ortho.get_orthologs()

with open("%s.txt" %(args.name), 'w') as seq_file:
    seq_file.write(orthologs)

alignment = aminoCons.build_alignment(os.getcwd() + os.sep + "%s.txt"%(args.name))

cons_dict = aminoCons.Rate4Site(msa= (alignment + os.sep+ "%s.aln"%(args.name)), identity=args.identity,
                                  score=args.score, qqint=args.qqint,std=args.std, gapped=args.gapped)
cons_dict = cons_dict.run()


with open("%s data"%(args.name), 'w') as mat_file:
    output = str(cons_dict)
    mat_file.write(output)
