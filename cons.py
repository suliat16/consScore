#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 10:50:45 2018

@author: suliat16

Hacky lil (script?) to get the conservation score of the residues in a protein

Step 1: Open fasta file, and get the identifier- use this to make a call to the
OrthoDB API. This should return a fasta file with the alignments back

Step 2: Enter this fasta file into the ClustalX wrapper inbuilt into Biopython
- return an alignment object 

Step 3: Get the information content for every position in the sequence, and input
it in an array, with index corresponding with aa posn

Step 4: Talk to a bioinformatician- convert the information content of position
from the position
"""

from Bio import SeqIO as sq
import json
import requests

filename = "lysozyme.fasta"
base_url = 'http://www.orthodb.org/'
headers = {'Content-Type': 'application/json'}
    
def get_response(file = ""):
    """
    """
    prerec = list(sq.parse(file,"fasta"))
    rec= prerec[0]
    st = rec.seq
    url = '{0}/blast?seq={1}&level=33208&limit=100'.format(base_url, st)
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        return json.loads(response.content.decode('utf-8'))

    else: 
        return None 
        # TODO: Wrap this in a try/catch block 
    


