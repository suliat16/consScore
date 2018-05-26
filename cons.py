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

base_url = 'http://www.orthodb.org/'
headers = {'Content-Type': 'application/json'}
    
def get_cluster(file = ""):
    """
    Take a fasta file and return a list of OrthoDB cluster IDs. Note that this 
    is only useful for animal proteins.
    
    Args:
        file (string): The name of the fasta file. Note that the fasta file must
            be present in the same directory as the script
    
    Returns: 
        A list of OrthoDB ortholog IDs which can then be used to retrieve a list 
        of orthologs in a fasta file.
    """
    prerec = list(sq.parse(file,"fasta"))
    rec= prerec[0]
    st = rec.seq
    url = '{0}/blast?seq={1}&level=33208&limit=100'.format(base_url, st)
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        response = json.loads(response.content.decode('utf-8'))
        return response['data']

    else: 
        return None 
        # TODO: Wrap this in a try/catch block 
    
def get_orthologs(IDlist= [], index=0, filename= 'default'):
    """
    Takes an OrthoDB ID, and returns a string consisting of sequences in fasta
    format
    
    Args: 
        IDlist (list): A list of strings, where each string is a OrthoDB cluster id
        index(int): Specifies the index of the desired OrthoDB id, to retrieve the 
            orthologs associated with that ID. Defaults to the first string in the
            because the API returns a list with the best matching cluster first
        filename (str): The name that will be given to the fasta file that is 
            written.
        
    Returns: A fasta file containing the relevant orthologs in the same directory
        as the program
        
    """
    ID = IDlist[index]
    url = '{0}/fasta?id={1}'.format(base_url, ID)
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        fast= response.content.decode('utf-8')
        with open(filename+'.fasta', 'w') as workfile:
            workfile.write(fast)
    else:
        return response.status_code
    

    
    
    
    
    
    
