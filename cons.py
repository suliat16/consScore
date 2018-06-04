#/usr/bin/env python3

"""
    Intakes a protein sequence in fasta format, and returns the conservation
    score of the residues in the protein
"""

#import os
import json
import requests
from Bio import SeqIO as sq
from Bio.Align.Applications import ClustalwCommandline as cline
from Bio import AlignIO as al
from Bio.Align import AlignInfo as ai
from Bio.SubsMat import FreqTable

ORTHODB_BASE_URL = 'http://www.orthodb.org/'
OMA_BASE_URL = 'http://omabrowser.org/api'
EBI_BASE_URL = 'www.ebi.ac.uk/proteins/api'
HEADERS = {'Content-Type': 'application/json'}
FREQUENCIES = dict(A=0.05, C=0.05, D=0.05, E=0.05, F=0.05, G=0.05, H=0.05,
                   I=0.05, K=0.05, L=0.05, M=0.05, N=0.05, P=0.05, Q=0.05, 
                   R=0.05, S=0.05, T=0.05, V=0.05, W=0.05, Y=0.05)
FREQ_TABLE = FreqTable.FreqTable(FREQUENCIES, 2)

def retrieve_OMA(sequence):
    """
    Takes a protein sequence and returns a json object with ortholog
    information
    
    Args:
        sequence (str): The single letter protein sequence to be queried
        
    Returns: 
        A parseable json object containing multiple dictionaries, each 
        corresponding to an ortholog
    """
    url = '{0}/sequence/?query={1}'.format(OMA_BASE_URL, sequence)
    response = requests.get(url, headers=HEADERS)
    
    if response.status_code == 200:
        response = json.loads(response.content.decode('utf-8'))
        #TODO: Save the next parameter as an object variable
        save = response['targets']
        return save[0]['canonicalid']
        
    else:
        if response.status_code == 500:
            print('[!][{0}] Server Error'.format(response.status_code))
        elif response.status_code == 404:
            print('[!] [{0}] URL not found: [{1}]'.format(response.status_code, url))
        elif response.status_code == 401:
            print('[!] [{0}] Authentication Failed'.format(response.status_code))
        elif response.status_code == 400:
            print('[!] [{0}] Bad Request'.format(response.status_code))
        elif response.status_code == 300:
            print('[!] [{0}] Unexpected Redirect'.format(response.status_code))
        else:
            print('[?] Unexpected Error: [HTTP {0}]: Content: {1}'.format(response.status_code, response.content))
        return None

def OMA_to_orthoID(omaid):
    """
    """
    url = '{0}/protein/{1}/orthologs/'.format(OMA_BASE_URL, omaid)
    response = requests.get(url, headers=HEADERS)
    if response. status_code == 200:
        intro = json.loads(response.content.decode('utf-8'))
        return intro['canonicalid']
    else:
        if response.status_code == 500:
            print('[!][{0}] Server Error'.format(response.status_code))
        elif response.status_code == 404:
            print('[!] [{0}] URL not found: [{1}]'.format(response.status_code, url))
        elif response.status_code == 401:
            print('[!] [{0}] Authentication Failed'.format(response.status_code))
        elif response.status_code == 400:
            print('[!] [{0}] Bad Request'.format(response.status_code))
        elif response.status_code == 300:
            print('[!] [{0}] Unexpected Redirect'.format(response.status_code))
        else:
            print('[?] Unexpected Error: [HTTP {0}]: Content: {1}'.format(response.status_code, response.content))
        return None

def OMA_to_seqdict(OMAlist, sequencelist):
    """
    Takes a list of Canonical protein IDs and a list of protein sequences and
    returns a dictionary mapping the IDs to the sequences.
    
    Args:
        OMAlist(list): A list of strings-canonical protein IDs, which can be used
        to query Uniprot, for example
        sequencelist(list): A list of strings, where the strings are the single
        letter protein sequences of the proteins referenced in the OMAlist. Note
        that the entries of this list must be in the same order as the entries in 
        the previous list.
        
    Returns:
        A dictionary, where each key is a member of OMAlist, and maps to a value
        from sequencelist, at the same index.
    """ 
    seq_dict = {}
    for idx, o in enumerate(OMAlist):
        seq_dict[o] = sequencelist[idx]
    return seq_dict
      
def get_orthoDBids(filepath):
    """
    Take a fasta file and return a list of OrthoDB cluster IDs. Note that this
    is only useful for animal proteins.

    Args:
        filepath (string): The absolute path leading to and including the fasta
        file of interest.
    Returns:
        A list of OrthoDB ortholog IDs which can then be used to retrieve a list
        of orthologs in a fasta file.
    """
    with open(filepath) as workfile:
        prerec = list(sq.parse(workfile, "fasta"))
    rec = prerec[0]
    st = rec.seq
    url = '{0}/blast?seq={1}&level=33208&limit=100'.format(ORTHODB_BASE_URL, st)
    # level in this case refers to the taxid from ncbi database- currently set
    # at metazoa
    # limit 100 means that at most 100 OrthoDB IDs will be returned- not general
    response = requests.get(url, headers=HEADERS)
    if response.status_code == 200:
        response = json.loads(response.content.decode('utf-8'))
        return response['data']
        
    # What exception should I raise? Do I have to make my own?
    else:
        if response.status_code == 500:
            print('[!][{0}] Server Error'.format(response.status_code))
        elif response.status_code == 404:
            print('[!] [{0}] URL not found: [{1}]'.format(response.status_code, url))
        elif response.status_code == 401:
            print('[!] [{0}] Authentication Failed'.format(response.status_code))
        elif response.status_code == 400:
            print('[!] [{0}] Bad Request'.format(response.status_code))
        elif response.status_code == 300:
            print('[!] [{0}] Unexpected Redirect'.format(response.status_code))
        else:
            print('[?] Unexpected Error: [HTTP {0}]: Content: {1}'.format(response.status_code, response.content))
        return None

def get_orthoDBlogs(ID, filename='default'):
    """
    Takes a list of OrthoDB ID strings, and returns a file containing sequences in fasta
    format

    Args: 
        ID(str): An OrthoDB cluster ID. Querying the API generally returns a list,
            so to call this method, feed the function the IDList[index]
        filename (str): The name that will be given to the fasta file that is 
            written.
        
    Returns: A fasta file containing the relevant orthologs in the same directory
        as the program
    """
    url = '{0}/fasta?id={1}'.format(ORTHODB_BASE_URL, ID)
    response = requests.get(url, headers=HEADERS)
    
    if response.status_code == 200:
        fast = response.content.decode('utf-8')
        with open(filename+'.fasta', 'w') as workfile:
            workfile.write(fast)
    #TODO: Call alignment software of choice within the with, so that the temp
            #file can be safely deleted        
            
    elif response.status_code == 500:
        print('[!][{0}] Server Error'.format(response.status_code))
    elif response.status_code == 404:
        print('[!] [{0}] URL not found: [{1}]'.format(response.status_code, url))
    elif response.status_code == 401:
        print('[!] [{0}] Authentication Failed'.format(response.status_code))
    elif response.status_code == 400:
        print('[!] [{0}] Bad Request'.format(response.status_code))
    elif response.status_code == 300:
        print('[!] [{0}] Unexpected Redirect'.format(response.status_code))
    else:
        print('[?] Unexpected Error: [HTTP {0}]: Content: {1}'.format(response.status_code, response.content))

def get_alignment(fastafile, path):
    """
    r"/weyr/software/clustalw2/v2.1-bin.app/bin/clustalw2"

    Takes a fasta file with multiple sequences as input, and outputs an alignment
    file generated by Clustalw2

    Args: 
        fastafile (str): The absolute path to the fasta file being aligned 
        path (str): The absolute path to the clustalw2 app
    
    Returns: A .aln and .dnd file, both showing the alignments of the input
        sequences
    """
    
    clus_cline = cline(path, infile=fastafile)
    #assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
  #  stdout, stderr = clus_cline()  # This is the way its done on the biopython website
    return clus_cline()

def get_information(filepath, fileformat="clustal"):
    """
    """
    prealign = list(al.parse(filepath, fileformat))
    align = prealign[0]
    summary = ai.SummaryInfo(align)
    return summary._get_column_info_content(obs_freq=0.05, e_freq_table=FREQ_TABLE, log_base=2, random_expected = 0)
    
    #return summary.information_content(start=89, end=89, e_freq_table=FREQ_TABLE, chars_to_ignore=['X'])
    #info_content = summary.information_content()
    #return info_content
    
    
    
    
