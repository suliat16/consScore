#/usr/bin/env python3

"""
    Intakes a protein sequence in fasta format, and returns the conservation
    score of the residues in the protein
"""

#import os
import json
import os
import requests
from Bio.Align.Applications import ClustalwCommandline as cline

class OrthologFinder:
    """
    """
    ORTHODB_BASE_URL = 'http://www.orthodb.org/'
    OMA_BASE_URL = 'http://omabrowser.org/api'
    HEADERS = {'Content-Type': 'application/json'}

    def __init__(self, sequence):
        self.sequence = sequence
        self.id = None
        self.ortholog_ids = []
        self.orthologs = ""
        self.has_run = False
        self.is_OMA = False
        self.is_ortho = False

    def retrieve_OMA(self):
        """
        Takes a protein sequence and returns the oma id of the best protein
        match

        Args:
            sequence (str): The single letter protein sequence to be queried

        Returns:
           A string containing the ID of the best protein match for the entered sequence
        """

        url = '{0}/sequence/?query={1}'.format(self.OMA_BASE_URL, self.sequence)
        response = requests.get(url, headers=self.HEADERS)

        if response.status_code == 200:
            response = json.loads(response.content.decode('utf-8'))
            save = response['targets']
            self.id = save[0]['omaid']
            self.is_OMA = True

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


    def update_OMA_orthoIDs(self):
        """
        Takes the OMA specific ID of a protein species, and returns a list of the
        canonical IDs of the orthologs of that protein

        Args:
            omaid(str): An OMA specific protein ID
        Returns:
            A list of strings, the canonical IDS for the orthologs of the protein
        """
        url = '{0}/protein/{1}/orthologs/'.format(self.OMA_BASE_URL, self.id)
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            intro = json.loads(response.content.decode('utf-8'))
            self.ortholog_ids = intro['canonicalid']

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


    def OMA_to_fasta(self):
        """
        Takes an OMA specific ID and returns a fasta string of the orthologs assciated
        with that protein

        Args:
            OMAid(str):  The OMA specific ID of a desired protein
        Returns:
            A single fasta string containing the orthologs of that protein, as
            dictated by OMA.Note that when the fasta file is parsed,
            the first id is the OMA ID, and the second is the canonical id.
        """
        url = 'https://omabrowser.org/oma/vps/{0}/fasta/'.format(self.id)
        response = requests.get(url)
        orthologs = str(response.text)
        self.orthologs = orthologs.replace(os.linesep, '')
        return self.orthologs

    def get_orthoDBids(self):
        """
        Take a fasta file and return an OrthoDB cluster ID. Note that this
        is only useful for animal proteins.

        Returns:
            An OrthoDB ortholog ID which can then be used to retrieve a list
            of orthologs in a fasta file. The website returns a list of IDs, with
            the best listed first- the first ID of the list is pulled.
        """

        url = '{0}/blast?seq={1}&level=33208&limit=100'.format(self.ORTHODB_BASE_URL, self.sequence)
        # level in this case refers to the taxid from ncbi database- currently set
        # at metazoa
        # limit 100 means that at most 100 OrthoDB IDs will be returned- not general
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            response = json.loads(response.content.decode('utf-8'))
            IDlist = response['data']
            self.is_ortho = True
            self.id = IDlist[0]

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


    def retrieve_orthoDBlogs(self):
        """
        Takes an OrthoDB ID from a list of IDs, and returns a file containing
        sequences in fasta format

        Returns: A string in fasta format containing all of the orthologs for
            that ID
        """
        url = '{0}/fasta?id={1}'.format(self.ORTHODB_BASE_URL, self.id)
        response = requests.get(url, headers=self.HEADERS)

        if response.status_code == 200:
            self.orthologs = response.content.decode('utf-8')
            return self.orthologs
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
    def get_orthologs(self):
        """
        """
        #I had return statements in every conditional block, but it was giving
        #me issues- why?
        output = None
        if self.has_run:
            output = self.orthologs
        else:
            self.has_run = True
            self.retrieve_OMA()
            if self.is_OMA:
                output = self.OMA_to_fasta()
            elif self.is_ortho:
                self.get_orthoDBids()
                output = self.retrieve_orthoDBlogs()
            else:
                output = "Could not determine the orthologs of your sequence."
        return output



#TODO: Move this into another class, jeez
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
