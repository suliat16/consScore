#/usr/bin/env python3

"""
Intakes a protein sequence in fasta format, and returns the conservation
score of the residues in the protein
"""

#import os
import json
import os
import requests

class OrthologFinder:

    """
    Queries OMA with a protein sequence to try and retrieve the
    orthologs of that protein. The ortholog data does not include the original
    sequence, but instead the sequence of the closest protein match in the database.
    """
    OMA_BASE_URL = 'https://omabrowser.org'
    HEADERS = {'Content-Type': 'application/json'}

    def __init__(self, sequence):
        self.sequence = sequence
        self.id = None
        self.ortholog_ids = []
        self.orthologs = ""
        self.has_run = False
        self.save_status = 0

    def retrieve_OMAid(self):
        """
        Takes a protein sequence and returns the oma id of the best protein
        match

        Args:
            sequence (str): The single letter protein sequence to be queried

        Returns:
           A string containing the ID of the best protein match for the entered sequence
        """
        url = self.build_url(tail='/sequence/?query={0}', variation=[self.sequence])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            response = json.loads(response.content.decode('utf-8'))
            save = response['targets']
            self.id = save[0]['omaid']

        else:
            self.save_status = response.status_code
            raise ImportError

    def update_OMA_orthoIDs(self):
        """
        Takes the OMA specific ID of a protein species, and returns a list of the
        canonical IDs of the orthologs of that protein

        Args:
            omaid(str): An OMA specific protein ID
        Returns:
            A list of strings, the canonical IDS for the orthologs of the protein
        """

        url = self.build_url(tail='/protein/{0}/orthologs/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            intro = json.loads(response.content.decode('utf-8'))
            self.ortholog_ids = intro['canonicalid']

        else:
            self.save_status = response.status_code
            raise ImportError

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
        url = self.build_url(tail='/oma/vps/{0}/fasta/', variation=[self.id])
        response = requests.get(url)

        if response.status_code == 200:
            orthologs = str(response.text)
            self.orthologs = orthologs.replace(os.linesep, '')
            return self.orthologs
        else:
            self.save_status = response.status_code
            raise ImportError

    def get_orthologs(self):
        """
        Returns the orthologous proteins to the sequence stored in the object
        """
        #I had return statements in every conditional block, but it was giving
        #me issues- why? Also IDK how to beef up that docstring

        try:
            output = None
            if self.has_run:
                output = self.orthologs
            else:
                self.has_run = True
                self.retrieve_OMAid()
                output = self.OMA_to_fasta()
            return output
        except ImportError:
            print('There was an issue querying the database.')
            print('Status code {0}'.format(self.save_status))

    def build_url(self, tail, variation, base_url=OMA_BASE_URL):
        """
        Takes the passed parameters and builds a URL to query the OMA database

        Args:
            tail(str): The path and REST parameters that returns the desired info
            variation(list): A list of strings that contain the parameters unique
                to the query
            base_url(str): The website that is being accessed, without any slashes
        """
        url = base_url + tail
        url = url.format(*variation)
        return url
