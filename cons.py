#/usr/bin/env python3

"""
Intakes a sequence in protein single letter alphabet, and returns the orthologs
of the closest protein to that sequence
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
        self.id = ""
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
        url = self.build_url(tail='/api/sequence/?query={0}', variation=[self.sequence])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_retOMA(response)
        if response.status_code == 504:
            self.save_status = response.status_code
            raise TimeoutError
        if response.status_code != 200:
            self.save_status = response.status_code
            raise ImportError("Status code:{0}".format(self.save_status))

    def read_resp_retOMA(self, response):
        response = json.loads(response.content.decode('utf-8'))
        save = response['targets']
        self.id = save[0]['omaid']

    def update_OMA_orthoIDs(self):
        """
        Takes the OMA specific ID of a protein species, and returns a list of the
        canonical IDs of the orthologs of that protein

        Args:
            omaid(str): An OMA specific protein ID
        Returns:
            A list of strings, the canonical IDS for the orthologs of the protein
        """
        url = self.build_url(tail='/api/protein/{0}/orthologs/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_upOMA(response)
        else:
            self.save_status = response.status_code
            raise ImportError("Status code:{0}".format(self.save_status))

    def read_resp_upOMA(self, response):
        """
        Retrieves the canonical IDs for the orthologs of the protein
        """
        intro = json.loads(response.content.decode('utf-8'))
        self.ortholog_ids = intro[0]['canonicalid']

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
            return self.read_resp_OMAfasta(response)
        else:
            self.save_status = response.status_code
            message = "Status code: {0}".format(self.save_status)
            raise ImportError(message)

    def read_resp_OMAfasta(self, response):
        """
        Takes the  ortholog response, converts it to a string and removes newlines
        """
        orthologs = str(response.text)
        self.orthologs = orthologs.replace(os.linesep, '')
        return self.orthologs

    def get_orthologs(self):
        """
        Returns the orthologous proteins to the sequence stored in the object
        """
        ##Can I put the return statements in a finally block?
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
            output = 'There was an issue querying the database. Status code {0}'.format(self.save_status)
            return output
        except TimeoutError:
            output = 'The database timed out. Could not determine the orthologs of your sequence. Status code {0}'.format(self.save_status)
            return output

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
