"""
Intakes a sequence in protein single letter alphabet, and returns the orthologs
of the closest protein to that sequence. Either a string in fasta format, or the
single letter alphabet sequence is required as input.
"""

import json
import os
import requests
from consScore import constool
from biskit.errors import BiskitError
from requests import exceptions


class SequenceError(BiskitError):
    pass


class OrthologFinder:
    """
    Queries OMA with a protein sequence or fasta to try and retrieve the
    orthologs of that protein. If get_HOGs is selected instead, the Hierarchical
    Orthologous Group of the protein will be retrieved instead of the orthologs, which provides
    a smaller group of more closely related orthologs than call_orthologs will return.
    """

    OMA_BASE_URL = 'https://omabrowser.org'
    HEADERS = {'Content-Type': 'application/json'}

    def __init__(self, fasta, rel_type=None):
        self.fasta = fasta
        self.sequence = ""
        self.id = ""
        self.rel_type = rel_type
        self.ortholog_ids = []
        self.orthologs = ""
        self.has_run = False
        self.has_run_hogs = False
        self.save_status = 0
        self.hog_level = ""
        self.HOGs = ""

    def retrieve_OMAid(self):
        """
        Takes a protein sequence and returns the oma id of the best protein
        match
        Returns:
           A string containing the ID of the best protein match for the entered sequence
        """
        url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/api/sequence/?query={0}', variation=[self.sequence])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_protID(response)
        if response.status_code == 504:
            self.save_status = response.status_code
            raise TimeoutError('The database timed out. Could not determine the orthologs of your sequence. Status code {0}'
                               .format(self.save_status))
        if response.status_code != 200:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'
                                              .format(self.save_status))

    def read_resp_protID(self, response):
        response = json.loads(response.content.decode('utf-8'))
        save = response['targets']
        self.id = save[0]['omaid']

    def retrieve_HOG_level(self, root=True):
        """
        Retrieve information on the taxonomic levels that the HOG spans through
        Args:
            root(Boolean): if true, return the deepest level of the HOG. If false, return a list
            of the alternative level that the HOG spans through
        Returns: The deepest level relating the HOG, or a list of all the levels
        """
        url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/api/hog/{0}/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            return self.read_HOGid(response, root)
        if response.status_code == 504:
            self.save_status = response.status_code
            raise TimeoutError('The database timed out. Could not determine the orthologs of your sequence. Status code {0}'
                               .format(self.save_status))
        if response.status_code != 200:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'
                                              .format(self.save_status))

    def read_HOGid(self, response, root):
        """
        If root is true, return the level of the retrieved HOG. If false, return the list of all the
        taxonomic levels that the HOG spans
        """
        response = json.loads(response.content.decode('utf-8'))
        if root:
            level = response[0]['level']
        else:
            level = response[0]['alternative_levels']
        self.hog_level = level
        return self.hog_level

    def update_orthoIDs(self):
        """
        Takes the OMA specific ID of a protein species, and returns a list of the
        canonical IDs of the orthologs of that protein
        Returns:
            A list of strings, the canonical IDS for the orthologs of the protein
        """
        if not self.rel_type:
            url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/api/protein/{0}/orthologs/rel_type={1}', variation=[self.id, self.rel_type])
        else:
            url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/api/protein/{0}/orthologs/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_orthoIDs(response)
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'
                                              .format(self.save_status))

    def read_resp_orthoIDs(self, response):
        """
        Retrieves the canonical IDs for the orthologs of the protein
        """
        intro = json.loads(response.content.decode('utf-8'))
        for i in intro:
            self.ortholog_ids.append(i['canonicalid'])

    def ID_to_fasta(self):
        """
        Takes an OMA specific ID and returns a fasta string of the orthologs assciated
        with that protein
        Returns:
            A single fasta string containing the orthologs of that protein, as
            dictated by OMA.Note that when the fasta file is parsed,
            the first id is the OMA ID, and the second is the canonical id.
        """
        url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/oma/vps/{0}/fasta/', variation=[self.id])
        response = requests.get(url)
        if response.status_code == 200:
            self.orthologs = str(response.text)
            if self.rel_type:
                self.orthologs = self.filter_rel(self.orthologs)
            return self.orthologs
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}')

    def filter_rel(self, orthologs):
        """Removes out all orthologs types except for the one that was specified"""
        ortho_list = constool.indv_block(orthologs)
        new_ortho = []
        for ortholog in ortho_list:
            if self.rel_type in ortholog:
                new_ortho.append(ortholog)
        filtered = ''.join(new_ortho)
        return filtered

    def HOG_to_fasta(self):
        """
        Retrieves the fasta file containing the sequences of the proteins in the HOG of the input protein
        """
        url = constool.build_url(base_url=self.OMA_BASE_URL, tail='/oma/hogs/{0}/{1}/fasta/',
                                 variation=[self.id, self.hog_level])
        response = requests.get(url)
        if response.status_code == 200:
            self.HOGs = str(response.text)
            return self.HOGs
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'
                                              .format(self.save_status))

    def get_HOGs(self):
        """
        Retrieves a fasta file containing the sequences of the proteins in the HOG to the input protein, based on the input
        parameters.
        """
        if not self.fasta:
            raise SequenceError("Input sequence is empty!")
        self.sequence = constool.get_fasta_sequence(fasta=self.fasta)
        if self.has_run_hogs:
            output = self.HOGs
        else:
            self.retrieve_OMAid()
            self.retrieve_HOG_level()
            output = self.HOG_to_fasta()
            output = constool.remove_protein(output, self.id)
            self.has_run_hogs = True
        return output

    def get_orthologs(self):
        """
        Retrieves a fasta file containing the sequences of the orthologous proteins, based on the input parameters
        """
        if not self.fasta:
            raise SequenceError("Input sequence is empty!")
        self.sequence = constool.get_fasta_sequence(fasta=self.fasta)
        if self.has_run:
            output = self.orthologs
        else:
            self.retrieve_OMAid()
            output = self.ID_to_fasta()
            output = constool.remove_first_protein(output)
            output = constool.seqnwl_strip(self.sequence) + os.linesep + output
            self.has_run = True
        return output
