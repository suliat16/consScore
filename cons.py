#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intakes a sequence in protein single letter alphabet, and returns the orthologs
of the closest protein to that sequence. Either a string in fasta format, or the
single letter alphabet sequence is required as input.
"""

import json
import requests
import re
import os
from biskit.errors import BiskitError
from requests import exceptions

class SequenceError(BiskitError):
    pass

class OrthologFinder:

    """
    Queries OMA with a protein sequence or fasta to try and retrieve the
    orthologs of that protein. The ortholog data includes the original
    sequence, which is first in the results, as well as the sequence of the
    closest protein match in the database, which is labelled with self. Note
    if the sequence entered is not in fasta format, then the resulting string will
    start with a sequence not in fasta format, followed by orthologs in fasta format
    """
    #Todo: Add something about what HOGS are and what they do

    OMA_BASE_URL = 'https://omabrowser.org'
    HEADERS = {'Content-Type': 'application/json'}

    def __init__(self, fasta):
        self.fasta = fasta
        self.sequence = ""
        self.id = ""
        self.ortholog_ids = []
        self.orthologs = ""
        self.has_run = False
        self.has_run_hogs = False
        self.save_status = 0
        self.hog_level = ""

    def retrieve_OMAid(self):
        """
        Takes a protein sequence and returns the oma id of the best protein
        match
        Args:
            sequence (str): The single letter protein sequence to be queried
        Returns:
           A string containing the ID of the best protein match for the entered sequence
        """
        url = OrthologFinder.build_url(tail='/api/sequence/?query={0}', variation=[self.sequence])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_protID(response)
        if response.status_code == 504:
            self.save_status = response.status_code
            raise TimeoutError('The database timed out. Could not determine the orthologs of your sequence. Status code {0}'.format(self.save_status))
        if response.status_code != 200:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'.format(self.save_status))

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
        url = OrthologFinder.build_url(tail='/api/hog/{0}/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            return self.read_HOGid(response, root)
        if response.status_code == 504:
            self.save_status = response.status_code
            raise TimeoutError ('The database timed out. Could not determine the orthologs of your sequence. Status code {0}'.format(self.save_status))
        if response.status_code != 200:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'.format(self.save_status))

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

        Args:
            omaid(str): An OMA specific protein ID
        Returns:
            A list of strings, the canonical IDS for the orthologs of the protein
        """
        url = OrthologFinder.build_url(tail='/api/protein/{0}/orthologs/', variation=[self.id])
        response = requests.get(url, headers=self.HEADERS)
        if response.status_code == 200:
            self.read_resp_orthoIDs(response)
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'.format(self.save_status))

    def read_resp_orthoIDs(self, response):
        """
        Retrieves the canonical IDs for the orthologs of the protein
        """
        intro = json.loads(response.content.decode('utf-8'))
        for i in intro:
            self.ortholog_ids.append(i['canonicalid'])

    def ortholog_to_fasta(self):
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
        url = OrthologFinder.build_url(tail='/oma/vps/{0}/fasta/', variation=[self.id])
        response = requests.get(url)
        if response.status_code == 200:
            self.orthologs = str(response.text)
            return self.orthologs
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'.format(self.save_status))

    def HOG_to_fasta(self):
        """
        Retrieves the fasta file containing the sequences of the proteins in the HOG of the input protein
        """
        url = OrthologFinder.build_url(tail='/oma/hogs/{0}/{1}/fasta/', variation=[self.id, self.hog_level])
        response = requests.get(url)
        if response.status_code == 200:
            self.HOGs = str(response.text)
            return self.HOGs
        else:
            self.save_status = response.status_code
            raise exceptions.RequestException('There was an issue querying the database. Status code {0}'.format(self.save_status))

    def get_HOGs(self):
        """
        Retrieves a fasta file containing the sequences of the proteins in the HOG to the input protein, based on the input
        parameters.
        """
        if not self.fasta:
            raise SequenceError("Input sequence is empty!")
        self.sequence = OrthologFinder.get_fasta_sequence(fasta=self.fasta)
        if self.has_run_hogs:
            output = self.HOGs
            return output
        else:
            self.retrieve_OMAid()
            self.retrieve_HOG_level()
            output = self.HOG_to_fasta()
            output = OrthologFinder.remove_protein(output, self.id)
            self.has_run_hogs = True
            return output

    def get_orthologs(self):
        """
        Retrieves a fasta file containing the sequences of the orthologous proteins, based on the input parameters
        """
        if not self.fasta:
            raise SequenceError("Input sequence is empty!")
        self.sequence = OrthologFinder.get_fasta_sequence(fasta=self.fasta)
        if self.has_run:
            output = self.orthologs
        else:
            self.retrieve_OMAid()
            output = self.ortholog_to_fasta()
            output = OrthologFinder.remove_first_protein(output)
            output = OrthologFinder.seqnwl_strip(self.sequence) + os.linesep + output
            self.has_run = True
        return output

    @classmethod
    def build_url(cls, tail, variation, base_url=OMA_BASE_URL):
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

    @classmethod
    def get_fasta_sequence(cls, fasta, index=0):
        """
        Given a fasta file, return the sequence at the given index
        Args:
            index(int): For a fasta file with multiple proteins, is the zero
                indexed position of the desired protein within the file. So
                for a file containing 5 proteins, and index of 3 would correspond
                to the 4th protein
        Returns:
            The sequence of the specified protein, as a single string, with newline
            characters removed.
        """
        fstr = OrthologFinder.indv_block(st=fasta)
        fstr = fstr[index]
        fstr = fstr.splitlines()
        for f in fstr:
            if f.startswith('>'):
                fstr.remove(f)
        fstr = "".join(fstr)
        return fstr

    @classmethod
    def indv_block(cls, st):
        """
        Return the header line and the sequence of individual constructs in a file
        Args:
            st(str): The text contained in a fasta file, as a string. Consists of a
                header, which is initiated by > and ends with a newline. Subsequent
                lines are sequence data, until another > is found.

        Returns:
            A list of strings, where each string is the construct header and sequence,
            as a single string. For example, a file containing 4 proteins would
            a list of 4 strings. Each string begins with >, and contains both the
            headers and the newline characters.
        """
        if st.startswith('>'):
            fstr = re.split('>', st)
            seq_list= []
            for f in fstr:
                if f:
                    f = '>' + f
                    f = f.rstrip()
                    seq_list.append(f)
            return seq_list
        else: return [st]

    @classmethod
    def seqnwl_strip(cls, string):
        """
        Removes the newline characters from within the sequences of the fasta
        string
        Args:
            string (str): The fasta sequence with excess newline characters
        Returns:
            A fasta string without the excess newline characters- Retains the
            newline character at the end of the header line
        """
        seqhead = OrthologFinder.header_check(string)
        newlist = seqhead.split(os.linesep)
        newlist = list(filter(None, newlist))
        header = newlist[0] + os.linesep + newlist[1]
        newlist.pop(0)
        newlist.insert(0, header)
        newstring = ''.join(newlist)
        return newstring

    @classmethod
    def header_check(cls, sequences):
        """
        Checks to see if input string has a fasta identification line- if not,
        supplies a default identification line
        Args:
            sequences (str): The sequence(s), in fasta format
        Returns:
            The string, with default identification line, ">Input Sequence\n"
            added if needed
        """
        ##Insert mutated sequence before aligning- make it the reference sequence :)
        if sequences.startswith(">"):
            seq = sequences
        elif not sequences:
            raise SequenceError("Empty Sequence entered.")
        elif sequences[0].isalpha():
            seq = ">Input Sequence" + os.linesep + sequences
        else:
            raise SequenceError("Not a sequence. Please try again")
        return seq

    @classmethod
    def remove_first_protein(cls, fasta):
        """
        Removes the first protein in a string containing proteins in fasta format
        Args:
            fasta(str): The proteins in fasta format
        Returns: A string containing the proteins in fasta format, without the first protein
        """
        fasta_list = OrthologFinder.indv_block(fasta)
        fasta_list.pop(0)
        fasta_string = ''
        for f in fasta_list:
            fasta_string = fasta_string + f
        return fasta_string

    @classmethod
    def remove_protein(cls, fasta, id):
        """
        Removes the protein matching the entered id from the fasta string
        Args:
            fasta(str): The identification line and sequences of all the proteins, in fasta format
            id(str): The OMA id of the closest protein match to the input sequence
        Returns:
             A string containing the proteins in fasta format, without the protein with the entered id
        """
        fasta_list = OrthologFinder.indv_block(fasta)
        for protein in fasta_list:
            if id in protein:
                fasta_list.remove(protein)
        fasta_string = ''
        for protein in fasta_list:
            fasta_string = fasta_string + protein
        return fasta_string
