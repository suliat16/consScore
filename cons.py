#/usr/bin/env python3

"""
Intakes a sequence in protein single letter alphabet, and returns the orthologs
of the closest protein to that sequence. Either a string in fasta format, or the
single letter alphabet sequence is required as input. 
"""

import json
import requests
import re
import os
from lxml import html

class OrthologFinder:

    """
    Queries OMA with a protein sequence or fasta to try and retrieve the
    orthologs of that protein. The ortholog data includes the original
    sequence, which is first in the results, as well as the sequence of the 
    closest protein match in the database, which is labelled with self. Note
    if the sequence entered is not in fasta format, then the resulting string will
    start with a sequence not in fasta format, followed by orthologs in fasta format
    """

    OMA_BASE_URL = 'https://omabrowser.org'
    HEADERS = {'Content-Type': 'application/json'}


    def __init__(self, fasta):
        self.fasta = fasta
        self.sequence = ""
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
        for i in intro:
            self.ortholog_ids.append(i['canonicalid'])

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
        Takes the ortholog response, and returns the fasta string
        """
        self.orthologs = str(response.text)
        return self.orthologs

    def get_orthologs(self):
        """
        Returns the orthologous proteins to the sequence stored in the object
        """
        ##Can I put the return statements in a finally block?
        try:
            output = None
            self.sequence =  self.get_fasta_sequence(fasta=self.fasta)
            if self.has_run:
                output = self.orthologs
                return output
            else:
                self.has_run = True
                self.retrieve_OMAid()
                output = self.OMA_to_fasta()
                output = self.fasta + output
                return output
        except ImportError:
            output = 'There was an issue querying the database. Status code {0}'.format(self.save_status)
            return output
        except TimeoutError:
            output = 'The database timed out. Could not determine the orthologs of your sequence. Status code {0}'.format(self.save_status)
            return output
        except IndexError: 
            output = "Input sequence is empty!"
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
    
    def get_fasta_sequence(self, fasta):
        """
        Given a fasta file, return the sequence at the given index
        
        Args:
            index(int): For a fasta file with multiple proteins, is the zero
                indexed position of the desired protein within the file
        
        Returns: 
            The sequence of the specified protein, as a single string, with newline
            characters removed.
        """
        fstr= self.indv_block(st=fasta)
        fstr=fstr[0]
        fstr= fstr.splitlines()
        for f in fstr:
            if f.startswith('>'):
                fstr.remove(f)
        fstr= "".join(fstr)
        return fstr
    
    def indv_block(self, st=""):
        """
        Return the header line and the sequence of individual constructs in a file
        
        Args:
            st(str): The text contained in a fasta file, as a string. Consists of a 
                header, which is inititated by > and ends with a newline. Subsequent 
                lines are sequence data, until another > is found.
            
        Returns: 
            A list of strings, where each string is the construct header and sequence,
            as a single string. For example, a file containing 4 proteins would
            a list of 4 strings. Each string begins with >, and contains both the
            headers and the newline characters. 
        """
        fstr=re.split('>',st)
        fstr = list(filter(None, fstr))
        fstr = ['>' + f for f in fstr]
        return fstr
    
    def seqnwl_strip(self, string):
        """
        Removes the newline characters from within the sequences of the fasta 
        string
        Args:
            string (str): The fasta sequence with excess newline characters
        Returns: 
            A fasta string without the excess newline characters- Retains the 
            newline character at the end of the header line
        """
        seqlist = self.indv_block(string)
        fasta = []
        for seq in seqlist:
            newlist = seq.split(os.linesep)
            header = newlist[0] + os.linesep + newlist [1]
            newlist.pop(0)
            newlist.pop(1)
            newlist.insert(0, header)
            newstring = ''.join(newlist)
            fasta.append(newstring)
        fasta = os.linesep.join(fasta)
        return fasta
    
    

