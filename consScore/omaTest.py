#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:49:39 2018

@author: suliat16

Test file for oma
"""
import biskit.test
import os
from consScore import oma
from consScore import constool
from requests import exceptions
from unittest.mock import patch, MagicMock


class TestCons(biskit.test.BiskitTest):
    """
    Test suite testing the behaviour of the oma module
    """

    TAGS = [biskit.test.NORMAL]

    def setUp(self):
        lysozyme = ("MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGY"
                    "NTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVAC"
                    "AKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
        self.lyz = oma.OrthologFinder(lysozyme)
        stragg = ('YYYYYYYYYYYVVVVVVVVVVVVBBBBBBBBBBBAAAAAAAAAAANNNNNNNNNNN'
                  'MMMMMMMMMMMWWWWWWWWWWWWWCCCCCCCCCCSSSSSSSSSS')
        self.aggregate = oma.OrthologFinder(stragg)
    @classmethod
    def setUpClass(cls):

        with open((os.getcwd()+os.sep+'example_data'+os.sep+'CDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = oma.OrthologFinder(arabidopsisCDC48A)

        with open((os.getcwd()+os.sep+'example_data'+os.sep+'OMA_get_id.txt'), 'r') as file:
            OMA_id_response = file.read()
        cls.response = bytes(OMA_id_response, 'utf-8')

        with open((os.getcwd() + os.sep + 'example_data' + os.sep + 'orhtolog_response.txt'), 'r') as file:
            ortho_response = file.read()
        cls.oresponse = bytes(ortho_response, 'utf-8')

        with open((os.getcwd() + os.sep + 'example_data' + os.sep + 'fasta_response.txt'), 'r') as file:
            fasta_response = file.read()
        cls.fresponse = bytes(fasta_response, 'utf-8')

        with open((os.getcwd() + os.sep + 'example_data' + os.sep + 'ATN1_HOGS.txt'), 'r') as file:
            HOG_response = file.read()
        cls.hresponse = bytes(HOG_response, 'utf-8')

        with open((os.getcwd() + os.sep + 'example_data' + os.sep + 'hog_content.txt'), 'r') as file:
            level_response = file.read()
        cls.lvlresponse = bytes(level_response, 'utf-8')

    def test_read_resp_retOMA(self):
        """tests that read_resp_protID correctly parses the JSON response and retrieves the protein id"""
        the_response = MagicMock(content = self.response)
        self.CDC48A.read_resp_protID(the_response)
        self.assertTrue(self.CDC48A.id == "ARATH09528")
        
    def test_read_resp_upOMA(self):
        """tests that read_resp_protID correctly parses the JSON response and saves a list of the ortholog ids"""
        the_response = MagicMock(content=self.oresponse)
        self.CDC48A.read_resp_orthoIDs(the_response)
        self.assertTrue("H3DFZ9" in self.CDC48A.ortholog_ids)

    def test_build_url_retOMA(self):
        """Tests that build URL inserts the variable portions of the url correcly into the base"""
        self.lyz.sequence = 'MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV'
        self.assertEqual(constool.build_url(base_url='https://omabrowser.org', tail='/api/sequence/?query={0}', variation=[self.lyz.sequence]), 'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')

    @patch('oma.requests.get')
    def test_ologs_stragg(self, requests_mock):
        """Tests that a request with a bad status code raises an exception with call_orthologs"""
        requests_mock.requests.get.return_value = None
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException) as cm:
            straggFasta = self.aggregate.get_orthologs()
        err = cm.exception
        self.assertTrue('There was an issue querying the database. Status code' in str(err))

    @patch('oma.requests.get')
    def test_ofasta_stragg(self, requests_mock):
        """Tests that a bad request raises an exception in ortholog_to_fasta"""
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.ortholog_to_fasta()

    @patch('oma.requests.get')
    def test_ofasta_cdc(self, requests_mock):
        """Tests that ortholog_to_fasta correctly parses the response data"""
        requests_mock().status_code = 200
        requests_mock().text = self.fresponse
        test = self.CDC48A.ortholog_to_fasta()
        self.assertTrue('[Arabis alpina]' in test)

    @patch('oma.requests.get')
    def test_orIDs_stragg(self, requests_mock):
        """Test that update_orthoIDs returns an exception given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.update_orthoIDs()

    @patch('oma.requests.get')
    def test_ret_hogs_stragg(self, requests_mock):
        """Tests that retrieve_HOG_level throws an exception when given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.retrieve_HOG_level()

    def test_empty_input(self):
        """Checks that the correct exception is raised when an empty sequence is entered"""
        fs = oma.OrthologFinder("")
        with self.assertRaises(oma.SequenceError) as cm:
            fs.get_orthologs()
        err = cm.exception
        self.assertEqual(str(err), 'Input sequence is empty!')

    def test_empty_input_HOG(self):
        """Checks that the correct exception is raised when an empty sequence is entered"""
        fs = oma.OrthologFinder("")
        with self.assertRaises(oma.SequenceError) as cm:
            fs.get_HOGs()
        err = cm.exception
        self.assertEqual(str(err), 'Input sequence is empty!')

    @patch('oma.requests.get')
    def test_hog_fasta(self, mock_request):
        """Tests that HOG_to_fasta correctly parses the request response"""
        self.lyz.hog_level = "Amniota"
        mock_request().status_code = 200
        mock_request().text = self.hresponse
        test = self.lyz.HOG_to_fasta()
        self.assertTrue("HOG:0377891.2a.2a" in test)

    @patch('oma.requests.get')
    def test_read_hog_roottrue(self, mock_request):
        """tests that read_HOGid retrieves the root ID when root=True"""
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=True)
        self.assertEqual(test, 'Amniota')

    @patch('oma.requests.get')
    def test_read_hog_rootfalse(self, mock_request):
        """tests that read_HOGid retrieves a list of the ids when root=False """
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=False)
        self.assertTrue('Hystricomorpha' in test)
        self.assertTrue('Caniformia' in test)
        self.assertTrue('Gorilla gorilla gorilla' in test)

if __name__ == '__main__':
    biskit.test.localTest()
