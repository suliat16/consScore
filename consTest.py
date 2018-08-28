#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:49:39 2018

@author: suliat16

Test file for cons
"""
import unittest
import os
import cons
from requests import exceptions
from requests.models import Response
from unittest.mock import patch, MagicMock
import numpy


class TestCons(unittest.TestCase):
    """
    """
    @classmethod
    def setUpClass(cls):
        lysozyme = ("MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGY"
                  "NTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVAC"
                  "AKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
        cls.lyz = cons.OrthologFinder(lysozyme)
        stragg = ('YYYYYYYYYYYVVVVVVVVVVVVBBBBBBBBBBBAAAAAAAAAAANNNNNNNNNNN'
                  'MMMMMMMMMMMWWWWWWWWWWWWWCCCCCCCCCCSSSSSSSSSS')
        cls.aggregate = cons.OrthologFinder(stragg)
        with open((os.getcwd()+os.sep+'example_data'+os.sep+'arabidopsisCDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = cons.OrthologFinder(arabidopsisCDC48A)

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
        the_response = MagicMock(content = self.response)
        self.CDC48A.read_resp_retOMA(the_response)
        self.assertTrue(self.CDC48A.id == "ARATH09528")
        
    def test_read_resp_upOMA(self):
        the_response = MagicMock(content=self.oresponse)
        self.CDC48A.read_resp_orthoIDs(the_response)
        self.assertTrue("H3DFZ9" in self.CDC48A.ortholog_ids)

    def test_build_url_retOMA(self):
        self.lyz.sequence = 'MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV'
        self.assertEqual(self.lyz.build_url(tail='/api/sequence/?query={0}', variation=[self.lyz.sequence]), 'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')

    def test_build_url_tofasta(self):
        self.assertEqual(self.lyz.build_url(tail='/oma/vps/{0}/fasta/', variation=[self.lyz.id]), 'https://omabrowser.org/oma/vps//fasta/')

    @patch('cons.requests.get')
    def test_ologs_stragg(self, requests_mock):
        requests_mock.requests.get.return_value = None
        requests_mock.requests.get.status_code = 400
        straggFasta = self.aggregate.get_orthologs()
        self.assertTrue('There was an issue querying the database. Status code' in straggFasta)

    @patch('cons.requests.get')
    def test_ofasta_stragg(self, requests_mock):
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.ortholog_to_fasta()

    @patch('cons.requests.get')
    def test_ofasta_cdc(self, requests_mock):
        requests_mock().status_code = 200
        requests_mock().text = self.fresponse
        test = self.CDC48A.ortholog_to_fasta()
        self.assertTrue('[Arabis alpina]' in test)

    @patch('cons.requests.get')
    def test_orIDs_stragg(self, requests_mock):
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.update_orthoIDs()

    @patch('cons.requests.get')
    def test_ret_hogs_stragg(self, requests_mock):
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.retrieve_HOG_level()
            
    def test_get_fasta_seq(self):
        tester = self.lyz.get_fasta_sequence(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertFalse('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester)
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester)

    def test_indv_blk(self):
         tester = self.lyz.indv_block(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
         self.assertTrue('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester[0])
         self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester[0])

    def test_empty_input(self):
        fs = cons.OrthologFinder("")
        with self.assertRaises(cons.SequenceError) as cm:
            fs.get_orthologs()
        err = cm.exception
        self.assertEqual(str(err), 'Input sequence is empty!')

    def test_ortho_chk(self):
        output = cons.OrthologFinder.header_check("Hello World")
        self.assertEqual(output, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        with self.assertRaises(cons.SequenceError) as cm:
            cons.OrthologFinder.header_check("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        with self.assertRaises(cons.SequenceError) as sm:
            cons.OrthologFinder.header_check("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a FASTA sequence. Please try again")

    def test_ortho_already_fasta(self):
        output = cons.OrthologFinder.header_check(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(output, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")


#The following tests work with internet disconnected- ask Arnaud
    @patch('cons.requests.get')
    def test_hog_fasta(self, mock_request):
        self.lyz.hog_level = "Amniota"
        mock_request().status_code = 200
        mock_request().text = self.hresponse
        test = self.lyz.HOG_to_fasta()
        self.assertTrue("HOG:0377891.2a.2a" in test)

    @patch('cons.requests.get')
    def test_read_hog_roottrue(self, mock_request):
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=True)
        self.assertEqual(test, 'Amniota')

    @patch('cons.requests.get')
    def test_read_hog_rootfalse(self, mock_request):
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=False)
        self.assertTrue('Hystricomorpha' in test)
        self.assertTrue('Caniformia' in test)
        self.assertTrue('Gorilla gorilla gorilla' in test)

if __name__ == '__main__':
    unittest.main()
