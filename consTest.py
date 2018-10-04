#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:49:39 2018

@author: suliat16

Test file for cons
"""
import biskit.test
import os
import oma
from requests import exceptions
from unittest.mock import patch, MagicMock


class TestCons(biskit.test.BiskitTest):
    """
    Test suite testing the behaviour of the cons module
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
        self.assertEqual(self.lyz.build_url(tail='/api/sequence/?query={0}', variation=[self.lyz.sequence]), 'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')

    @patch('cons.requests.get')
    def test_ologs_stragg(self, requests_mock):
        """Tests that a request with a bad status code raises an exception with get_orthologs"""
        requests_mock.requests.get.return_value = None
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException) as cm:
            straggFasta = self.aggregate.get_orthologs()
        err = cm.exception
        self.assertTrue('There was an issue querying the database. Status code' in str(err))

    @patch('cons.requests.get')
    def test_ofasta_stragg(self, requests_mock):
        """Tests that a bad request raises an exception in ortholog_to_fasta"""
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.ortholog_to_fasta()

    @patch('cons.requests.get')
    def test_ofasta_cdc(self, requests_mock):
        """Tests that ortholog_to_fasta correctly parses the response data"""
        requests_mock().status_code = 200
        requests_mock().text = self.fresponse
        test = self.CDC48A.ortholog_to_fasta()
        self.assertTrue('[Arabis alpina]' in test)

    @patch('cons.requests.get')
    def test_orIDs_stragg(self, requests_mock):
        """Test that update_orthoIDs returns an exception given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.update_orthoIDs()

    @patch('cons.requests.get')
    def test_ret_hogs_stragg(self, requests_mock):
        """Tests that retrieve_HOG_level throws an exception when given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.retrieve_HOG_level()
            
    def test_get_fasta_seq(self):
        """Tests that get_fasta_seq only gets the sequence of the fasta string, not the identifying line"""
        tester = oma.OrthologFinder.get_fasta_sequence(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertFalse('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester)
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester)

    def test_get_fasta_multi(self):
        """tests that get_fasta returns a list of sequences given a fasta file"""
        tester = oma.OrthologFinder.get_fasta_sequence(""">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE
>LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEASTPKVSKQGRSEEISESE
>ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQSAKKARVEEASTPKVNKQSRSEXETSAP""", index=1)
        self.assertFalse("[Loxodonta africana]" in tester)
        self.assertFalse(">PROCA12070 | ENSPCAG00000012030" in tester)
        self.assertEqual(tester, "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEASTPKVSKQGRSEEISESE")

    def test_indv_blk(self):
        """tests that indv_block can pull out individual fasta sequences with their identifying line"""
        tester = oma.OrthologFinder.indv_block(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertTrue('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester[0])
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester[0])

    def test_indv_blk_multi(self):
        """tests that indv_block can pull out individual fasta sequences from a string with multiple"""
        tester = oma.OrthologFinder.indv_block(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"  
"MKTRQNKDSMSMRSGRKKEAPGPREELRS")
        self.assertEqual(len(tester), 3)
        self.assertEqual(tester[1],(">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                    "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA" ))

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

    def test_header_chk(self):
        """test that header_check adds a header to a header-free string """
        output = oma.OrthologFinder.header_check("Hello World")
        self.assertEqual(output, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        """tests that header_check raises an exception if an empty sequence is entered"""
        with self.assertRaises(oma.SequenceError) as cm:
            oma.OrthologFinder.header_check("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        """tests that header_check raises an exception if given a non alphabetic sequence"""
        with self.assertRaises(oma.SequenceError) as sm:
            oma.OrthologFinder.header_check("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a sequence. Please try again")

    def test_ortho_already_fasta(self):
        """tests that header_check doesn't change fasta sequences that already have a header"""
        output = oma.OrthologFinder.header_check(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(output, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")

    @patch('cons.requests.get')
    def test_hog_fasta(self, mock_request):
        """Tests that HOG_to_fasta correctly parses the request response"""
        self.lyz.hog_level = "Amniota"
        mock_request().status_code = 200
        mock_request().text = self.hresponse
        test = self.lyz.HOG_to_fasta()
        self.assertTrue("HOG:0377891.2a.2a" in test)

    @patch('cons.requests.get')
    def test_read_hog_roottrue(self, mock_request):
        """tests that read_HOGid retrieves the root ID when root=True"""
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=True)
        self.assertEqual(test, 'Amniota')

    @patch('cons.requests.get')
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
