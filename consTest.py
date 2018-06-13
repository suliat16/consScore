#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:49:39 2018

@author: suliat16

Test file for cons
"""
import unittest
from unittest.mock import Mock
from requests.models import Response
import cons

class TestCons(unittest.TestCase):
    """
    """
    def setUp(self):
        strlyz = ("MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGY"
                  "NTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVAC"
                  "AKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
        self.lyz = cons.OrthologFinder(strlyz)
        
        stragg = ('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
                  'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
                  'BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
        self.aggregate = cons.OrthologFinder(stragg)
        
        arabidopsisCDC48A = ("MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPAT"
                              "MEKLQLFRGDTILIKGKKRKDTVCIALADETCEEPKIRMNKVVRSNLRVR"
                              "LGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLKPYFLEAYRP"
                              "VRKGDLFLVRGGMRSVEFKVIETDPAEYCVVAPDTEIFCEGEPVKREDEE"
                              "RLDDVGYDDVGGVRKQMAQIRELVELPLRHPQLFKSIGVKPPKGILLYGP"
                              "PGSGKTLIARAVANETGAFFFCINGPEIMSKLAGESESNLRKAFEEAEKN"
                              "APSIIFIDEIDSIAPKREKTNGEVERRIVSQLLTLMDGLKSRAHVIVMGA"
                              "TNRPNSIDPALRRFGRFDREIDIGVPDEIGRLEVLRIHTKNMKLAEDVDLE"
                              "RISKDTHGYVGADLAALCTEAALQCIREKMDVIDLEDDSIDAEILNSMAVTN"
                              "EHFHTALGNSNPSALRETVVEVPNVSWNDIGGLENVKRELQETVQYPVEHP"
                              "EKFEKFGMSPSKGVLFYGPPGCGKTLLAKAIANECQANFISVKGPELLTMWF"
                              "GESEANVREIFDKARQSAPCVLFFDELDSIATQRGGGSGGDGGGAADRVL"
                              "NQLLTEMDGMNAKKTVFIIGATNRPDIIDSALLRPGRLDQLIYIPLPDED"
                              "SRLNIFKAALRKSPIAKDVDIGALAKYTQGFSGADITEICQRACKYAIR"
                              "ENIEKDIEKEKRRSENPEAMEEDGVDEVSEIKAAHFEESMKYARRSVSDA"
                              "DIRKYQAFAQTLQQSRGFGSEFRFENSAGSGATTGVADPFATSAAAAGDDD"
                              "DLYN")
        self.CDC48A = cons.OrthologFinder(arabidopsisCDC48A)

        humanPYK2 = ("MSGVSEPLSRVKLGTLRRPEGPAEPMVVVPVDVEKEDVRILKVCFYSNSFNPGKNF"
                      "KLVKCTVQTEIREIITSILLSGRIGPNIRLAECYGLRLKHMKSDEIHWLHPQMTVGE"
                      "VQDKYECLHVEAEWRYDLQIRYLPEDFMESLKEDRTTLLYFYQQLRNDYMQRYASKV"
                      "SEGMALQLGCLELRRFFKDMPHNALDKKSNFELLEKEVGLDLFFPKQMQENLKPKQF"
                      "RKMIQQTFQQYASLREEECVMKFFNTLAGFANIDQETYRCELIQGWNITVDLVIGPK"
                      "GIRQLTSQDAKPTCLAEFKQIRSIRCLPLEEGQAVLQLGIEGAPQALSIKTSSLAEA"
                      "ENMADLIDGYCRLQGEHQGSLIIHPRKDGEKRNSLPQIPMLNLEARRSHLSESCSIE"
                      "SDIYAEIPDETLRRPGGPQYGIAREDVVLNRILGEGFFGEVYEGVYTNHKGEKINVAV"
                      "KTCKKDCTLDNKEKFMSEAVIMKNLDHPHIVKLIGIIEEEPTWIIMELYPYGELGHY"
                      "LERNKNSLKVLTLVLYSLQICKAMAYLESINCVHRDIAVRNILVASPECVKLGDFGL"
                      "SRYIEDEDYYKASVTRLPIKWMSPESINFRRFTTASDVWMFAVCMWEILSFGKQPFF"
                      "WLENKDVIGVLEKGDRLPKPDLCPPVLYTLMTRCWDYDPSDRPRFTELVCSLSDVYQM"
                      "EKDIAMEQERNARYRTPKILEPTAFQEPPPKPSRPKYRPPPQTNLLAPKLQFQVPEG"
                      "LCASSPTLTSPMEYPSPVNSLHTPPLHRHNVFKRHSMREEDFIQPSSREEAQQLWEA"
                      "EKVKMRQILDKQQKQMVEDYQWLRQEEKSLDPMVYMNDKSPLTPEKEVGYLEFTGPP"
                      "QKPPRLGAQSIQPTANLDRTDDLVYLNVMELVRAVLELKNELCQLPPEGYVVVVKNV"
                      "GLTLRKLIGSVDDLLPSLPSSSRTEIEGTQKLLNKDLAELINKMRLAQQNAVTSLSE"
                      "ECKRQMLTASHTLAVDAKNLLDAVDQAKVLANLAHPPAE")
        self.PYK2 = cons.OrthologFinder(humanPYK2)
        
        self.fake_sequery = Mock()
        
        self.fake_sequery.json.return_value= str({
                'query':'ABC',
                'identified by': 'exact match',
                'targets':[
                    {
                        'omaid':'8e8e',
                        'chromosome': '12',
                        "locus": {
                            "start": 69348409,
                            "end": 69353219,
                            "strand": 1
                        },
                'sequence_length': 22,
                'random value number list':[1,2,3,4,45,5]
                    }
                ],
            'orthologs': 'https://omabrowser.org/api/protein/7655148/orthologs/'
        })
        self.fake_sequery.status_code = 200
        
    
        
    ## The following tests check the functions of the methods using a mock response object
    def test_read_resp_retOMA(self):
        self.lyz.read_resp_retOMA(self.fake_sequery)
        self.assertEqual(self.lyz.id, '8e8e')



## The following tests all query the API directly
#    def test_build_url_retOMA(self):
#        self.assertEqual(self.lyz.build_url(tail='/api/sequence/?query={0}', variation=[self.lyz.sequence]), 'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')
#
#    def test_build_url_tofasta(self):
#        self.assertEqual(self.lyz.build_url(tail='/oma/vps/{0}/fasta/', variation=[self.lyz.id]), 'https://omabrowser.org/oma/vps//fasta/')
#
#    def test_ologs_lyz(self):
#        lyzFasta = self.lyz.get_orthologs()
#        self.assertTrue('LYSC_HUMAN' in lyzFasta)
#
#    def test_ologs_stragg(self):
#        ##SUUUUUUUPER slow when querying the database - database times out.
#        straggFasta = self.aggregate.get_orthologs()
#        self.assertEqual(straggFasta, 'There was an issue querying the database. Status code 400')
#
#    def test_ologs_arab(self):
#        arabFasta = self.CDC48A.get_orthologs()
#        self.assertTrue('Arabidopsis' in arabFasta)
#
#    def test_ologs_pyk(self):
#        pykFasta = self.PYK2.get_orthologs()
#        self.assertTrue('FAK2_HUMAN' in pykFasta)
#
#    def test_OMAid_lyz(self):
#        self.lyz.retrieve_OMAid()
#        self.assertEqual(self.lyz.id,'HUMAN06786')
#
#    def test_OMAid_arab(self):
#        self.CDC48A.retrieve_OMAid()
#        self.assertEqual(self.CDC48A.id, 'ARATH09528')
#
#    def test_OMAid_pyk(self):
#        self.PYK2.retrieve_OMAid()
#        self.assertEqual(self.PYK2.id,'HUMAN27610')
#
#    def test_ofasta_lyz(self):
#        self.lyz.id = 'HUMAN06786'
#        lyzFasta = self.lyz.OMA_to_fasta()
#        self.assertTrue('LYSC_HUMAN' in lyzFasta)
#
#    def test_ofasta_stragg(self):
#        with self.assertRaises(ImportError):
#            self.aggregate.OMA_to_fasta()
#
#    def test_ofasta_arab(self):
#        self.CDC48A.id = 'ARATH09528'
#        arabFasta = self.CDC48A.OMA_to_fasta()
#        self.assertTrue('Arabidopsis' in arabFasta)
#
#    def test_ofasta_pyk(self):
#        self.PYK2.id = 'HUMAN27610'
#        pykFasta = self.PYK2.OMA_to_fasta()
#        self.assertTrue('FAK2_HUMAN' in pykFasta)
#        
#    def test_orIDs_lyz(self):
#        self.lyz.id = 'HUMAN06786'
#        self.lyz.update_OMA_orthoIDs()
#        self.assertTrue(len(self.lyz.ortholog_ids))
#    
#    def test_orIDs_stragg(self):
#        with self.assertRaises(ImportError):
#            self.aggregate.update_OMA_orthoIDs()
#
#    def test_orIDs_arab(self):
#        self.CDC48A.id = 'ARATH09528'
#        self.CDC48A.update_OMA_orthoIDs()
#        self.assertTrue(len(self.CDC48A.ortholog_ids))
#        
#    def test_orIDs_pyk(self):
#        self.PYK2.id = 'HUMAN27610'
#        self.PYK2.update_OMA_orthoIDs()
#        self.assertTrue(len(self.PYK2.ortholog_ids))

if __name__ == '__main__':
    unittest.main()
