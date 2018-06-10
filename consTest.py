#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:49:39 2018

@author: suliat16

Test file for cons
"""
import unittest
import cons

class TestCons(unittest.TestCase):
    """
    """
    def setUp(self):
        strlyz = ("MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGY"
                  "NTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVAC"
                  "AKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
        self.lyz = cons.OrthologFinder(strlyz)
        stragg = ('YYYYYYYYYYYVVVVVVVVVVVVBBBBBBBBBBBAAAAAAAAAAANNNNNNNNNNN'
                  'MMMMMMMMMMMWWWWWWWWWWWWWCCCCCCCCCCSSSSSSSSSS')
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

    def test_build_url_retOMA(self):
        self.assertEqual(self.lyz.build_url(tail='/api/sequence/?query={0}', variation=[self.lyz.sequence]), 'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')

    def test_build_url_tofasta(self):
        self.assertEqual(self.lyz.build_url(tail='/oma/vps/{0}/fasta/', variation=[self.lyz.id]), 'https://omabrowser.org/oma/vps/None/fasta/')

#    def test_ologs_lyz(self):
#        lyzFasta = self.lyz.get_orthologs()
#        self.assertTrue('LYSC_HUMAN' in lyzFasta)
#
#    def test_ologs_stragg(self):
#        ##SUUUUUUUPER slow when querying the database - database times out.
#        straggFasta = self.aggregate.get_orthologs()
#        self.assertEqual(straggFasta, "Could not determine the orthologs of your sequence")
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

    def test_ofasta_lyz(self):
        lyzFasta = self.lyz.OMA_to_fasta()
        self.assertTrue('LYSC_HUMAN' in lyzFasta)

    def test_ofasta_stragg(self):
        straggFasta = self.aggregate.OMA_to_fasta()
        self.assertEqual(straggFasta, "Could not determine the orthologs of your sequence")

    def test_ofasta_arab(self):
        arabFasta = self.CDC48A.OMA_to_fasta()
        self.assertTrue('Arabidopsis' in arabFasta)

    def test_ofasta_pyk(self):
        pykFasta = self.PYK2.OMA_to_fasta()
        self.assertTrue('FAK2_HUMAN' in pykFasta)

if __name__ == '__main__':
    unittest.main()
