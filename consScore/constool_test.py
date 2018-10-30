#!/usr/bin/env python3

from consScore import aminoCons as am
import biskit.test
from consScore import constool

class test_constool(biskit.test.BiskitTest):

    """Test suite testing the behaviour of the constool module"""

    TAGS = [biskit.test.NORMAL]

    def test_get_num_good(self):
        """Tests that get num can retrieve a single float from a string"""
        digit = constool.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        """Tests that get_num can retrieve multiple floats and integers from a string """
        digit = constool.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        """Tests that a string containing no digits returns no digits"""
        digit = constool.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])

    def test_get_num_negative(self):
        """Tests that get_num can retrieve negative and positive numbers from within a string"""
        digit = constool.get_num("Hello, how are all 1234.56 and -42 of you today")
        self.assertEqual(digit, [1234.56, -42])

    def test_header_chk(self):
        """test that header_check adds a header to a header-free string """
        output = constool.header_check("Hello World")
        self.assertEqual(output, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        """tests that header_check raises an exception if an empty sequence is entered"""
        with self.assertRaises(constool.SequenceError) as cm:
            constool.header_check("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        """tests that header_check raises an exception if given a non alphabetic sequence"""
        with self.assertRaises(constool.SequenceError) as sm:
            constool.header_check("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a sequence. Please try again")

    def test_ortho_already_fasta(self):
        """tests that header_check doesn't change fasta sequences that already have a header"""
        output = constool.header_check(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(output, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")

    def test_get_fasta_seq(self):
        """Tests that get_fasta_seq only gets the sequence of the fasta string, not the identifying line"""
        tester = constool.get_fasta_sequence(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertFalse('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester)
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester)

    def test_get_fasta_multi(self):
        """tests that get_fasta returns a list of sequences given a fasta file"""
        tester = constool.get_fasta_sequence(""">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]
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
        tester = constool.indv_block(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertTrue('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester[0])
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester[0])

    def test_indv_blk_multi(self):
        """tests that indv_block can pull out individual fasta sequences from a string with multiple"""
        tester = constool.indv_block(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"  
"MKTRQNKDSMSMRSGRKKEAPGPREELRS")
        self.assertEqual(len(tester), 3)
        self.assertEqual(tester[1],(">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                    "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA" ))
    def test_remove_protein(self):
        """tests that remove_protein can accurately remove a  protein with the given id"""
        tester = constool.remove_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS", "LOXAF14113")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_remove_protein_nothing(self):
        """Tests that remove protein deosnt remove anythign if the entry ID is not in string"""
        tester = constool.remove_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS", "PROXY90210")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_remove_first_protein(self):
        """Tests that remove_first_protein removes the first protein in a fasta format file"""
        tester = constool.remove_first_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")
        self.assertEqual(tester, ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_seqnwl_strip(self):
        """Tests that seqnwl_strip removes the newlines from within the fasta string, but keeps the newline at the end"""
        tester = constool.seqnwl_strip(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE")

    def test_seqnwl_strip_messy(self):
        """Tests that seqnwl_strip removes the newlines from within the fasta string"""
        tester = constool.seqnwl_strip(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNK\nDSMSMRSGRKKEAPGPREEL\nRSRGRASPGGVSTSSSDGKAEKSRQTAK\nKARVEEVSAPKVSKQGRGEEIS\nESE\n")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE")

if __name__ == '__main__':
    biskit.test.localTest()
