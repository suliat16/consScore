"""

"""
import consLogo
import os
import biskit.test
from unittest.mock import patch

class ConsLogoTest(biskit.test.BiskitTest):

    """Test suite for the Conservation score logo generator"""

    @classmethod
    def setUpClass(cls):
        cls.path = os.getcwd() + os.sep + 'example_data' + os.sep
        #This path is subject to change, based on where example data is stored for this file
        cls.test_logo = consLogo.OrthoLogo(cls.path + "atn1seq.fa", "PHHHQHSHIHSHLHLHQ")

    def test_get_sequence_seqinput(self):
        """Tests that given a pure sequence input, that get_sequence saves just the sequence information"""
        seq = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
        test = consLogo.OrthoLogo(seq, "UUU")
        test.get_sequence()
        self.assertEqual(seq, test.sequence)
        self.assertEqual('Protein_Sequence', test.name)

    def test_get_sequence_fileinput(self):
        """Tests that given a file as input, get_sequence retrieves the sequence inside the file-"""
        test = consLogo.OrthoLogo(self.path + 'atn1seq.fa', "UUU")
        test.get_sequence()
        self.assertTrue("PASSSAPAPPMRFPYSSSSSSSAAASSSSSSSSSSASPFPASQALPSYPHSFPPPTSLSV" in test.sequence)
        self.assertEqual('atn1seq', test.name)

    @patch('consLogo.oma.OrthologFinder.get_HOGs')
    def test_call_orthologs(self, HOG_mock):
        HOG_mock.return_value = self.path + 'atn1seq.orth'
        test = self.test_logo.call_orthologs()
        self.assertTrue(os.path.isfile(test))
        self.assertTrue(HOG_mock.called)

    @patch('consLogo.aminoCons.build_alignment')
    def test_call_alignment(self, mock_aln):
        """Tests that build_alignment is properly called"""
        mock_aln.return_value = self.path + 'atn1seq.aln'
        test = self.test_logo.call_alignment( self.path + 'atn1seq.orth')
        self.assertTrue(os.path.isfile(test))
        self.assertTrue(mock_aln.called)

    def test_find_motif_present(self):
        """tests that find motif returns the correct index if the motif is present in the alignment"""
        index = self.test_logo.find_motif(self.path + 'atn1seq.aln', "PHHHQHSHIHSHLHLHQ")
        self.assertEqual(1254, index)

    def test_find_motif_absent(self):
        """tests that find_motif returns -1 if the motif cannot be found in the alignment"""
        index = self.test_logo.find_motif(self.path + 'atn1seq.aln', "THEMOTIFISINTHEMOTIVE")
        self.assertEqual(-1, index)

    @patch('consLogo.OrthoLogo.find_motif')
    def test_run_seq2logo(self, mock_motif):
        mock_motif.return_value = 1254
        test = self.test_logo.run_seq2logo(self.path + "atn1seq.aln")
        self.assertTrue(os.path.isdir(test))
        self.assertTrue(os.path.isfile(os.getcwd() + os.sep + 'Seq2Logo' + os.sep + 'atn1seq_freq.mat'))
        self.assertTrue(mock_motif.called)

    @patch('consLogo.OrthoLogo.get_sequence')
    @patch('consLogo.OrthoLogo.call_orthologs')
    @patch('consLogo.OrthoLogo.call_alignment')
    @patch('consLogo.OrthoLogo.run_seq2logo')
    def test_get_logo(self, mock_run, mock_aln, mock_orth, mock_seq):
        """Tests that get_logo calls the right methods in sequence"""
        mock_seq.return_value = self.path + 'atn1seq.fa'
        mock_orth.return_value = self.path + 'atn1seq.orth'
        mock_aln.return_value = self.path + 'atn1seq.aln'
        mock_run.return_value = os.getcwd() + os.sep + 'Seq2Logo'
        test = consLogo.OrthoLogo(self.path + 'atn1seq.fa', 'PHHHQHSHIHSHLHLHQ')
        ret = test.get_logo()
        self.assertTrue(mock_run.called)
        self.assertTrue(mock_aln.called)
        self.assertTrue(mock_orth.called)
        self.assertTrue(mock_run.called)
        self.assertTrue(os.path.isdir(ret))

if __name__== '__main__':
    biskit.test.localTest()
