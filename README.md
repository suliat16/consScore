# ConsScore

Pipeline which retrieves the conservation score of the amino acids in sequence
given an input sequence or fasta file. Intakes the sequence or fasta string of a protein using the standard, single
letter alphabet, gets the orthologs from an external (online) database, the OMA browser. It then generates a
multiple sequence alignment (MSA) using Mafft, and calculates the conservation score of each
amino acid at each position using Rate4Site, with respect to the entered sequence.

### Usage

To run the basic, no frills pipe, going from sequence-> conservation scores

```python
import consScore ## import package
from seq2conservation import ConservationPipe ##Import pipe module

pipe = ConservationPipe(arguments) ## Create a ConservationPipe object with arguments of choice
pipe.pipe() ## Run the pipeline, returns either a dictionary or a biskit.ProfileCollection containing the conservation score for each amino acid in the sequence
```
The arguments of the pipe are as follows:
 - sequence (str): The sequence of the protein of interest, or the filepath of the fasta file containing
    the sequence of interest
 - name (str): The name of the output alignment file. If a file is given, the name of the fasta file is used.
    Name defaults to Protein Sequence
 - cache(boolean): When true, generated MSA is saved in a folder called protein sequences. When false, all
    files are deleted.
    
If the following parameters are true, the output will contain that information, in the order of the following arguments.

 - identity (boolean): The identity of the amino acid (Single letter code) at each position
 - score(boolean): The conservation scores. lower value = higher conservation.
 - qqint(boolean): QQ-INTERVAL, the confidence interval for the rate estimates. The default interval is 25-75 percentiles
 - std(boolean): The standard deviation of hte posterior rate distribution
 - gapped(boolean): MSA DATA, the number of aligned sequences having an amino acid (non-gapped) from the overall
    number of sequences at each position.

### Setup

Dependencies:
 - Requests, http://docs.python-requests.org/en/master/
 - Biskit, http://biskit.pasteur.fr/
 - Mafft, https://mafft.cbrc.jp/alignment/software/
 - Rate4Site, https://www.tau.ac.il/~itaymay/cp/rate4site.html
 
Requests and Biskit are python modules, Mafft and Rate4Site are external programs. All external programs must be on the PATH (UNIX) in order to execute. 

To install consScore:

```UNIX
$ git clone https://github.com/suliat16/consScore/
$ cd consScore
$ python setup.py install
```
Its recommended to run the test files to ensure that everything is installed and running properly. 

## Citations
OMA database:
Altenhoff A et al., The OMA orthology database in 2018: retrieving evolutionary relationships among
all domains of life through richer web and programmatic interfaces Nucleic Acids Research, 2018,
46 (D1): D477-D485 (doi:10.1093/nar/gkx1019).

Mafft:
*Katoh, Misawa, Kuma and Miyata (Nucleic Acids Res. 30:3059-3066, 2002) MAFFT: a novel method for rapid
multiple sequence alignment based on fast Fourier transform (describes the FFT-NS-1, FFT-NS-2 and FFT-NS-i
strategies)


Rate4Site:
Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference methods:
Bayesian methods are superior. Mol Biol Evol 21: 1781-1791.
"""
