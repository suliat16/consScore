# ConsScore

Pipeline which retrieves the conservation score of the amino acids in sequence
given an input sequence or fasta file. Intakes the sequence or fasta string of a protein using the standard, single
letter alphabet, gets the orthologs from an external (online) database, the OMA browser. It then generates a
multiple sequence alignment (MSA) using Mafft, and calculates the conservation score of each
amino acid at each position using Rate4Site, with respect to the entered sequence.

## Usage

To run the basic, no frills pipe, going from sequence-> conservation scores

```python
import conscore ## import package
from seq2conservation import ConservationPipe ##Import pipe module

pipe = ConservationPipe(arguments) ## Create a ConservationPipe object with arguments of choice
pipe.run() ## Run the pipeline, returns either a dictionary or a biskit.ProfileCollection containing the conservation score for each amino acid in the sequence
```



## Setup



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
