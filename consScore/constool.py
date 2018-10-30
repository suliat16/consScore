"""
A variety of module level functions useful for handling sequence data
"""
import os
import re
from biskit.errors import BiskitError


class SequenceError(BiskitError):
    pass


def remove_protein(fasta, iden):
    """
    Removes the protein matching the entered id from the fasta string
    Args:
        fasta(str): The identification line and sequences of all the proteins, in fasta format
        iden(str): The OMA id of the closest protein match to the input sequence
    Returns:
         A string containing the proteins in fasta format, without the protein with the entered id
    """
    fasta_list = indv_block(fasta)
    for protein in fasta_list:
        if iden in protein:
            fasta_list.remove(protein)
    fasta_string = ''
    for protein in fasta_list:
        fasta_string = fasta_string + protein + os.linesep
    return fasta_string.strip()


def remove_first_protein(fasta):
    """
    Removes the first protein in a string containing proteins in fasta format
    Args:
        fasta(str): The proteins in fasta format
    Returns: A string containing the proteins in fasta format, without the first protein
    """
    fasta_list = indv_block(fasta)
    fasta_list.pop(0)
    fasta_string = ''
    for f in fasta_list:
        fasta_string = fasta_string + f + os.linesep
    return fasta_string.strip()


def header_check(sequences):
    """
    Checks to see if input string has a fasta identification line- if not,
    supplies a default identification line
    Args:
        sequences (str): The sequence(s), in fasta format
    Returns:
        The string, with default identification line, ">Input Sequence\n"
        added if needed
    """
    ##Insert mutated sequence before aligning- make it the reference sequence :)
    if sequences.startswith(">"):
        seq = sequences
    elif not sequences:
        raise SequenceError("Empty Sequence entered.")
    elif sequences[0].isalpha():
        seq = ">Input Sequence" + os.linesep + sequences
    else:
        raise SequenceError("Not a sequence. Please try again")
    return seq


def seqnwl_strip(string):
    """
    Removes the newline characters from within the sequences of the fasta
    string- Can only take a single protein as input (in fasta format)
    Args:
        string (str): The fasta sequence with excess newline characters
    Returns:
        A fasta string without the excess newline characters- Retains the
        newline character at the end of the header line
    """
    seqhead = header_check(string)
    newlist = seqhead.split(os.linesep)
    newlist = list(filter(None, newlist))
    header = newlist[0] + os.linesep
    newlist.pop(0)
    newlist.insert(0, header)
    newstring = ''.join(newlist)
    return newstring


def indv_block(st):
    """
    Return the header line and the sequence of individual constructs in a file
    Args:
        st(str): The text contained in a fasta file, as a string. Consists of a
            header, which is initiated by > and ends with a newline. Subsequent
            lines are sequence data, until another > is found.

    Returns:
        A list of strings, where each string is the construct header and sequence,
        as a single string. For example, a file containing 4 proteins would
        a list of 4 strings. Each string begins with >, and contains both the
        headers and the newline characters.
    """
    if st.startswith('>'):
        fstr = re.split('>', st)
        seq_list = []
        for f in fstr:
            if f:
                f = '>' + f
                f = f.rstrip()
                seq_list.append(f)
        return seq_list
    else:
        return [st]


def get_fasta_sequence(fasta, index=0):
    """
    Given a string in fasta format, return the sequence at the given index
    Args:
        fasta (str):The input sequence in fasta format
        index(int): For a fasta file with multiple proteins, is the zero
            indexed position of the desired protein within the file. So
            for a file containing 5 proteins, and index of 3 would correspond
            to the 4th protein
    Returns:
        The sequence of the specified protein, as a single string, with newline
        characters removed.
    """
    fstr = indv_block(st=fasta)
    fstr = fstr[index]
    fstr = fstr.splitlines()
    for f in fstr:
        if f.startswith('>'):
            fstr.remove(f)
    fstr = "".join(fstr)
    return fstr


def build_url(tail, variation, base_url):
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


def get_num(string):
    """
    Takes a string and returns the numbers in that string, decimals included
    Args:
        string(str): The string containing a number
    Returns:
        A list of all the numbers contained in that string, as floats
    """
    digits = r"-?[0-9]*\.?[0-9]+"
    parameter = re.findall(digits, string)
    return list(map(float, parameter))
