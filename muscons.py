import random
import editdistance
import sys
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_sequences_from_file(filename):
    with open(filename, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    return sequences


if len(sys.argv) == 1:
    print("sim_reads.py mismatch_rate insertion_rate deletion_rate")
    sys.exit(1)

def printaln(aln):
    aligned_seq1 = aln[0]
    aligned_seq2 = aln[1]

    alignment_string = ""

    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
           alignment_string += aligned_seq2[i]
        elif aligned_seq1[i] == '-':
           s=1
        elif aligned_seq2[i] == '-':
           alignment_string += aligned_seq2[i] 
        else:
           alignment_string += aligned_seq2[i]

    print(alignment_string)
    return alignment_string

from collections import Counter

def most_common_characters(strings):
    # Get the length of the longest string
    max_length = max(len(s) for s in strings)
    # Initialize a list to store the most common characters at each position
    most_common_chars = []

    # Iterate over each position
    for i in range(max_length):
        # Get the characters at the current position for all strings
        chars_at_position = [s[i] for s in strings if i < len(s)]
        #c = [item for item in chars_at_position if item != '-']
        #chars_at_position = c

        if len(chars_at_position) != 0:
            # Count occurrences of each character
            char_counts = Counter(chars_at_position)
            # Get the most common character at the current position
            most_common_char = char_counts.most_common(1)[0][0] if char_counts else None
            # Append the most common character to the list
            most_common_chars.append(most_common_char)
        else:
           most_common_chars.append('-')


    return ''.join(most_common_chars)

def hamming_distance(str1, str2):
    # Check if the strings are of equal length

    # Initialize the distance counter
    distance = 0

    # Iterate over the characters of the strings and compare them
    for char1, char2 in zip(str1, str2):
        if char1 != char2:
            distance += 1

    return distance

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tempfile import NamedTemporaryFile

from Bio.Align.Applications import MuscleCommandline

from Bio import SeqIO
from tempfile import NamedTemporaryFile
import os
import subprocess

def perform_alignment(sequences):
    # Write sequences to a temporary FASTA file
    with NamedTemporaryFile(mode="w") as temp_file:
        for i, seq in enumerate(sequences, start=1):
            temp_file.write(f">Seq{i}\n{seq}\n")
        temp_file.flush()

        # Create a temporary file to store aligned sequences
        aligned_file = NamedTemporaryFile(mode="w", delete=True)
        aligned_file.close()
        
        print("Input: ",temp_file.name)
        print("Output: ",aligned_file.name)

        params = ['/usr/bin/muscle','-gapopen', '-1', '-gapextend', '-1', '-matrix', './matrix.txt', '-in',temp_file.name,'-out',aligned_file.name]
        print(params)
        subprocess.call(params)

        # Read the aligned sequences from the temporary file
        aligned_sequences = list(SeqIO.parse(aligned_file.name, "fasta"))

        # Remove the temporary aligned file
        os.unlink(aligned_file.name)

        # Return the aligned sequences
        return aligned_sequences


sequences = []
with open(sys.argv[1], "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        sequences.append(record.seq)

print(sequences)

r = perform_alignment(sequences)

print(r)
# Print the aligned sequences
allstr = []
for record in r:
   print(record.seq)
   allstr.append(record.seq)

cons = most_common_characters(allstr)
cons = cons.replace('-','')
print("c: ",cons)

with open(sys.argv[2], 'w') as file:
    file.write("@read")
    file.write("\n")
    file.write(cons)
    file.write("\n")
    file.write("+read")
    file.write("\n")
    file.write('0' * len(cons))
    file.write("\n")
