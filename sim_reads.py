import random
import sys

def read_sequences_from_file(filename):
    with open(filename, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    return sequences


def introduce_mismatch(sequence, mismatch_rate):
    mutated_sequence = ""
    for base in sequence:
        if random.random() < mismatch_rate:
            mutated_base = random.choice('ATGC'.replace(base, ''))
            mutated_sequence += mutated_base
        else:
            mutated_sequence += base
    return mutated_sequence

def introduce_insertion(sequence, insertion_rate):
    mutated_sequence = ""
    for base in sequence:
        mutated_sequence += base
        if random.random() < insertion_rate:
            mutated_sequence += random.choice('ATGC')
    return mutated_sequence

def introduce_deletion(sequence, deletion_rate):
    mutated_sequence = ""
    for i in range(len(sequence)):
        if random.random() >= deletion_rate:
            mutated_sequence += sequence[i]
    return mutated_sequence

def introduce_errors(sequence, mismatch_rate, insertion_rate, deletion_rate):
    mutated_sequence = sequence
    mutated_sequence = introduce_mismatch(mutated_sequence, mismatch_rate)
    mutated_sequence = introduce_insertion(mutated_sequence, insertion_rate)
    mutated_sequence = introduce_deletion(mutated_sequence, deletion_rate)
    return mutated_sequence

# Example usage
filename = sys.argv[1]
reads = read_sequences_from_file(filename)

if len(sys.argv) == 1:
    print("sim_reads.py mismatch_rate insertion_rate deletion_rate")
    sys.exit(1)

mismatch_rate  = float(sys.argv[2]) 
insertion_rate = float(sys.argv[3]) 
deletion_rate  = float(sys.argv[4])   
subreads       = int(sys.argv[5])

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = dna_sequence[::-1]
    complement_sequence = ''.join([complement_dict[base] for base in reverse_sequence])
    return complement_sequence


reverse = False

def print_mutated_sequence(r, sequence, mismatch_rate, insertion_rate, deletion_rate):
    total_length = 0
    adapter_sequence = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"
    for i in range(subreads):
        mutated_sequence = introduce_errors(sequence, mismatch_rate, insertion_rate, deletion_rate)

        cxbits = 0b00000011 # sets presence of bother forward and reverse adapter (which is I assume detected in subread extraction on insutrment


        if i % 2 == 0:
            mutated_sequence = mutated_sequence
            cxbits = cxbits | (1 << 5) # forward read
            #+ adapter_sequence
        if i % 2 == 1:
            if reverse:
                mutated_sequence = reverse_complement(mutated_sequence)
                cxbits = cxbits | (1 << 6) # reverse read
            else:
                mutated_sequence = mutated_sequence
                cxbits = cxbits | (1 << 5) # forward read

        mutated_length = len(mutated_sequence)+len(adapter_sequence)
        print(f"simulateddata_1_1/{r}/{total_length}_{total_length+mutated_length}\t4\t*\t0\t255\t*\t*\t0\t0\t{mutated_sequence}\t{'!' * len(mutated_sequence)}\tcx:i:{int(cxbits)}\t{'ip:B:C,58,20,10' + ',10' * (mutated_length-3)}\tnp:i:1\tpw:B:C,11,9,10{',10' * (mutated_length-3)}\tqe:i:299243\tqs:i:292072\trq:f:0.8\tsn:B:f,8.05448,11.8154,3.03123,5.28106\twe:i:10798436\tws:i:10426011\tzm:i:{r}\tRG:Z:6a5a1569")
        total_length += len(mutated_sequence)

r=1
for sequence in reads:
    print_mutated_sequence(r, sequence, mismatch_rate, insertion_rate, deletion_rate)
    r += 1

