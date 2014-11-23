#!/usr/bin/python

# This program was created by Cedric Lood as part of an assignment
# for the Practical computing for bioinformatics at KU Leuven (Fall 2014)

import itertools

def relative_nucleotides_freq(seq):
    a = seq.count('A')
    c = seq.count('C')
    g = seq.count('G')
    t = seq.count('T')
    seq_len = float(len(seq))

    return 100*(a/seq_len), 100*(c/seq_len), 100*(g/seq_len), 100*(t/seq_len)

def melting_temperature_basic(seq):
    return 4*(seq.count('G') + seq.count('C')) + 2*(seq.count('A') + seq.count('T'))

def melting_temperature_extended(seq):
    return 64.9 + 41*(seq.count('G') + seq.count('C') - 16.4) / len(seq)

def generate_codon_dict():
    return {x+y+z:0 for x in "ACGT" for y in "ACGT" for z in "ACGT"}
    
def count_codon(codon, seq):
    # seq length is not a multiple of 3 -> ignore last uncompleted triplet
    end = len(seq) - len(seq) % 3

    i = 0
    while i < end:
        codon[seq[i:i+3]] += 1
        i = i + 3
    return codon
        

seq1 = 'ACTAATGCCT'
seq2 = 'ATGAGTGAACGTCTGAGCATTACCCCGCTGGGG'
a, c, g, t = relative_nucleotides_freq('ACTAATGCCT')

print "Relative frequencies of 'A': " + str(a)
print "Relative frequencies of 'C': " + str(c)
print "Relative frequencies of 'G': " + str(g)
print "Relative frequencies of 'T': " + str(t) + "\n"

print "Melting temperature of ACTAATGCCT (Basic Formula): " + \
    str(melting_temperature_basic('ACTAATGCCT'))
print "Melting temperature of ACTAATGCCT (Extended Formula): " + \
    str(melting_temperature_extended('ACTAATGCCT'))

codon = generate_codon_dict()
print codon

codon = count_codon(codon, seq2)
print codon
