#!/usr/bin/python
# encoding: utf-8

# This program was created by Cedric L. and is part of an assignment
# for the Practical computing for bioinformatics (Fall 2014)

########################
## Required Functions ##
########################

def relative_nucleotides_freq(seq):
    if len(seq) == 0 :
        raise ValueError
    A, C, G, T = seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')
    seq_len = float(len(seq))
    return A/seq_len, C/seq_len, G/seq_len, T/seq_len

def melting_temperature_basic(seq):
    if len(seq) == 0 :
        raise ValueError
    return 4*(seq.count('G') + seq.count('C')) + 2*(seq.count('A') + seq.count('T'))

def melting_temperature_extended(seq):
    if len(seq) == 0 :
        raise ValueError
    return 64.9 + 41*(seq.count('G') + seq.count('C') - 16.4) / len(seq)

def generate_codon_dict():
    base = "ACGT"
    return {x+y+z:0 for x in base for y in base for z in base}
    
def count_codon(seq):
    if len(seq) == 0 :
        raise ValueError
    # seq length is not a multiple of 3 -> ignore last uncompleted triplet
    end = len(seq) - len(seq) % 3
    codon = generate_codon_dict()
    codon_count = 0 #how many codons is there in the sequence?

    # Splicing of the sequence into codons and counting absolute number of appearances
    i = 0
    while i < end:
        codon[seq[i:i+3]] += 1
        codon_count += 1
        i = i + 3

    # Frequencies of appearance of each codon
    key = codon.keys()
    for k in key:
        codon[k] /= float(codon_count)

    return codon
        


##################
## Main program ##
##################

seq1 = 'ACTAATGCCT'
seq2 = 'ATGAGTGAACGTCTGAGCATTACCCCGCTGGGGCCGTATATCGGCGCACAATAA'
a, c, g, t = relative_nucleotides_freq(seq1)

print "For the sequence: " + seq1 + "\n"

print "\tRelative frequencies of 'A': " + str(a)
print "\tRelative frequencies of 'C': " + str(c)
print "\tRelative frequencies of 'G': " + str(g)
print "\tRelative frequencies of 'T': " + str(t) + "\n"

print "Melting temperature of ACTAATGCCT (Basic Formula): " + \
    str(melting_temperature_basic('ACTAATGCCT'))
print "Melting temperature of ACTAATGCCT (Extended Formula): " + \
    str(melting_temperature_extended('ACTAATGCCT'))

codon_empty = generate_codon_dict()
key = codon_empty.keys()
key.sort()
print "\nHere is the initialized codons dictionary\n{ \n"
for k in key:
    print "\t '" + k + "' : " + str(codon_empty[k])
print "}"

codon_freq = count_codon(seq2)
key = codon_freq.keys()
key.sort()
print "\nHere is the codons frequencies dictionary\n{ \n"
for k in key:
    if codon_freq[k] == 0.0:
        continue
    print "\t '" + k + "' : " + str(codon_freq[k])
print "}"
