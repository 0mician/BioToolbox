#!/usr/bin/python
# encoding: utf-8

from bioinfotools import *
import pytest

# This program was created by Cedric Lood as part of an assignment
# for the Practical computing for bioinformatics (Fall 2014)

def test_relative_nucleotides_freq():
    A, C, G, T = relative_nucleotides_freq("AAAAACCCCCTTTTTGGGGG")
    assert A == 0.25 and C == 0.25 and G == 0.25 and T == 0.25 
    A, C, G, T = relative_nucleotides_freq("AAAAATTTTT")
    assert A == 0.5 and T == 0.5 and C == 0.0 and G == 0.0
    with pytest.raises(ValueError):
        A, C, G, T = relative_nucleotides_freq("")

def test_melting_temperature_basic():
    melt_40 = melting_temperature_basic("AACCAACCTTGCTA")
    assert melt_40 == 40
    with pytest.raises(ValueError):
        melt = melting_temperature_basic("")

def test_melting_temperature_extended():
    melt_34 = melting_temperature_extended("AACCAACCTTGCTA")
    assert abs(melt_34 - 34.442857) < 0.001
    melt_noGC = melting_temperature_extended("ATATATATATTTTTT")
    assert abs(melt_noGC - 20.073333) < 0.001
    with pytest.raises(ValueError):
        melt = melting_temperature_extended("")

def test_generate_codon_dict():
    codon_table = { 'AAA' : 0, 'AAC' : 0, 'AAG' : 0, 'AAT' : 0, \
                    'ACA' : 0, 'ACC' : 0, 'ACG' : 0, 'ACT' : 0, \
                    'AGA' : 0, 'AGC' : 0, 'AGG' : 0, 'AGT' : 0, \
                    'ATA' : 0, 'ATC' : 0, 'ATG' : 0, 'ATT' : 0, \
                    'CAA' : 0, 'CAC' : 0, 'CAG' : 0, 'CAT' : 0, \
                    'CCA' : 0, 'CCC' : 0, 'CCG' : 0, 'CCT' : 0, \
                    'CGA' : 0, 'CGC' : 0, 'CGG' : 0, 'CGT' : 0, \
                    'CTA' : 0, 'CTC' : 0, 'CTG' : 0, 'CTT' : 0, \
                    'GAA' : 0, 'GAC' : 0, 'GAG' : 0, 'GAT' : 0, \
                    'GCA' : 0, 'GCC' : 0, 'GCG' : 0, 'GCT' : 0, \
                    'GGA' : 0, 'GGC' : 0, 'GGG' : 0, 'GGT' : 0, \
                    'GTA' : 0, 'GTC' : 0, 'GTG' : 0, 'GTT' : 0, \
                    'TAA' : 0, 'TAC' : 0, 'TAG' : 0, 'TAT' : 0, \
                    'TCA' : 0, 'TCC' : 0, 'TCG' : 0, 'TCT' : 0, \
                    'TGA' : 0, 'TGC' : 0, 'TGG' : 0, 'TGT' : 0, \
                    'TTA' : 0, 'TTC' : 0, 'TTG' : 0, 'TTT' : 0 } 
    gen_codon_dict = generate_codon_dict()
    assert codon_table == gen_codon_dict

def test_count_codon():
    with pytest.raises(ValueError):
        empty_seq = melting_temperature_basic("")
    codon_freq6A = count_codon("AAAAAA")
    assert codon_freq6A['AAA'] == 1.0
    assert codon_freq6A['TCA'] == 0.0
    codon_freq4T = count_codon("TTTT")
    assert codon_freq4T['TTT'] == 1.0
    assert codon_freq4T['GCA'] == 0.0
    codon_freq3A3C = count_codon("AAACCC")
    assert codon_freq3A3C['AAA'] == 0.5
    assert codon_freq3A3C['CCC'] == 0.5
