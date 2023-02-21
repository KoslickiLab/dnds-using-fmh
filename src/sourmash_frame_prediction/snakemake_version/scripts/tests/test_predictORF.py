
import pytest
import unittest


from dNdS.predictORF import get_kmer_argument

def test_get_kmer_argument():
    """Test get_kmer_argument function"""
    #For search command
    kmer_size = '7,14,21,28,35,42,49,56,63,70'
    assert get_kmer_argument('7,14,21,28,35,42,49,56,63,70') == "k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70"
    assert get_kmer_argument('7,14,21,28,35,42,49,56,63,70').split(",") == ["k=7","k=14","k=21","k=28","k=35","k=42","k=49","k=56","k=63","k=70"]
    kmer_size = '7'
    assert get_kmer_argument('7') == "k=7"
    #For prefetch command
    kmer_size = '7,14,21,28,35,42,49,56,63,70'
    assert get_kmer_argument('7,14,21,28,35,42,49,56,63,70').replace("k=","").split(",") == ["7","14","21","28","35","42","49","56","63","70"]
    kmer_size = '7'
    assert get_kmer_argument('7').replace("k=","").split(",") == ["7"]

