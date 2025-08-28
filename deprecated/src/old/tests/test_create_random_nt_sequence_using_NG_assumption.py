import pytest
import unittest

from dNdS.create_random_nt_sequence_using_NG_assumption import positive_selection_based_on_mutation_rate_p, get_coding_sequence_from_nucleotide_sequence, negative_selection_based_on_mutation_rate_p, mutate_with_nucleotides

def test_mutate_with_nucleotides():
    assert mutate_with_nucleotides('A') == ['G','C','T']
    assert mutate_with_nucleotides('G') == ['A','G','T']
    assert mutate_with_nucleotides('T') == ['A','G','C']
    assert mutate_with_nucleotides('C') == ['A','G','T']

def test_positive_selection_based_on_mutation_rate_p():
    """Test mutated_sequence_based_on_mutation_rate_p"""
    seq = 'ATGCCCCCCCCCTTT'

    mutate_1 = positive_selection_based_on_mutation_rate_p(seq,0.2)
    mutate_2 = positive_selection_based_on_mutation_rate_p(seq,0.4)
    mutate_3 = positive_selection_based_on_mutation_rate_p(seq,0.6)
    mutate_4 = positive_selection_based_on_mutation_rate_p(seq,0.6)
    mutate_5 = positive_selection_based_on_mutation_rate_p(seq,0.8)

    assert seq != mutate_1
    assert seq != mutate_2
    assert seq != mutate_3
    assert seq != mutate_4
    assert seq != mutate_5

    assert mutate_1 != mutate_2
    assert mutate_1 != mutate_3
    assert mutate_1 != mutate_4
    assert mutate_1 != mutate_5

    assert mutate_2 != mutate_3
    assert mutate_2 != mutate_4
    assert mutate_2 != mutate_5

    assert mutate_3 != mutate_4
    assert mutate_3 != mutate_5

    assert mutate_4 != mutate_5

def test_get_coding_sequence_from_nucleotide_sequence():

    assert get_coding_sequence_from_nucleotide_sequence('ATGCCCCCCCCCTTT') == ['ATG','CCC','CCC','CCC','TTT']
    assert get_coding_sequence_from_nucleotide_sequence('ATGACTCCGGGGCCCTTT') == ['ATG','ACT','CCG','GGG','CCC','TTT']
    assert get_coding_sequence_from_nucleotide_sequence('ATGACGCCGGGCCTCTTT') == ['ATG','ACG','CCG','GGC','CTC','TTT']

def test_negative_selection_based_on_mutation_rate_p():
    """Test mutated_sequence_based_on_mutation_rate_p"""
    seq = 'ATGCCCCCCCCCTTT'

    mutate_1 = negative_selection_based_on_mutation_rate_p(seq,0.2)
    mutate_2 = negative_selection_based_on_mutation_rate_p(seq,0.4)
    mutate_3 = negative_selection_based_on_mutation_rate_p(seq,0.6)
    mutate_4 = negative_selection_based_on_mutation_rate_p(seq,0.6)
    mutate_5 = negative_selection_based_on_mutation_rate_p(seq,0.8)

    assert seq != mutate_1
    assert seq != mutate_2
    assert seq != mutate_3
    assert seq != mutate_4
    assert seq != mutate_5

    assert mutate_1 != mutate_2
    assert mutate_1 != mutate_3
    assert mutate_1 != mutate_4
    assert mutate_1 != mutate_5

    assert mutate_2 != mutate_3
    assert mutate_2 != mutate_4
    assert mutate_2 != mutate_5

    assert mutate_3 != mutate_4
    assert mutate_3 != mutate_5

    assert mutate_4 != mutate_5    