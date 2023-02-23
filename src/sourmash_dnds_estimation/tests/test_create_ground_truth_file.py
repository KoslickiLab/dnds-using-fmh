import pytest
import unittest

from dNdS.create_ground_truth_file import mutated_sequence_based_on_mutation_rate_p, get_coding_sequence_from_nucleotide_sequence

def test_mutated_sequence_based_on_mutation_rate_p():
    """Test mutated_sequence_based_on_mutation_rate_p"""
    seq = 'ATGCCCCCCCCCTTT'

    mutate_1 = mutated_sequence_based_on_mutation_rate_p(seq,0.2)
    mutate_2 = mutated_sequence_based_on_mutation_rate_p(seq,0.4)
    mutate_3 = mutated_sequence_based_on_mutation_rate_p(seq,0.6)
    mutate_4 = mutated_sequence_based_on_mutation_rate_p(seq,0.6)
    mutate_5 = mutated_sequence_based_on_mutation_rate_p(seq,0.8)

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

    seq = 'ATGCCCCCCCCCTTT'

    assert get_coding_sequence_from_nucleotide_sequence(seq) == ['ATG','CCC','CCC','CCC','TTT']

    