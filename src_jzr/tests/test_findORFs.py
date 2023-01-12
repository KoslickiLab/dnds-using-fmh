import pytest
import unittest
import numpy.testing as npt

@pytest.mark.parametrize(
    "test, expected",
    [
        ('CCCCCCCCCCGA', 'TCGGGGGGGGGG')
    ])
def test_reverse_complement(test, expected):
    """Test min function works for zeroes, positive integers, mix of positive/negative integers."""
    from dNdS.findORFs import reverse_complement
    npt.assert_allclose(reverse_complement(test), expected)


