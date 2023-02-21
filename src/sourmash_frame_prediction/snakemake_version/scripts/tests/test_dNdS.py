import pytest
import unittest

from dNdS.reportdNdS import calc_Pnt, calc_PdN, calc_nomutation, calc_PdS, dNdS_ratio

def test_calc_Pnt():
    """Test calculation for P_nt"""
    assert calc_Pnt(0,7) == 1

def test_calc_PdN():
    """Test calculation for P_dN"""
    assert calc_PdN(0,7) == 1

def test_calc_nomutation():
    """Test calculation for P_no_mutation"""
    assert calc_nomutation(0,7) == 0

def test_calc_PdS():
    """Test calculation for P_dS"""
    assert calc_PdS(0,0,7) == 1

def test_dNdS_ratio():
    """Test calculation for dN/dS ratio using containment indexes"""
    assert dNdS_ratio(0,0,7) == 1
