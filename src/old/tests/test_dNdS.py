import pytest
import unittest
import logging

from src.CfracdNdS import calc_Pnt, calc_PdN, calc_nomutation, calc_PdS, dNdS_ratio, dNdS_ratio_constant

#Test 1: Undefined dN/dS where containments for both protein and nucleotide are 0
#Test 2:
#Test 3: Neutral Selection
#Test 4:


def test_calc_Pnt():
    """Test calculation for P_nt"""
    #Test 1 
    assert calc_Pnt(0,7) == 1
    #Test 2 
    assert calc_Pnt(3,7) == -0.053707472093069475

def test_calc_PdN():
    """Test calculation for P_dN"""
    #Test 1
    assert calc_PdN(0,7) == 1
    #Test 2
    assert calc_PdN(2,7) == -0.10408951367381225

def test_calc_nomutation():
    """Test calculation for P_no_mutation"""
    #Test 1
    assert calc_nomutation(0,7) == 0
    #Test 2
    assert calc_nomutation(3,7) == 1.1699308127586872

def test_calc_PdS():
    """Test calculation for P_dS"""
    #Test 1
    assert calc_PdS(0,0,7) == 0
    #Test 2
    assert calc_PdS(2,3,7) == -0.0658412990848749

def test_dNdS_ratio():
    """Test calculation for dN/dS ratio using containment indexes"""
    #Test 1
    with pytest.raises(ZeroDivisionError, match='float division by zero'):
        dNdS_ratio(1,1,7)
    #Test 2
    assert dNdS_ratio(2,3,7) == 1.5809152480365283

def test_dNdS_ratio_constant():
    """Test calculation for dN/dS ratio using containment indexes and constant"""
    #Test 2
    assert dNdS_ratio_constant(2,3,7) == 0.5458765654655277


