import pytest
import unittest
import pandas as pd

#from src.reportCI import grab_containment_from_mat

""" deprecated test for deprecated function
def test_grab_containment_from_mat():
    #Test grab_containment_from_mat function
    mat = pd.DataFrame(columns=['ref','q1'])
    mat = pd.concat([mat, pd.DataFrame.from_records([{'ref':1.0,'q1':0.0}])])
    mat = pd.concat([mat, pd.DataFrame.from_records([{'ref':0.0,'q1':1.0}])])

    output=[['q1','ref', 0.0, 7]]
    df = pd.DataFrame(output,columns=['A','B','containment','ksize'])

    assert grab_containment_from_mat(mat,7).equals(df)
"""