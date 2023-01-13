import pytest
import unittest

from dNdS.findORFs import reverse_complement, frame_cds, listORFs, translate_ORFs, CDS_from_six_reading_frames_ORFs

def test_reverse_complement():
    """Test reverse complement function"""
    assert reverse_complement('ATGaataatgcgaagTTttc') == "GAAAACTTCGCATTATTCAT"
    
def test_frame_cds():
    """Test frame cds function"""
    assert frame_cds('ATGaataatgcgaagTTttcc', 1) == ["ATG","AAT","AAT","GCG","AAG","TTT","TCC"]
    assert frame_cds('ATGaataatgcgaagTTttcc', 2) == ["TGA","ATA","ATG","CGA","AGT","TTT","CC"]
    assert frame_cds('ATGaataatgcgaagTTttcc', 3) == ["GAA","TAA","TGC","GAA","GTT","TTC","C"]
    assert frame_cds('ATGaataatgcgaagTTttc', 4) == ["GAA","AAC","TTC","GCA","TTA","TTC","AT"]
    assert frame_cds('ATGaataatgcgaagTTttc', 5) == ["AAA","ACT","TCG","CAT","TAT","TCA","T"]
    assert frame_cds('ATGaataatgcgaagTTttc', 6) == ["AAA","CTT","CGC","ATT","ATT","CAT"]    

def test_listORFs():
    """test listORFs function"""
    assert listORFs(["ATG","GCG","AAG","TTT","TAA"]) == [["ATG","GCG","AAG","TTT"]]
    assert listORFs(["ATG","AAG","TTT","TAA","GCG","AAG","ATG","GCG","TAG","GCG"]) == [["ATG","AAG","TTT"],["ATG","GCG"]]
    assert listORFs(["ATG","ATG","AAG","TTT","TAA","GCG","AAG","ATG","GCG","TAG","GCG"]) == [["ATG","ATG","AAG","TTT"],["ATG","GCG"]]
    assert listORFs(["GCG","ATG","GCG","AAG","TTT","TAA"]) == [["ATG","GCG","AAG","TTT"]]
    assert listORFs(["ATG","AAG","TTT","TAA","ATG","GCG","TAG","GCG"]) == [["ATG","AAG","TTT"],["ATG","GCG"]]

def test_translate_ORFs():
    """test ORF translation function"""
    assert translate_ORFs(["ATG","GCG","AAG","TTT"]) == "MAKF"

def test_CDS_from_six_reading_frames_ORFs():
    """test CDS_from_six_reading_frames_ORFs function"""
    assert CDS_from_six_reading_frames_ORFs("ATGaataatgcgaagTTttcctaa") == {1:["MNNAKFS"],2:[],3:[],4:[],5:[],6:[]}
    assert CDS_from_six_reading_frames_ORFs("ATGAAGTTTTAAGCGAAGATGGCGTAGGCG") == {1:["MKF","MA"],2:[],3:[],4:[],5:[],6:[]}
    assert CDS_from_six_reading_frames_ORFs("ATGATGAAGTTTTAAGCGAAGATGGCGTAGGCG") == {1:["MMKF","MA"],2:[],3:[],4:[],5:[],6:[]}
