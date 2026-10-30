import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import reference
from analyze_martini import rg, debye_proxy
import numpy as np

def test_lengths_and_reference():
    reference.validate_sequence()
    assert len(reference.MATURE)==169
    assert len(reference.TAIL)==64

def test_tail_matches_full_precursor():
    assert reference.TAIL == reference.PRECURSOR[126:]

def test_reference_phosphosite_positions():
    assert reference.PRECURSOR[147]=='S'
    assert reference.MATURE[126]=='S'

def test_debye_proxy_q0():
    x=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
    assert np.isclose(debye_proxy(x,np.array([0]))[0],1)
    assert 0 < debye_proxy(x,np.array([1]))[0] <= 1

def test_rg_invariant_to_translation():
    coords=np.array([[0,0,0],[2,0,0],[4,0,0]],float)
    assert np.isclose(rg(coords),np.sqrt(8/3))
    assert np.isclose(rg(coords),rg(coords+10))
