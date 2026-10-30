from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from reference import PRECURSOR
from prepare_structure import extract, AMINO
from wire_topology import wire, rename_ions

INVERSE={v:k for k,v in AMINO.items()}

def test_extraction_from_verified_synthetic_pdb(tmp_path):
    # This is a parser-only synthetic PDB test, not an actual AlphaFold structure.
    inp=tmp_path/'mock.pdb'
    lines=[]
    for n,a in enumerate(PRECURSOR,1):
        lines.append(f'ATOM   {n:4d}   CA {INVERSE[a]} A{n:4d}    {1.0*n:8.3f}{2.0:8.3f}{3.0:8.3f}{1.00:6.2f}{60.00:6.2f}           C')
    inp.write_text('\n'.join(lines)+'\nEND\n')
    result=extract(inp,tmp_path/'out')
    assert result['mature_169']['residues']==169
    assert result['tail_64']['residues']==64
    assert sum(l.startswith('ATOM') for l in (tmp_path/'out'/'tail_64.pdb').read_text().splitlines())==64

def test_topology_wiring_and_ion_names(tmp_path):
    top=tmp_path/'topol.top';top.write_text('#include "protein.itp"\n[ system ]\nT\n[ molecules ]\nProtein 1\nNA+ 1\n')
    gro=tmp_path/'system.gro';gro.write_text('title\n1\n    1NA+     NA+    1   0.000   0.000   0.000\n  4 4 4\n')
    wire(top);wire(top)
    assert top.read_text().count('martini_v3.0.0.itp')==1
    rename_ions(top,gro)
    assert 'NA+' not in top.read_text()
    assert 'NA+' not in gro.read_text()
