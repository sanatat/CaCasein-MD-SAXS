from pathlib import Path
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from wire_topology import build_system


def fixtures(p):
    (p/'topol.top').write_text('#include "martini.itp"\n#include "molecule_0.itp"\n\n[ system ]\nprotein\n[ molecules ]\nmolecule_0 1\n')
    (p/'molecule_0.itp').write_text('[ moleculetype ]\n; name nrexcl\nmolecule_0 1\n[ atoms ]\n')
    for file in ('martini_v3.0.0.itp','martini_v3.0.0_ions_v1.itp','martini_v3.0.0_solvents_v1.itp'):
        (p/file).write_text('; placeholder for unit test, not a physical force field!\n')
    (p/'insane.top').write_text('#include "martini.itp"\n[ system ]\nInsane\n[ molecules ]\n; name count\nProtein   1\nW 1234\nNA+ 12\nCL- 13\n')
    (p/'system.gro').write_text('test\n2\n    1NA+    NA+    1   0.0  0.0  0.0\n    1CL-    CL-    2   0.0  0.0  0.0\n   7 7 7\n')


def test_assemble_preserves_model_and_replaces_legacy(tmp_path):
    fixtures(tmp_path)
    f=tmp_path
    report = build_system(f/'topol.top', f/'insane.top', f/'system.top', f/'system.gro')
    result=(f/'system.top').read_text()
    assert '#include "martini.itp"' not in result
    assert '#include "martini_v3.0.0.itp"' in result
    assert '#include "molecule_0.itp"' in result
    assert result.index('martini_v3.0.0.itp') < result.index('molecule_0.itp')
    assert any(line.split() == ['molecule_0', '1'] for line in result.splitlines())
    assert 'NA+ ' not in result and 'CL- ' not in result
    assert '[ molecules ]' in result and 'W                  1234' in result
    assert 'NA+' not in (f/'system.gro').read_text()
    assert report['molecules'] == [('molecule_0',1),('W',1234),('NA',12),('CL',13)]


def test_assemble_can_repair_insane_top_overwritten_in_place(tmp_path):
    fixtures(tmp_path)
    f=tmp_path
    (f/'system.top').write_text((f/'insane.top').read_text())
    build_system(f/'topol.top', f/'system.top', f/'system.top', f/'system.gro')
    assert 'martini.itp' not in (f/'system.top').read_text()
    assert 'molecule_0.itp' in (f/'system.top').read_text()


def test_missing_real_ff_fails_closed(tmp_path):
    fixtures(tmp_path)
    (tmp_path/'martini_v3.0.0.itp').unlink()
    with pytest.raises(FileNotFoundError):
        build_system(tmp_path/'topol.top',tmp_path/'insane.top',tmp_path/'system.top',tmp_path/'system.gro')
