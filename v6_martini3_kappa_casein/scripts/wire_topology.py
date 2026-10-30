"""Assemble GROMACS Martini 3 system topology without INSANE's legacy alias.

Martinize2 generates the protein ITP and topol.top; Insane generates the
coordinates and approximate molecule counts. Neither file alone is a complete
Martini 3 system topology. We explicitly combine and validate both.
"""
import argparse
from pathlib import Path
import re

FORCEFIELD = (
    'martini_v3.0.0.itp',
    'martini_v3.0.0_ions_v1.itp',
    'martini_v3.0.0_solvents_v1.itp',
)
LEGACY = 'martini.itp'


def section_lines(text, name):
    match = re.search(r'^\s*\[\s*' + re.escape(name) + r'\s*\]\s*(?:;.*)?$', text, re.M | re.I)
    if not match:
        raise ValueError('Missing [ ' + name + ' ] section')
    rest = text[match.end():]
    end = re.search(r'^\s*\[\s*[^]]+\]\s*', rest, re.M)
    return rest[:end.start()] if end else rest


def molecule_type(itp):
    lines = section_lines(Path(itp).read_text(), 'moleculetype').splitlines()
    for line in lines:
        body = line.split(';', 1)[0].strip()
        if body and not body.startswith('#'):
            return body.split()[0]
    raise ValueError('Protein ITP lacks a [ moleculetype ] entry')


def build_system(protein_top, insane_top, output_top, gro):
    protein_top, insane_top, output_top, gro = map(Path, (protein_top, insane_top, output_top, gro))
    ptxt = protein_top.read_text()
    itps = re.findall(r'^\s*#include\s+"([^"]+)"', ptxt, re.M)
    proteins = [x for x in itps if x not in (*FORCEFIELD, LEGACY) and x.endswith('.itp')]
    if len(proteins) != 1:
        raise ValueError(f'Expected one Martinize2 protein ITP, got: {proteins}')
    protein_itp = proteins[0]
    path = protein_top.parent / protein_itp
    if not path.is_file():
        raise FileNotFoundError(f'Protein ITP not found: {path}')
    if any(not (output_top.parent / f).is_file() for f in FORCEFIELD):
        missing = [f for f in FORCEFIELD if not (output_top.parent / f).is_file()]
        raise FileNotFoundError('Missing Martini 3 parameters: ' + ', '.join(missing))
    molname = molecule_type(path)
    ins = insane_top.read_text()  # read before replacing output when paths coincide
    entries = []
    for line in section_lines(ins, 'molecules').splitlines():
        body = line.split(';', 1)[0].strip()
        if not body:
            continue
        fields = body.split()
        if len(fields) != 2 or not fields[1].isdigit():
            raise ValueError(f'Unrecognized INSANE molecule entry: {line}')
        name, count = fields
        if name == 'Protein':
            name = molname
        name = {'NA+':'NA', 'CL-':'CL'}.get(name, name)
        if name not in (molname, 'W', 'NA', 'CL'):
            raise ValueError(f'Unexpected molecule type {name} in INSANE prototype; inspect manually')
        entries.append((name, int(count)))
    if sum(n for name, n in entries if name == molname) != 1:
        raise ValueError(f'Expected exactly one protein {molname}, got {entries}')
    if not any(name == 'W' and n > 0 for name, n in entries):
        raise ValueError('Expected some Martini water molecules')
    if len(entries) != len(set(name for name, count in entries)):
        raise ValueError('Repeated molecule type entries require manual checking')

    include_order = [FORCEFIELD[0], FORCEFIELD[1], FORCEFIELD[2], protein_itp]
    header = '; Explicit Martini 3 system topology assembled from Martinize2 + Insane\n'
    header += ''.join(f'#include "{f}"\n' for f in include_order)
    header += '\n[ system ]\nUnmodified bovine kappa-casein 64-residue C-terminal tail in Martini 3 water\n\n'
    header += '[ molecules ]\n; molecule  count\n'
    header += ''.join(f'{name:<18s} {count}\n' for name, count in entries)
    output_top.write_text(header)
    # Insane legacy labels: fixed-width GRO is preserved using same-width replacement.
    original_gro = gro.read_text()
    fixed_gro = original_gro.replace('NA+', 'NA ').replace('CL-', 'CL ')
    if fixed_gro != original_gro:
        gro.write_text(fixed_gro)
    return {'protein_moleculetype': molname, 'protein_itp': protein_itp,
            'system_top': str(output_top), 'molecules': entries,
            'legacy_martini_include_removed': LEGACY not in output_top.read_text()}


def wire(top):
    """Legacy standalone helper retained for tests/backward compatibility."""
    top = Path(top)
    data = top.read_text()
    section_lines(data, 'molecules')
    section_lines(data, 'system')
    original = [line for line in data.splitlines() if line.strip().startswith('#include')]
    skip = set(FORCEFIELD) | {LEGACY}
    lines = [line for line in data.splitlines() if not any(re.search(r'\b' + re.escape(ff) + r'\b', line) for ff in skip)]
    required = ''.join(f'#include "{filename}"\n' for filename in FORCEFIELD)
    top.write_text('; Martini3 FF references\n' + required + '\n' + '\n'.join(lines) + '\n')
    return original


def rename_ions(top, gro):
    top, gro = Path(top), Path(gro)
    orig_top, orig_gro = top.read_text(), gro.read_text()
    top.write_text(orig_top.replace('NA+','NA').replace('CL-','CL'))
    gro.write_text(orig_gro.replace('NA+','NA ').replace('CL-','CL '))
    return {'top_changed':orig_top!=top.read_text(), 'gro_changed':orig_gro!=gro.read_text()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='action', required=True)
    a = sub.add_parser('assemble', help='Assemble final system.top from Martinize2 and INSANE topologies')
    a.add_argument('protein_top'); a.add_argument('insane_top'); a.add_argument('output_top'); a.add_argument('gro')
    b = sub.add_parser('wire'); b.add_argument('top')
    c = sub.add_parser('normalize'); c.add_argument('top'); c.add_argument('gro')
    args = parser.parse_args()
    if args.action == 'assemble':
        print(build_system(args.protein_top,args.insane_top,args.output_top,args.gro))
    elif args.action == 'wire':
        print(wire(args.top))
    else:
        print(rename_ions(args.top,args.gro))
