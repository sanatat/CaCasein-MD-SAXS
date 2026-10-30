"""Bovine kappa-casein canonical reference (UniProt P02668), variant A.

Primary sequence: https://pubchem.ncbi.nlm.nih.gov/protein/P02668
Mature chain: precursor positions 22--190 (169 amino acids).
C-terminal caseinomacropeptide: mature positions 106--169 (64 amino acids).
"""
from collections import Counter
import csv
import json
from pathlib import Path

ACCESSION = 'P02668'
PRECURSOR = (
    'MMKSFFLVVTILALTLPFLGAQEQNQEQPIRCEKDERFFSDKIAKYIPIQYVLSRYPSYGLNYYQQKPVALINNQ'
    'FLPYPYYAKPAAVRSPAQILQWQVLSNTVPAKSCQAQPTTMARHPHPHLSFMAIPPKKNQDKTEIPTINTIAS'
    'GEPTSTPTTEAVESTVATLEDSPEVIESPPEINTVQVTSTAV'
)
MATURE = PRECURSOR[21:]
TAIL = MATURE[105:]
ACIDIC = set('DE')
BASIC = set('KR')


def validate_sequence():
    assert len(PRECURSOR) == 190, len(PRECURSOR)
    assert len(MATURE) == 169, len(MATURE)
    assert len(TAIL) == 64, len(TAIL)
    assert MATURE[104:106] == 'FM', 'Expect chymosin site F105/M106 in mature numbering'
    assert PRECURSOR[147] == 'S', 'Expect UniProt phosphosite S148, precursor numbering'
    assert MATURE[126] == 'S', 'Corresponding mature phosphosite S127'
    assert set(PRECURSOR) <= set('ACDEFGHIKLMNPQRSTVWY')


def report(folder):
    validate_sequence()
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for name, seq in [('precursor_190', PRECURSOR), ('mature_169', MATURE), ('c_tail_64', TAIL)]:
        with (folder / (name + '.fasta')).open('w') as handle:
            handle.write('>P02668|bovine_kappa_casein|' + name + '\n')
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i+70] + '\n')
    with (folder / 'mature_residue_annotation.csv').open('w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['mature_position', 'precursor_position', 'residue', 'sequence_region', 'side_chain_charge_at_neutral_pH_approx', 'note'])
        for j, aa in enumerate(MATURE, 1):
            note = ('chymosin cleavage N-terminal F105' if j == 105 else
                    'chymosin cleavage C-terminal M106' if j == 106 else
                    'documented phosphoserine in canonical precursor S148' if j == 127 else '')
            w.writerow([j, j+21, aa, 'C-terminal tail' if j>=106 else 'N-terminal segment',
                        (-1 if aa in ACIDIC else 1 if aa in BASIC else 0), note])
    dat = {'accession':ACCESSION, 'precursor_length':len(PRECURSOR), 'mature_length':len(MATURE),
           'tail_length':len(TAIL), 'signal_peptide_range_precursor':[1,21],
           'mature_range_precursor':[22,190], 'chymosin_site_mature':'Phe105-Met106',
           'phosphoserine_reference_precursor':148, 'phosphoserine_reference_mature':127,
           'counts_mature':dict(Counter(MATURE)), 'counts_tail':dict(Counter(TAIL)),
           'model_scope':'unmodified, unglycosylated, unphosphorylated isolated tail only',
           'validation_status':'sequence-level only; no Martini 3 force-field run performed'}
    (folder / 'reference_report.json').write_text(json.dumps(dat, indent=2) + '\n')
    return dat

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results')
    args = ap.parse_args()
    print(json.dumps(report(args.out), indent=2))
