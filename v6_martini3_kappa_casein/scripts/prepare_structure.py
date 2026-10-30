"""Download a prediction from the AlphaFold DB API and extract mature kappa-casein.

An AlphaFold single-structure prediction is *not* an equilibrium ensemble,
particularly for disordered regions. This extraction does not validate the fold.
"""
import argparse
import json
from pathlib import Path
from urllib.request import urlopen, Request

from reference import PRECURSOR, MATURE, TAIL, validate_sequence

AMINO = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}

def download_predicted_pdb(target):
    url='https://alphafold.ebi.ac.uk/api/prediction/P02668'
    req=Request(url, headers={'User-Agent':'CaCasein-MD-SAXS-research-pilot/0.6'})
    with urlopen(req, timeout=45) as handle:
        models=json.load(handle)
    if not models:
        raise RuntimeError('No AlphaFold DB prediction returned for P02668.')
    entry=models[0] if isinstance(models,list) else models
    pdburl=entry.get('pdbUrl')
    if not pdburl or not pdburl.startswith('https://'):
        raise RuntimeError('AlphaFold API response did not provide an HTTPS pdbUrl.')
    with urlopen(Request(pdburl, headers={'User-Agent':'CaCasein-MD-SAXS-research-pilot/0.6'}),timeout=80) as handle:
        raw=handle.read()
    if b'ATOM ' not in raw or len(raw)>10000000:
        raise ValueError('Unexpected downloaded PDB data')
    target.write_bytes(raw)
    return {'api':url,'pdb_download_url':pdburl,'bytes':len(raw)}

def extract(raw_pdb, outdir):
    validate_sequence()
    outdir=Path(outdir)
    outdir.mkdir(parents=True,exist_ok=True)
    residues={}
    for line in Path(raw_pdb).read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21:22] not in ('A',' '):
            continue
        if line[16:17] not in (' ','A'):
            continue
        try:
            number=int(line[22:26])
        except ValueError:
            continue
        residues.setdefault(number,[]).append(line)
    missing=[i for i in range(1,191) if i not in residues]
    if missing:
        raise ValueError(f'Missing precursor positions in predicted coordinates: {missing[:12]}')
    sequence=''.join(AMINO[residues[i][0][17:20].strip()] for i in range(1,191))
    if sequence!=PRECURSOR:
        idx=next((i for i,(a,b) in enumerate(zip(sequence,PRECURSOR),1) if a!=b), None)
        raise ValueError('Predicted PDB sequence differs from documented P02668 variant A. '
                         f'First mismatch at precursor {idx}. No silent substitution allowed.')
    outputs={}
    for name, start, end in [('mature_169',22,190),('tail_64',127,190)]:
        path=outdir/(name+'.pdb')
        with path.open('w') as fd:
            fd.write('REMARK Predicted input structure, not an experimentally measured conformation\n')
            fd.write('REMARK P02668 precursor residue offset retained in provenance, output renumbered from 1\n')
            for num in range(start,end+1):
                for line in residues[num]:
                    newnum=num-start+1
                    fd.write(line[:21]+'A'+f'{newnum:4d}'+line[26:]+'\n')
            fd.write('TER\nEND\n')
        ca=[float(line[60:66]) for num in range(start,end+1) for line in residues[num] if line[12:16].strip()=='CA']
        if len(ca)!=end-start+1:
            raise ValueError('Missing CA atoms')
        outputs[name]={'path':str(path),'residues':len(ca),'mean_predicted_plddt':round(sum(ca)/len(ca),2),
                       'low_confidence_fraction_below_70':round(sum(v<70 for v in ca)/len(ca),3)}
    return outputs

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--pdb',default=None,help='Optional local P02668 full precursor PDB (avoid network)')
    parser.add_argument('--out',default='inputs')
    args=parser.parse_args()
    target=Path(args.out);target.mkdir(parents=True,exist_ok=True)
    if args.pdb:
        pdb=args.pdb;origin={'local_input':args.pdb}
    else:
        pdb=target/'P02668_predicted_full.pdb'
        origin=download_predicted_pdb(pdb)
    record={'source':origin,'structures':extract(pdb,target)}
    (target/'prediction_provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
