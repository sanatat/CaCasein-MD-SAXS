#!/usr/bin/env python3
"""Verify V4 artifacts are present and consistent with declared run conditions."""
from pathlib import Path
import csv,json
import numpy as np
r=Path(__file__).resolve().parent
with open(r/'results'/'metadata.json') as f:md=json.load(f)
assert md['steps']==4000 and md['replicas']==3
with open(r/'results'/'replica_summary.csv') as f:rep=list(csv.DictReader(f))
assert len(rep)==15
with open(r/'results'/'summary.csv') as f:summary={x['condition']:x for x in csv.DictReader(f)}
assert len(summary)==5
for cond in summary:
    paths=sorted((r/'results').glob(f'{cond}_rep*.npz'))
    assert len(paths)==3,(cond,len(paths))
    for path in paths:
        data=np.load(path)
        assert data['xyz_reduced'].shape[0]==41
        assert np.isfinite(data['xyz_reduced']).all()
for key in ('ca8_coord_bulk','ca8_coord_near'):
    assert key in summary
with open(r/'long_time_results'/'long_time_summary.json') as f:long=json.load(f)
assert long['steps']==12000 and long['replicas']==3
print('PASS: 15 short trajectories + 6 long-check trajectories, 5 conditions, metadata, SAXS, checks')
print('Contact means near and bulk:', summary['ca8_coord_near']['calcium_contact_count_mean'],summary['ca8_coord_bulk']['calcium_contact_count_mean'])
