#!/usr/bin/env python3
"""Independent-replica (not frame pseudo-replicate) uncertainty report."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np


def analyze(result_dir):
    result_dir=Path(result_dir)
    with (result_dir/'frame_metrics.csv').open() as f: rows=list(csv.DictReader(f))
    grouped=defaultdict(list)
    for row in rows: grouped[(row['scenario'],int(row['replica']))].append(row)
    names=('largest_cluster_chains','cluster_count','largest_cluster_rg','weight_average_aggregate_chains')
    rep_data=[]
    for (scenario,rep),records in grouped.items():
        late=records[max(1,int(len(records)*.55)):]
        entry={'scenario':scenario,'replica':rep,'late_frame_count':len(late)}
        for name in names:
            entry[f'mean_{name}']=float(np.mean([float(r[name]) for r in late]))
        rep_data.append(entry)
    report=[]
    for scenario in ('weak','medium','strong'):
        entries=[r for r in rep_data if r['scenario']==scenario]
        if not entries:continue
        item={'scenario':scenario,'independent_replicas':len(entries)}
        for name in names:
            arr=np.array([r[f'mean_{name}'] for r in entries])
            item[f'mean_{name}']=float(np.mean(arr))
            item[f'sd_across_replicas_{name}']=float(np.std(arr,ddof=1)) if len(arr)>1 else None
        report.append(item)
    with (result_dir/'per_replica_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rep_data[0]);w.writeheader();w.writerows(rep_data)
    with (result_dir/'independent_replica_uncertainty.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=report[0]);w.writeheader();w.writerows(report)
    (result_dir/'independent_replica_uncertainty.json').write_text(json.dumps(report,indent=2)+'\n')
    for row in report:
        sd=row['sd_across_replicas_largest_cluster_chains']
        sd_text=f' ± {sd:.2f}' if sd is not None else ' (SD unavailable with one replica)'
        print(f"{row['scenario']} largest cluster: {row['mean_largest_cluster_chains']:.2f}{sd_text} chains; replicas={row['independent_replicas']}")
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',nargs='?',default='results')
    args=parser.parse_args()
    analyze(args.directory)
