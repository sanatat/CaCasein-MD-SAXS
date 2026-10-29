#!/usr/bin/env python3
"""V4 robustness study for an exploratory explicit-Ca bead model.

NO Martini / no experimental validation / no biological quantitative prediction.
Reduced time and length units throughout. Compare ion initialization protocols,
co-ordination control, and analyze synthetic SAXS and independent replicas.
"""
from __future__ import annotations
import argparse, csv, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import t as t_dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'v3_explicit_ca'))
from explicit_ca_md import (BOX,DT,COORD_EPS,NPROT,NBEADS,NCHAINS, initial_state,
                            integrate,calcium_metrics,phospho_interchain_contacts,rg,sequence_features,
                            cluster_groups,unwrapped_component,forces)
# Different simulation designs, not actual concentrations or calcium chemical potentials.
CONDITIONS={
  'no_ca':(0,0.0,'none'),
  'ca8_inert_bulk':(8,0.0,'bulk'),
  'ca8_coord_bulk':(8,COORD_EPS,'bulk'),
  'ca8_coord_near':(8,COORD_EPS,'near'),
  'ca16_coord_bulk':(16,COORD_EPS,'bulk'),
}

def bulk_initial_state(seed,n_ca,box=BOX):
    """Protein initial geometry identical to legacy v3; uniformly randomized ions.

    Rejection avoids initial overlaps; placement uniformly samples available volume.
    This is *not* an equilibrated concentration or adequate binding free energy run.
    """
    protein=initial_state(seed,0,box)
    if n_ca==0:return protein
    rng=np.random.default_rng(seed+1971)
    coords=list(protein)
    for _ in range(n_ca):
        for _attempt in range(60000):
            candidate=rng.uniform(0.,box,3)
            dr=np.array(coords)-candidate
            dr-=box*np.rint(dr/box)
            if np.min(np.linalg.norm(dr,axis=1))>1.15:
                coords.append(candidate);break
        else:raise RuntimeError('ion placement failed')
    return np.array(coords)

def form_factor_debye(xyz,weights,q,box=None):
    """Debye pair formula for an isolated finite object (toy bead contrast).

    No PBC replication, interparticle structure factor or solvent subtraction.
    Caller must supply coordinates from a single unwrapped finite assembly.
    """
    xyz=np.asarray(xyz,float);weights=np.asarray(weights,float);q=np.asarray(q,float)
    dist=np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=-1)
    ii,jj=np.triu_indices(len(xyz),k=1)
    rij=dist[ii,jj]
    amplitude=np.sum(weights*weights)
    out=np.full(q.shape,amplitude,float)
    pair=2*weights[ii]*weights[jj]
    for start in range(0,len(rij),6000):
        w=pair[start:start+6000];r=rij[start:start+6000]
        out+=np.sum(w[None,:]*np.sinc(q[:,None]*r[None,:]/np.pi),axis=1)
    return out

def scattering_unwrapped_clusters(xyz,q,weights):
    """Incoherent sum over independently scattering clusters; matches v3 design."""
    comps,shifts=cluster_groups(xyz[:NPROT],NCHAINS,BOX)
    curve=np.zeros_like(q)
    for group in comps:
        xyz_clus=unwrapped_component(xyz[:NPROT],group,shifts,NBEADS)
        wi=np.concatenate([weights[c*NBEADS:(c+1)*NBEADS] for c in group])
        curve+=form_factor_debye(xyz_clus,wi,q)
    return curve

def summarize_replica_means(data,cond,field):
    x=np.array([d[field] for d in data if d['condition']==cond],float)
    m=float(x.mean());sd=float(x.std(ddof=1)) if len(x)>1 else 0.
    ci=float(t_dist.ppf(.975,len(x)-1)*sd/np.sqrt(len(x))) if len(x)>1 else None
    return m,sd,ci

def save_csv(path,rows):
    with open(path,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def run(out,steps=4000,replicas=3,stride=100,seed=20261029):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    frags,h,z,p,w,pc=sequence_features()
    h=np.tile(h,NCHAINS);z=np.tile(z,NCHAINS);p=np.tile(p,NCHAINS);w=np.tile(w,NCHAINS)
    q=np.geomspace(.045,3.2,68)
    rep_rows=[];curve_rows=[];time_rows=[];percurve={}
    for name,(nca,eps,mode) in CONDITIONS.items():
        h0=np.r_[h,np.zeros(nca)];z0=np.r_[z,np.full(nca,2.)];p0=np.r_[p,np.zeros(nca)]
        percurve[name]=[]
        for rep in range(replicas):
            sd=seed+rep*271
            x0=initial_state(sd,nca) if mode!='bulk' else bulk_initial_state(sd,nca)
            tr,energies=integrate(x0,h0,z0,p0,NPROT,NBEADS,eps,sd+31,steps,stride)
            if not np.isfinite(tr).all() or not np.isfinite(energies).all():
                raise ValueError('Non-finite simulated trajectory')
            np.savez_compressed(out/f'{name}_rep{rep+1}.npz',xyz_reduced=tr,
                potential_reduced=energies,box_reduced=BOX,dt_reduced=DT,stride=stride,
                mode=mode,coordination_proxy_kbt=eps)
            nfr=len(tr);half=max(1,nfr//2)
            xs=[];curves=[]
            for fi,x in enumerate(tr):
                occ,bridge,frac=calcium_metrics(x)
                contact=phospho_interchain_contacts(x)
                comps,shifts=cluster_groups(x[:NPROT],NCHAINS,BOX)
                main=unwrapped_component(x[:NPROT],comps[0],shifts,NBEADS)
                row=dict(condition=name,replica=rep+1,frame=fi,time_reduced=fi*stride*DT,
                    calcium_contact_count=occ,calcium_bridge_count=bridge,
                    pser_interchain_contact_count=contact,largest_cluster_chains=len(comps[0]),
                    largest_cluster_rg_reduced=rg(main),potential_reduced=float(energies[fi]))
                time_rows.append(row)
                if fi>=half:
                    xs.append(row)
                    if fi%2==0:
                        curves.append(scattering_unwrapped_clusters(x,q,w))
            calc=np.mean(curves,axis=0)
            percurve[name].append(calc/calc[0]);
            for qi,I in zip(q,calc):curve_rows.append(dict(condition=name,replica=rep+1,
                                          q_reduced=float(qi),I_synthetic=float(I),
                                          I_normalized=float(I/calc[0])))
            summary=dict(condition=name,replica=rep+1,n_calcium=nca,
                         ion_start_mode=mode,epsilon_coord_reduced=eps,
                         num_frames_late=len(xs))
            for field in ('calcium_contact_count','calcium_bridge_count',
                          'pser_interchain_contact_count','largest_cluster_chains',
                          'largest_cluster_rg_reduced'):
                summary[field+'_late_mean']=float(np.mean([d[field] for d in xs]))
                summary[field+'_early_mean']=float(np.mean([d[field] for d in time_rows if d['condition']==name and d['replica']==rep+1 and d['frame']<half]))
            rep_rows.append(summary)
            print(f'{name} rep={rep+1}: bound={summary["calcium_contact_count_late_mean"]:.3f} '
                  f'bridges={summary["calcium_bridge_count_late_mean"]:.3f}',flush=True)
    save_csv(out/'replica_summary.csv',rep_rows)
    save_csv(out/'time_metrics.csv',time_rows)
    save_csv(out/'synthetic_saxs_replicas.csv',curve_rows)
    # Evaluate bootstrap-equivalent independent-replica mean with conservative Student t intervals.
    stats=[]
    for cond in CONDITIONS:
        rec={'condition':cond,'replicas':replicas}
        for field in ('calcium_contact_count','calcium_bridge_count','pser_interchain_contact_count',
                      'largest_cluster_chains','largest_cluster_rg_reduced'):
            m,sd,hw=summarize_replica_means(rep_rows,cond,field+'_late_mean')
            rec[field+'_mean']=m;rec[field+'_sd']=sd
            rec[field+'_ci95_halfwidth']=hw
        early=np.mean([r['calcium_contact_count_early_mean'] for r in rep_rows if r['condition']==cond])
        late=rec['calcium_contact_count_mean']
        rec['early_late_delta_contacts']=late-early
        stats.append(rec)
    save_csv(out/'summary.csv',stats)
    # Pre-defined scalar distance between normalized curves; no noise model / no hypothesis test.
    contrast=[]
    a=np.array(percurve['ca8_coord_bulk'])
    for cond in CONDITIONS:
        b=np.array(percurve[cond]); d=np.sqrt(np.mean((np.log(np.maximum(a,1e-12))-np.log(np.maximum(b,1e-12)))**2,axis=1))
        contrast.append({'reference':'ca8_coord_bulk','condition':cond,
                         'rms_log_shape_diff_matched_seed_mean':float(d.mean()),
                         'rms_log_shape_diff_matched_seed_sd':float(d.std(ddof=1) if len(d)>1 else 0)})
    save_csv(out/'saxs_shape_distances.csv',contrast)
    md={
        'model':'V4 exploratory implicit-solvent bead Langevin MD with explicit mobile Ca proxy',
        'status':'NOT Martini 3; not calibrated; no experimental fit performed',
        'main_new_control':'uniform bulk vs local Ca initial placement; same protein initial structure; matched random seeds',
        'n_chains':NCHAINS,'sequence':'bovine beta-casein A2 209 aa','aa_per_bead':'~10',
        'total_protein_beads':NPROT,'reduced_box':BOX,'reduced_dt':DT,
        'steps':steps,'stride':stride,'replicas':replicas,'seed':seed,
        'sampling':'last half of each run; this does not establish equilibrium',
        'condition_definitions':CONDITIONS,
        'saxs':'incoherent sum of independently scattered unwrapped clusters; protein-only toy weights; no solvent or Ca contrasts; q in inverse reduced units',
        'uncertainty':'between-replica SD and illustrative Student-t intervals (n small, no proof of convergence)',
        'physical_units':False,
        'experimental_data_used':False,
        'limitations':['short trajectories, not equilibrated','protein bundle initialized artificially',
            'no counterions or explicit water','coarse-grain phosphate chemistry not atomistically resolved',
            'not a real milk micelle','not Martini','no physical calcium molarity','synthetic SAXS only']}
    (out/'metadata.json').write_text(json.dumps(md,indent=2,ensure_ascii=False)+'\n')
    plot(out,stats,percurve,q,rep_rows,time_rows,contrast)
    return stats,contrast

def plot(out,stats,curves,q,rep_rows,time_rows,contrast):
    names=list(CONDITIONS)
    labels={'no_ca':'0 Ca','ca8_inert_bulk':'8 Ca: no specific binding, bulk start',
            'ca8_coord_bulk':'8 Ca: specific binding, bulk start',
            'ca8_coord_near':'8 Ca: specific binding, near-site start',
            'ca16_coord_bulk':'16 Ca: specific binding, bulk start'}
    plt.rcParams.update({'figure.dpi':135,'font.size':10,'axes.spines.right':False,'axes.spines.top':False})
    fig,ax=plt.subplots(figsize=(9,5.3))
    for key in names:
        x=np.array(curves[key]);m=x.mean(axis=0);sd=x.std(axis=0,ddof=1) if len(x)>1 else np.zeros(len(q))
        ax.loglog(q,m,label=labels[key]);ax.fill_between(q,np.maximum(1e-4,m-sd),m+sd,alpha=.13)
    ax.set(xlabel='q (inverse reduced length; NOT nm^-1)',ylabel='Synthetic I(q)/I(q_min)',
           title='V4: Debye scattering of independently diffusing protein clusters')
    ax.legend(fontsize=7.8);fig.tight_layout();fig.savefig(out/'01_synthetic_saxs.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.4))
    for ax,key,title in zip(axs,['calcium_contact_count','calcium_bridge_count'],
                            ['Ca near phosphate-sensitive beads','Ca touching phosphate sites on two chains']):
        vals=[s[key+'_mean'] for s in stats];sd=[s[key+'_sd'] for s in stats]
        ax.bar(range(len(names)),vals,yerr=sd,capsize=4)
        ax.set(xticks=list(range(len(names))),xticklabels=['0 Ca','8 inert\nbulk','8 bind\nbulk','8 bind\nnear','16 bind\nbulk'],ylabel='Late-time replica mean ± SD',title=title)
    fig.suptitle('V4 robustness test — ion starting positions and specific attraction')
    fig.tight_layout();fig.savefig(out/'02_binding_control.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.6,4.5))
    for key in ('ca8_inert_bulk','ca8_coord_bulk','ca8_coord_near','ca16_coord_bulk'):
        time=[d for d in time_rows if d['condition']==key and d['replica']==1]
        ax.plot([d['time_reduced'] for d in time],
                [d['calcium_contact_count'] for d in time],label=labels[key],alpha=.85)
    ax.set(xlabel='Reduced time',ylabel='Number of Ca beads within cutoff',title='Initial ion placement sensitivity: replica 1')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'03_binding_timeseries.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.4))
    c=[d for d in contrast if d['condition']!='ca8_coord_bulk'];names2=[d['condition'] for d in c]
    ax.bar(range(len(c)),[d['rms_log_shape_diff_matched_seed_mean'] for d in c],
           yerr=[d['rms_log_shape_diff_matched_seed_sd'] for d in c],capsize=4)
    ax.set(xticks=list(range(len(c))),xticklabels=[k.replace('_','\n') for k in names2],
           ylabel='RMS(log normalized SAXS) ± replica SD',title='SAXS shape differences vs 8 Ca bulk-binding reference')
    fig.tight_layout();fig.savefig(out/'04_saxs_distances.png');plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps',type=int,default=4000)
    parser.add_argument('--stride',type=int,default=100)
    parser.add_argument('--replicas',type=int,default=3)
    parser.add_argument('--seed',type=int,default=20261029)
    parser.add_argument('--out',type=Path,default=HERE/'results')
    a=parser.parse_args()
    if a.steps<=0 or a.stride<=0 or a.steps%a.stride or a.replicas<2:
        parser.error('steps positive divisible by stride; replicas at least 2')
    run(a.out,a.steps,a.replicas,a.stride,a.seed)
