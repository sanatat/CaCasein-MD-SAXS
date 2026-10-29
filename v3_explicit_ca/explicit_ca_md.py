#!/usr/bin/env python3
"""Exploratory explicit-Ca2+ bead Langevin MD -> synthetic SAXS for beta-casein A2.

IMPORTANT: uncalibrated, reduced-unit, implicit-solvent, sequence-informed toy model.
NOT Martini 3, NOT quantitatively predictive of Ca binding or real milk micelles.
"""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'v2_sequence_calcium'))
from casein_v2 import sequence_features, independent_cluster_saxs, cluster_groups, unwrapped_component, rg

# 4x 209-aa beta-casein A2, each coarse grained in sequence order to 21 ~10-aa beads.
NCHAINS=4
NBEADS=len(sequence_features()[0]); NPROT=NCHAINS*NBEADS
BOX=42.0  # reduced, not nm
DT=.002
BOND_K=60.0
SIGMA_WCA=1.0
WCA=2**(1/6)
CUTOFF=4.5
PHOS_R=1.36
PHOS_WIDTH=.38
COORD_EPS=3.0 # kBT per effective phosphate-site weight, NOT measured affinity
HYDRO_R=1.4
HYDRO_WIDTH=.44
ELECTRO=0.065 # screened Coulomb in effective units
SCREEN=2.2
CONDITIONS={
  'no_ca':{'n_ca':0,'coord_eps':0.0},
  'ca8_no_coord':{'n_ca':8,'coord_eps':0.0},
  'ca8_coord':{'n_ca':8,'coord_eps':COORD_EPS},
  'ca16_coord':{'n_ca':16,'coord_eps':COORD_EPS},
}

@njit(cache=True)
def forces(x,h,z,site,nprot,nbeads,coord_eps,box=BOX):
    """Energy and analytic forces for fully differentiated toy potentials.

    bonded chain neighbors: harmonic bond, all other pairs WCA repulsion,
    inter-protein nonbonded sequence-dependent attraction, screened Yukawa,
    Ca-phosphate coordination Gaussian (not a calibrated force field).
    """
    n=len(x)
    f=np.zeros_like(x)
    U=0.
    for i in range(n-1):
        for j in range(i+1,n):
            bonded=(i<nprot and j<nprot and (i//nbeads)==(j//nbeads) and j==i+1)
            dx=x[i,0]-x[j,0];dx-=box*np.rint(dx/box)
            dy=x[i,1]-x[j,1];dy-=box*np.rint(dy/box)
            dz=x[i,2]-x[j,2];dz-=box*np.rint(dz/box)
            r2=max(dx*dx+dy*dy+dz*dz,1e-9)
            r=np.sqrt(r2)
            fr=0. # F_i = fr * (x_i - x_j)
            if bonded:
                d=r-1.0
                U+=.5*BOND_K*d*d
                fr-=BOND_K*d/r
            elif r<CUTOFF:
                if r<WCA:
                    r6=(1/r2)**3
                    r12=r6*r6
                    U+=4*(r12-r6)+1
                    fr+=24*(2*r12-r6)/r2
                if i<nprot and j<nprot and i//nbeads!=j//nbeads:
                    eps=.20+1.2*h[i]*h[j]
                    v=(r-HYDRO_R)/HYDRO_WIDTH
                    G=np.exp(-.5*v*v)
                    U-=eps*G
                    fr-=eps*G*(r-HYDRO_R)/(HYDRO_WIDTH**2*r)
                if (i<nprot and j>=nprot and site[i]>0) and coord_eps>0:
                    v=(r-PHOS_R)/PHOS_WIDTH
                    G=np.exp(-.5*v*v)
                    eps=coord_eps*site[i]
                    U-=eps*G
                    fr-=eps*G*(r-PHOS_R)/(PHOS_WIDTH**2*r)
                if z[i]!=0 and z[j]!=0:
                    coupling=ELECTRO*z[i]*z[j]
                    e=np.exp(-r/SCREEN)
                    U+=coupling*e/r
                    fr+=coupling*e*(1/r2+1/(SCREEN*r))/r
            f[i,0]+=fr*dx;f[j,0]-=fr*dx
            f[i,1]+=fr*dy;f[j,1]-=fr*dy
            f[i,2]+=fr*dz;f[j,2]-=fr*dz
    return f,U


def initial_state(seed,n_ca,box=BOX):
    """Reproducible initial partially associated protein bundle, ions sampled locally.

    Calcium ions start in a local neighbourhood of phosphosites (not bulk-equilibrated).
    This construction is biased towards binding: we explicitly disclose this.
    """
    rng=np.random.default_rng(seed)
    x=np.empty((NPROT+n_ca,3))
    centers=[(-1.9,-1.9),(-1.9,1.9),(1.9,-1.9),(1.9,1.9)]
    for c,(cx,cy) in enumerate(centers):
        for b in range(NBEADS):
            x[c*NBEADS+b]=[cx,cy,(b-(NBEADS-1)/2)*1.0]
    x[:NPROT]+=rng.normal(0,.035,(NPROT,3))
    _,_,_,p,_,_=sequence_features()
    site_ids=np.flatnonzero(np.tile(p,NCHAINS)>0)
    for k in range(n_ca):
        # place near a randomly selected Ca-susceptible coarse bead, without
        # collisions with nearby ions / protein; not a physical concentration.
        for trial in range(2500):
            center=x[int(rng.choice(site_ids))]
            vec=rng.normal(size=3);vec/=np.linalg.norm(vec)
            candidate=center+vec*rng.uniform(1.28,2.5)
            prev=x[:NPROT+k]-candidate
            prev-=box*np.rint(prev/box)
            if (np.linalg.norm(prev,axis=1)>1.0).all():
                x[NPROT+k]=candidate
                break
        else:
            raise RuntimeError('Could not place non-overlapping calcium bead')
    x+=box/2
    return x

@njit(cache=True)
def integrate(x0,h,z,site,nprot,nbeads,coord_eps,seed,steps,stride,box=BOX):
    np.random.seed(seed)
    x=x0.copy()
    v=np.random.normal(0.,1.,x.shape)
    for dim in range(3):
        v[:,dim]-=np.mean(v[:,dim])
    f,U=forces(x,h,z,site,nprot,nbeads,coord_eps,box)
    nframes=steps//stride+1
    tr=np.empty((nframes,len(x),3))
    en=np.empty(nframes)
    tr[0]=x;en[0]=U
    damp=np.exp(-DT)
    thermal=np.sqrt(1.-damp*damp)
    index=1
    for step in range(1,steps+1):
        v+=.5*DT*f
        x+=.5*DT*v
        v=damp*v+thermal*np.random.normal(0.,1.,x.shape)
        x+=.5*DT*v
        f,U=forces(x,h,z,site,nprot,nbeads,coord_eps,box)
        v+=.5*DT*f
        if step%stride==0:
            tr[index]=x
            en[index]=U
            index+=1
    return tr,en


def calcium_metrics(x, cutoff=1.65):
    """Coordination contacts are geometric only: ions within cutoff of pSer bins."""
    _,_,_,p,_,_=sequence_features()
    phosph=np.flatnonzero(np.tile(p,NCHAINS)>0)
    if len(x)==NPROT:
        return 0,0,0.0
    dr=x[NPROT:,None,:]-x[None,phosph,:]
    dr-=BOX*np.rint(dr/BOX)
    dsq=(dr*dr).sum(axis=2)
    bound=dsq<cutoff**2
    occupied=int(np.any(bound,axis=1).sum())
    bridged=0
    for hits in bound:
        if len(set((phosph[np.flatnonzero(hits)]//NBEADS).tolist()))>=2:
            bridged+=1
    return occupied,bridged,occupied/(len(x)-NPROT)


def phospho_interchain_contacts(x,cutoff=1.65):
    _,_,_,p,_,_=sequence_features()
    inds=np.flatnonzero(np.tile(p,NCHAINS)>0)
    contact=0
    for i in range(len(inds)):
        for j in range(i+1,len(inds)):
            a,b=int(inds[i]),int(inds[j])
            if a//NBEADS==b//NBEADS:continue
            dr=x[a]-x[b];dr-=BOX*np.rint(dr/BOX)
            if np.dot(dr,dr)<cutoff**2:contact+=1
    return contact


def intracluster_pair_distances(x):
    """Pair-distance histogram for protein only; cluster decomposition handles PBC."""
    comp,shifts=cluster_groups(x[:NPROT],NCHAINS,BOX)
    distances=[]
    for group in comp:
        xyz=unwrapped_component(x[:NPROT],group,shifts,NBEADS)
        ii,jj=np.triu_indices(len(xyz),1)
        distances.append(np.linalg.norm(xyz[ii]-xyz[jj],axis=1))
    return np.concatenate(distances)


def save_rows(path,rows):
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def run(out,steps=2000,replicas=3,stride=80,seed=20261029):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    _,hydro,charge,sites,scatt,_=sequence_features()
    H0=np.tile(hydro,NCHAINS);Z0=np.tile(charge,NCHAINS)
    S0=np.tile(sites,NCHAINS);W=np.tile(scatt,NCHAINS)
    q=np.geomspace(.045,3.5,78)
    time_rows=[];saxs_rows=[];pr_rows=[];rep_rows=[]
    all_saxs={}; all_times={}; all_traj={}
    for cond,props in CONDITIONS.items():
        n_ca=props['n_ca'];coord=props['coord_eps']
        H=np.r_[H0,np.zeros(n_ca)]
        Z=np.r_[Z0,np.full(n_ca,2.)]
        S=np.r_[S0,np.zeros(n_ca)]
        all_saxs[cond]=[]
        print(f'{cond}: n_explicit_Ca={n_ca}, coordination_proxy={coord} kBT',flush=True)
        for rep in range(replicas):
            s=seed+rep*241
            x0=initial_state(s,n_ca)
            tr,U=integrate(x0,H,Z,S,NPROT,NBEADS,coord,s+41,steps,stride)
            if not np.isfinite(tr).all() or not np.isfinite(U).all():
                raise RuntimeError('Unstable MD: non-finite trajectory / potential')
            np.savez_compressed(out/f'{cond}_rep{rep+1}_trajectory.npz',xyz_reduced=tr,energy_reduced=U,
                t_reduced=np.arange(len(tr))*stride*DT,n_ca_explicit=n_ca,
                ion_indices_start=NPROT,box_reduced=BOX,coordination_proxy_kBT=coord)
            if rep==0:all_traj[cond]=tr
            nframe=len(tr);t=np.arange(nframe)*stride*DT
            late_start=max(1,int(nframe*.50))
            late_saxs=[];late_d=[]
            for frame,x in enumerate(tr):
                occ,bridge,frac=calcium_metrics(x)
                contacts=phospho_interchain_contacts(x)
                comps,shifts=cluster_groups(x[:NPROT],NCHAINS,BOX)
                largest=unwrapped_component(x[:NPROT],comps[0],shifts,NBEADS)
                r=rg(largest)
                row={'condition':cond,'replica':rep+1,'t_reduced':float(t[frame]),
                     'bound_Ca_count':occ,'Ca_bridge_two_chains_count':bridge,
                     'Ca_occupancy_fraction':frac,'interchain_phospho_contacts':contacts,
                     'largest_cluster_chains':len(comps[0]),'largest_cluster_rg_reduced':r,
                     'energy_reduced':float(U[frame])}
                time_rows.append(row)
                if frame>=late_start and (frame-late_start)%2==0:
                    curve,_,_=independent_cluster_saxs(x[:NPROT],W,q,NCHAINS,BOX)
                    late_saxs.append(curve)
                    late_d.append(intracluster_pair_distances(x))
            y=np.mean(late_saxs,axis=0)
            all_saxs[cond].append(y)
            for qi,I in zip(q,y):
                saxs_rows.append({'condition':cond,'replica':rep+1,'q_reduced':float(qi),
                       'I_reduced':float(I),'I_shape_normalized_at_q_min':float(I/y[0])})
            counts,bins=np.histogram(np.concatenate(late_d),bins=np.linspace(0,40,101),density=True)
            for di,dens in zip((bins[:-1]+bins[1:])/2,counts):
                pr_rows.append({'condition':cond,'replica':rep+1,'r_reduced':float(di),
                                'pair_distance_probability_density':float(dens)})
            last=[r for r in time_rows if r['condition']==cond and r['replica']==rep+1]
            late=last[late_start:]
            rep_rows.append({'condition':cond,'replica':rep+1,'n_Ca':n_ca,
                'coordination_proxy_kBT':coord,
                'mean_bound_Ca_count':float(np.mean([z['bound_Ca_count'] for z in late])),
                'mean_bridge_count':float(np.mean([z['Ca_bridge_two_chains_count'] for z in late])),
                'mean_phospho_contacts':float(np.mean([z['interchain_phospho_contacts'] for z in late])),
                'mean_largest_cluster_chains':float(np.mean([z['largest_cluster_chains'] for z in late])),
                'mean_largest_cluster_rg_reduced':float(np.mean([z['largest_cluster_rg_reduced'] for z in late]))})
    save_rows(out/'time_metrics.csv',time_rows)
    save_rows(out/'replica_summary.csv',rep_rows)
    save_rows(out/'saxs_curve_replicas.csv',saxs_rows)
    save_rows(out/'pair_distance_P_r_replicas.csv',pr_rows)
    summary=[]
    for cond,props in CONDITIONS.items():
        rec={'condition':cond,'n_Ca':props['n_ca'],'coordination_proxy_kBT':props['coord_eps'],'replicas':replicas}
        rows=[row for row in rep_rows if row['condition']==cond]
        for key in ('mean_bound_Ca_count','mean_bridge_count','mean_phospho_contacts','mean_largest_cluster_chains','mean_largest_cluster_rg_reduced'):
            z=np.array([row[key] for row in rows])
            rec[key.replace('mean_','')+'_replica_mean']=float(np.mean(z))
            rec[key.replace('mean_','')+'_replica_sd']=float(np.std(z,ddof=1)) if replicas>1 else None
        summary.append(rec)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    save_rows(out/'summary.csv',summary)
    metadata={
        'model':'Explicit divalent Ca2+ bead with screened electrostatic and toy local phosphate Gaussian; bead spring beta-casein A2',
        'status':'NONCALIBRATED EXPERIMENTAL CODE; no experimental validation, no actual calcium concentration, not Martini 3',
        'biological_model':'4 beta-casein A2 chains, 21 beads/chain, 5 phosphoserines per chain mapped into 2 sites',
        'coarse_graining':'~10 AA per bead, sequence-informed; cannot resolve calcium coordination geometry',
        'force_field':'homemade phenomenological; screened charge, WCA, harmonic bonds, sequence hydrophobicity, Ca-phosphate attraction',
        'calcium':'explicit mobile point beads with formal +2 charge (screened potential); no counterions or explicit water',
        'conditions':CONDITIONS,'box_reduced':BOX,'dt_reduced':DT,'steps':steps,'stride':stride,'replicas':replicas,'seed':seed,
        'initialization':'protein chains prepared in common aligned bundle; ions randomly initialized locally near phosphosite beads; may bias binding, not equilibrated',
        'saxs':'protein-only independently scattered clusters via Debye formula, normalized shape; no solvent excess electron density, no structure factor, no absolute nm^-1 q',
        'P_r':'pair distance distribution for protein beads only, NOT rigorous SAXS P(r) reconstruction',
        'time':'reduced units only; late 50 percent of short trajectories are not demonstrated equilibrium',
        'physical_mapping_to_experimental_units':False,'comparison_to_experimental_saxs_performed':False,
        'DO_NOT_CLAIM':['quantitative calcium concentration','full casein micelle','Martini3','SAXS experimental agreement','physical binding affinity']}
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    draw(out,q,all_saxs,rep_rows,time_rows,pr_rows,all_traj)
    return summary


def draw(out,q,all_saxs,rep_rows,time_rows,pr_rows,trajs):
    plt.rcParams.update({'figure.dpi':135,'font.size':10,'axes.spines.right':False,'axes.spines.top':False})
    labels={'no_ca':'0 Ca','ca8_no_coord':'8 Ca, no specific attraction',
            'ca8_coord':'8 Ca, coordinated (toy)','ca16_coord':'16 Ca, coordinated (toy)'}
    fig,ax=plt.subplots(figsize=(8.5,5))
    for cond,curves in all_saxs.items():
        a=np.array(curves);norm=a/a[:,[0]]
        ax.loglog(q,norm.mean(axis=0),label=labels[cond]);
        if len(norm)>1:
            sd=norm.std(axis=0,ddof=1)
            ax.fill_between(q,np.maximum(norm.mean(axis=0)-sd,1e-5),norm.mean(axis=0)+sd,alpha=.13)
    ax.set(xlabel='q (reduced inverse distance, NOT nm^-1)',ylabel='I(q) / I(q_min)',
        title='Synthetic SAXS: explicit Ca toy model; protein contrast only')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'01_synthetic_SAXS.png');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    metrics=['mean_bound_Ca_count','mean_bridge_count','mean_phospho_contacts']
    titles=['Ca beads near phosphosites','Ca beads bridging two chains','Interchain phosphosite contacts']
    for ax,key,title in zip(axs,metrics,titles):
        means=[];stds=[]
        for cond in CONDITIONS:
            arr=np.array([r[key] for r in rep_rows if r['condition']==cond]);means.append(arr.mean());stds.append(arr.std(ddof=1) if len(arr)>1 else 0)
        ax.bar(np.arange(len(CONDITIONS)),means,yerr=stds,capsize=4)
        ax.set_xticks(range(len(CONDITIONS)),['0 Ca','8 inert','8 coord','16 coord'],rotation=28)
        ax.set(title=title,ylabel='Mean per frame (replica mean ± SD)')
    fig.suptitle('Late-trajectory structural diagnostics (not equilibrated)')
    fig.tight_layout();fig.savefig(out/'02_ion_occupancy_and_bridges.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.5))
    for cond in CONDITIONS:
        rows=[r for r in time_rows if r['condition']==cond and r['replica']==1]
        ax.plot([r['t_reduced'] for r in rows],[r['Ca_bridge_two_chains_count'] for r in rows],label=labels[cond])
    ax.set(xlabel='Reduced time (not ns)',ylabel='Ca beads bridging ≥2 chains (geometric)',title='Calcium bridging time traces — replica 1')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'03_bridges_time.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.5))
    for cond in CONDITIONS:
        rows=[r for r in pr_rows if r['condition']==cond and r['replica']==1]
        ax.plot([r['r_reduced'] for r in rows],[r['pair_distance_probability_density'] for r in rows],label=labels[cond])
    ax.set(xlabel='r (reduced length)',ylabel='Normalized protein bead pair-distance density',title='Intracluster pair-distance distribution, replica 1')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'04_pair_distance.png');plt.close(fig)
    fig,axs=plt.subplots(1,len(CONDITIONS),figsize=(12,3.8),sharex=True,sharey=True)
    for ax,cond in zip(axs,CONDITIONS):
        x=trajs[cond][-1]
        for chain in range(NCHAINS):
            block=x[chain*NBEADS:(chain+1)*NBEADS]
            ax.plot(block[:,0],block[:,2],lw=1.2,alpha=.8)
        if len(x)>NPROT:ax.scatter(x[NPROT:,0],x[NPROT:,2],s=15,marker='x',color='black',label='Ca')
        ax.set(title=labels[cond],xlabel='x (reduced)')
    axs[0].set_ylabel('z (reduced)')
    fig.suptitle('Final 2D projections — conceptual coordinates; replica 1')
    fig.tight_layout();fig.savefig(out/'05_snapshots.png');plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps',type=int,default=2000)
    p.add_argument('--replicas',type=int,default=3)
    p.add_argument('--stride',type=int,default=80)
    p.add_argument('--seed',type=int,default=20261029)
    p.add_argument('--out',type=Path,default=Path('results'))
    a=p.parse_args()
    if a.steps< a.stride or a.steps%a.stride: p.error('--steps must be >= --stride and divisible by stride')
    if a.replicas<1: p.error('--replicas must be >=1')
    summary=run(a.out,a.steps,a.replicas,a.stride,a.seed)
    print('\nSummary — exploratory toy model, reduced units:')
    for row in summary:
        print(f"{row['condition']}: Ca bridges {row['bridge_count_replica_mean']:.2f} +/- {row['bridge_count_replica_sd']}, "
              f"phosphosite contacts {row['phospho_contacts_replica_mean']:.2f} +/- {row['phospho_contacts_replica_sd']}")
if __name__=='__main__':main()
