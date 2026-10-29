#!/usr/bin/env python3
"""CaCasein-MD-SAXS: transparent, *generic* coarse-grained polymer pilot.

This is Langevin molecular dynamics of bead-spring chains with a tunable
inter-chain Lennard-Jones attraction, followed by single-cluster Debye SAS.
It is NOT a calibrated simulation of casein, calcium chemistry or measured SAXS.
Units are reduced: sigma=1 length unit (no nanometer calibration),
epsilon=kBT=mass=1, time=sigma*sqrt(m/kBT) (arbitrary time unit).
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit

DEFAULT_EPS = {'weak':0.10,'medium':0.65,'strong':1.50}
NCHAINS=12
BEADS=6
BOX=30.0
BOND_K=100.0
SIGMA=1.0
R_CUT=2.5
R_WCA=2.0**(1.0/6.0)
CLUSTER_CUTOFF=1.5

@njit(cache=True)
def pair_forces(pos, eps, box, nbeads, bond_k):
    """Force from harmonic bonds, WCA intra-chain, shifted LJ inter-chain."""
    n=pos.shape[0]
    force=np.zeros_like(pos)
    potential=0.0
    cutoff2=R_CUT**2
    wca2=R_WCA**2
    shift=4*eps*((1/R_CUT)**12-(1/R_CUT)**6)
    for i in range(n):
        for j in range(i+1,n):
            dx=pos[i,0]-pos[j,0]; dx-=box*np.rint(dx/box)
            dy=pos[i,1]-pos[j,1]; dy-=box*np.rint(dy/box)
            dz=pos[i,2]-pos[j,2]; dz-=box*np.rint(dz/box)
            r2=dx*dx+dy*dy+dz*dz
            if r2 < 1e-12:
                r2=1e-12
            same=i//nbeads==j//nbeads
            adjacent=same and j==i+1
            # Bonded pair also has no excluded volume, as their bond is stiff.
            if adjacent:
                r=np.sqrt(r2)
                delta=r-SIGMA
                potential+=0.5*bond_k*delta*delta
                factor=-bond_k*delta/r
            elif same:
                if r2>=wca2: continue
                inv2=1.0/r2
                inv6=inv2*inv2*inv2
                inv12=inv6*inv6
                potential+=4*(inv12-inv6)+1
                factor=24*(2*inv12-inv6)*inv2
            else:
                if r2>=cutoff2: continue
                inv2=1.0/r2
                inv6=inv2*inv2*inv2
                inv12=inv6*inv6
                potential+=4*eps*(inv12-inv6)-shift
                factor=24*eps*(2*inv12-inv6)*inv2
            force[i,0]+=factor*dx; force[j,0]-=factor*dx
            force[i,1]+=factor*dy; force[j,1]-=factor*dy
            force[i,2]+=factor*dz; force[j,2]-=factor*dz
    return force, potential


def initial_positions(seed, n_chains=NCHAINS, beads=BEADS, box=BOX):
    """Identical loose bundle by seed, with small displacements, no overlaps."""
    rng=np.random.default_rng(seed)
    if n_chains!=12:
        raise ValueError('pilot uses exactly 12 chains')
    points=[]
    for a in range(3):
        for b in range(4):
            points.append(((a-1)*1.48,(b-1.5)*1.48))
    xyz=np.zeros((n_chains*beads,3))
    for c,(x,y) in enumerate(points):
        center=np.array([x,y,0.])
        for t in range(beads):
            xyz[c*beads+t]=center+np.array([0,0,(t-(beads-1)/2)*1.0])
    xyz += rng.normal(0,0.015,xyz.shape)
    xyz += box/2
    return xyz


@njit(cache=True)
def integrate_baoab(pos_initial, epsilon, seed, steps, stride, dt, gamma, temperature, box, beads):
    """BAOAB Langevin MD: B(force/2), A(position/2), O(OU), A, B.

    Positions are unwrapped. Minimum image applies only to forces/contacts.
    """
    np.random.seed(seed)
    x=pos_initial.copy()
    v=np.random.normal(0,np.sqrt(temperature),x.shape)
    for k in range(3):
        avg=0.0
        for i in range(v.shape[0]):
            avg+=v[i,k]
        avg/=v.shape[0]
        for i in range(v.shape[0]):
            v[i,k]-=avg
    f, u=pair_forces(x,epsilon,box,beads,BOND_K)
    save_count=steps//stride+1
    traj=np.empty((save_count,x.shape[0],3),dtype=np.float64)
    times=np.empty(save_count,dtype=np.float64)
    pot=np.empty(save_count,dtype=np.float64)
    traj[0]=x; times[0]=0.; pot[0]=u
    save_idx=1
    a=np.exp(-gamma*dt)
    b=np.sqrt(temperature*(1-a*a))
    for step in range(1,steps+1):
        v+=0.5*dt*f
        x+=0.5*dt*v
        v=a*v+b*np.random.normal(0,1.0,x.shape)
        x+=0.5*dt*v
        f,u=pair_forces(x,epsilon,box,beads,BOND_K)
        v+=0.5*dt*f
        if step%stride==0:
            traj[save_idx]=x
            times[save_idx]=step*dt
            pot[save_idx]=u
            save_idx+=1
    return times,traj,pot


def contacts_clusters(pos, box=BOX, beads=BEADS, cutoff=CLUSTER_CUTOFF):
    """Connected chain components; edge if a non-bonded interchain bead
    pair is within cutoff (minimum image). Return groups and PBC shifts.
    """
    n=len(pos)//beads
    edges={i:[] for i in range(n)}
    for a in range(n):
        for b in range(a+1,n):
            vec=pos[a*beads:(a+1)*beads,None,:]-pos[None,b*beads:(b+1)*beads,:]
            nearest=vec-box*np.rint(vec/box)
            d2=np.sum(nearest*nearest,axis=2)
            if np.min(d2)<cutoff**2:
                ia,ib=np.unravel_index(np.argmin(d2),d2.shape)
                correction=box*np.rint(vec[ia,ib]/box)
                # b must shift by correction relative to a.
                edges[a].append((b,correction))
                edges[b].append((a,-correction))
    comps=[]; visited=set(); shifts={}
    for start in range(n):
        if start in visited: continue
        group=[]; stack=[start]; visited.add(start); shifts[start]=np.zeros(3)
        while stack:
            node=stack.pop(); group.append(node)
            for nxt,correction in edges[node]:
                if nxt not in visited:
                    visited.add(nxt);shifts[nxt]=shifts[node]+correction;stack.append(nxt)
        comps.append(sorted(group))
    comps.sort(key=lambda g:(-len(g),g[0]))
    return comps,shifts


def component_coords(pos, chains, shifts, beads=BEADS):
    return np.concatenate([pos[ch*beads:(ch+1)*beads]+shifts[ch] for ch in chains],axis=0)


def radius_gyration(x):
    return float(np.sqrt(((x-np.mean(x,axis=0))**2).sum(axis=1).mean()))


def debye_curve(x, qs, weights=None):
    """Orientation-averaged point-scatterer intensity (not absolute SAXS)."""
    x=np.asarray(x,float)
    w=np.ones(len(x)) if weights is None else np.asarray(weights,float)
    d=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=-1)
    ii,jj=np.triu_indices(len(x),k=1)
    products=w[ii]*w[jj]
    intensity=np.full(len(qs), np.dot(w,w),dtype=float)
    for k,q in enumerate(qs):
        intensity[k]+=2*np.dot(products,np.sinc(q*d[ii,jj]/np.pi))
    return intensity


def cluster_scattering(pos, qs, box=BOX, beads=BEADS, dr=0.25, maxdist=25):
    """Independent-cluster/dilute-mixture approximation: no intercluster
    interference. Return I(q)/N_beads, P(r) within connected clusters.
    """
    groups,shifts=contacts_clusters(pos,box,beads)
    intens=np.zeros(len(qs))
    hist=np.zeros(int(maxdist/dr))
    bins=np.linspace(0,maxdist,len(hist)+1)
    for grp in groups:
        xyz=component_coords(pos,grp,shifts,beads)
        intens+=debye_curve(xyz,qs)
        pairs=np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=-1)
        ii,jj=np.triu_indices(len(xyz),1)
        hist+=np.histogram(pairs[ii,jj],bins=bins)[0]
    return intens/len(pos),hist,groups,shifts


def eval_frame(pos, box=BOX, beads=BEADS):
    groups,shifts=contacts_clusters(pos,box,beads)
    largest=component_coords(pos,groups[0],shifts,beads)
    return {'largest_cluster_chains':len(groups[0]),
            'cluster_count':len(groups),
            'largest_cluster_rg':radius_gyration(largest),
            'weight_average_aggregate_chains':sum(len(g)**2 for g in groups)/sum(map(len,groups))}


def simulate_all(out, steps=12000, replicas=3, seed=20261029, stride=200, eps=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    eps=DEFAULT_EPS if eps is None else eps
    q=np.geomspace(0.06,6.0,125)
    frame_data=[];summary=[];curves={};pr={};snapshots={}
    for scenario,epsilon in eps.items():
        all_curves=[];all_pairs=[];statrows=[]
        for replica in range(replicas):
            rep_seed=seed+replica*71
            initial=initial_positions(rep_seed)
            times,traj,pot=integrate_baoab(initial,epsilon,rep_seed+17,steps,stride,0.004,1.0,1.0,BOX,BEADS)
            if not np.all(np.isfinite(traj)):
                raise RuntimeError('simulation diverged; lower dt')
            np.savez_compressed(out/f'{scenario}_rep{replica+1}_trajectory.npz',time_reduced=times,xyz_reduced=traj,potential_reduced=pot,box=BOX,epsilon=epsilon)
            final_idx=max(1,int(len(times)*0.55))
            for frame_idx in range(len(times)):
                metrics=eval_frame(traj[frame_idx]);rec={'scenario':scenario,'epsilon':epsilon,'replica':replica+1,'time_reduced':times[frame_idx],'potential_reduced':pot[frame_idx],**metrics}
                frame_data.append(rec)
                if frame_idx>=final_idx: statrows.append(metrics)
            # Time-average late snapshots; treat components as independently
            # tumbling noninteracting objects in dilute solution.
            for frame_idx in range(final_idx,len(times),2):
                i,p,g,s=cluster_scattering(traj[frame_idx],q)
                all_curves.append(i);all_pairs.append(p)
            snapshots[(scenario,replica)]=traj[-1].copy()
        curves[scenario]=np.mean(all_curves,axis=0)
        pr[scenario]=np.mean(all_pairs,axis=0)
        summary.append({'scenario':scenario,'epsilon_interchain_kbt':epsilon,
                        'replicas':replicas,'late_frames':len(statrows),
                        'mean_largest_cluster_chains':float(np.mean([r['largest_cluster_chains'] for r in statrows])),
                        'mean_cluster_count':float(np.mean([r['cluster_count'] for r in statrows])),
                        'mean_largest_cluster_rg_reduced':float(np.mean([r['largest_cluster_rg'] for r in statrows])),
                        'mean_weight_average_aggregate_chains':float(np.mean([r['weight_average_aggregate_chains'] for r in statrows]))})
    # I(0)/N = beads * mean weight-average number of chains (not q=qmin).
    for line in summary:
        line['I0_per_bead_from_late_frame_mean']=BEADS*line['mean_weight_average_aggregate_chains']
    with (out/'frame_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=frame_data[0].keys());writer.writeheader();writer.writerows(frame_data)
    with (out/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=summary[0].keys());writer.writeheader();writer.writerows(summary)
    with (out/'saxs_curves.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['q_reduced_inv_length']+[f'I_per_bead_{s}' for s in eps])
        for i,qq in enumerate(q):writer.writerow([qq]+[curves[s][i] for s in eps])
    with (out/'pair_distance.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['r_reduced']+[f'pair_count_{s}' for s in eps]);r=(np.arange(len(next(iter(pr.values()))))+.5)*.25
        for i,rr in enumerate(r):writer.writerow([rr]+[pr[s][i] for s in eps])
    metadata={'title':'CaCasein-MD-SAXS Pilot','date_utc':'2026-10-29','status':'in-silico toy model, not calibrated casein',
              'integrator':'BAOAB Langevin MD','chains':NCHAINS,'beads_per_chain':BEADS,'box_reduced':BOX,
              'dt_reduced':0.004,'gamma_reduced':1.0,'temperature_kbt':1.0,'bond_k_reduced':BOND_K,
              'epsilon_interchain':eps,'steps':steps,'save_stride_steps':stride,'replicas':replicas,'seed_base':seed,
              'late_analysis_start_fraction':0.55,'SAXS':'isotropic Debye, point contrast=1, ensemble dilute no interference',
              'length_units':'reduced sigma, NOT physical nanometers','time_units':'reduced, NOT physical seconds',
              'caveat':'epsilon is attraction proxy; no explicit Ca ions, casein sequences, Martini, solvent, solvent contrast or experimental data'}
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    make_plots(out,frame_data,summary,q,curves,r,pr,snapshots,eps)
    from analyze_replicas import analyze
    analyze(out)
    return summary


def make_plots(out,rows,summary,q,curves,r,pr,snapshots,eps):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
    colors=dict(zip(eps,plt.rcParams['axes.prop_cycle'].by_key()['color']))
    fig,ax=plt.subplots(figsize=(9,5))
    for scenario in eps:
        arr=np.array([[r['time_reduced'],r['largest_cluster_rg']] for r in rows if r['scenario']==scenario and r['replica']==1])
        ax.plot(arr[:,0],arr[:,1],label=f'{scenario} (eps={eps[scenario]:g})',color=colors[scenario],lw=1.8)
    ax.set(xlabel='Reduced simulation time',ylabel='Largest cluster Rg (reduced length)',title='Largest connected aggregate radius of gyration — replica 1')
    ax.legend();fig.tight_layout();fig.savefig(out/'01_rg_time.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for scenario in eps:
        arr=np.array([[r['time_reduced'],r['cluster_count'],r['largest_cluster_chains']] for r in rows if r['scenario']==scenario and r['replica']==1])
        ax.plot(arr[:,0],arr[:,2],label=f'{scenario} (eps={eps[scenario]:g})',color=colors[scenario],lw=1.8)
    ax.set(xlabel='Reduced simulation time',ylabel='Chains in largest cluster (of 12)',ylim=(0,13),title='Cluster retention — replica 1')
    ax.legend();fig.tight_layout();fig.savefig(out/'02_clusters.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for scenario in eps:
        ax.loglog(q,curves[scenario],label=f'{scenario} (eps={eps[scenario]:g})',lw=2,color=colors[scenario])
    ax.set(xlabel='q (inverse reduced length)',ylabel='I(q) / bead (arbitrary units)',title='Synthetic Debye SAS: dilute ensemble of connected clusters')
    ax.legend();fig.tight_layout();fig.savefig(out/'03_debye_saxs.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for scenario in eps:
        vals=pr[scenario]
        ax.plot(r,vals/np.sum(vals) if vals.sum()>0 else vals,label=scenario,color=colors[scenario],lw=2)
    ax.set(xlabel='Pair separation (reduced length)',ylabel='Normalized P(r) histogram (pair counts)',title='Within-cluster pair-distance distribution (late snapshots)')
    ax.legend();fig.tight_layout();fig.savefig(out/'04_pair_distance.png');plt.close(fig)
    fig=plt.figure(figsize=(12,4.5))
    for i,scenario in enumerate(eps):
        ax=fig.add_subplot(1,3,i+1,projection='3d')
        x=snapshots[(scenario,0)]
        gr,sh=contacts_clusters(x)
        aligned=np.zeros_like(x)
        for ch in range(NCHAINS):aligned[ch*BEADS:(ch+1)*BEADS]=x[ch*BEADS:(ch+1)*BEADS]+sh[ch]
        for grp in gr:
            for ch in grp:
                seq=aligned[ch*BEADS:(ch+1)*BEADS]
                ax.plot(seq[:,0],seq[:,1],seq[:,2],marker='o',ms=2.5,lw=.8)
        xyz=aligned
        cc=xyz.mean(axis=0)
        radius=max(5,0.57*np.ptp(xyz,axis=0).max())
        ax.set(xlim=(cc[0]-radius,cc[0]+radius),ylim=(cc[1]-radius,cc[1]+radius),zlim=(cc[2]-radius,cc[2]+radius),title=f'{scenario}: {len(gr)} clusters')
        ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z')
    fig.suptitle('Final-snapshot configurations (replica 1; arbitrary length units)')
    fig.tight_layout();fig.savefig(out/'05_snapshot.png');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description='Generic polymer Langevin MD -> Debye SAS pilot (not calibrated casein)')
    parser.add_argument('--steps',type=int,default=12000)
    parser.add_argument('--replicas',type=int,default=3)
    parser.add_argument('--stride',type=int,default=200)
    parser.add_argument('--seed',type=int,default=20261029)
    parser.add_argument('--out',type=str,default='results')
    args=parser.parse_args()
    if args.steps<args.stride or args.steps%args.stride or args.stride<1 or args.replicas<1:
        parser.error('steps must be positive, >= stride and exactly divisible by stride; replicas >=1')
    res=simulate_all(args.out,args.steps,args.replicas,args.seed,args.stride)
    for rec in res:
        print(f"{rec['scenario']:>6}: largest={rec['mean_largest_cluster_chains']:.2f}/12 chains; clusters={rec['mean_cluster_count']:.2f}; Rg={rec['mean_largest_cluster_rg_reduced']:.2f}; I(0)/bead={rec['I0_per_bead_from_late_frame_mean']:.2f}")
    print('Saved to',args.out)

if __name__=='__main__':main()
