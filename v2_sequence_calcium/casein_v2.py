#!/usr/bin/env python3
"""Sequence-informed implicit-calcium beta-casein coarse-grained MD -> synthetic SAXS.

Research pilot only. 209-residue bovine beta-casein A2 mature sequence from
BIOPEP-UWM (source in README), coarsened to 21 contiguous ~10-AA blocks.
Five annotated phosphoserines provide localized Ca-sensitive attraction sites.
No explicit ions, no Martini force field, no absolute length/time, no experimental fit.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Mature Bos taurus beta-casein A2, UniProt P02666; BIOPEP-UWM record 1843.
MATURE_BETA_CASEIN_A2=(
    'RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKE'
    'MPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV'
)
# 1-indexed phosphoserine indices in mature A2 from inline modified sequence
PHOSPHOSERINES=(15, 17, 18, 19, 35)
AA_PER_BEAD=10
NCHAINS=8
BOX=54.0
BOND_K=80.0
R0=1.0
RWCA=2**(1/6)
CUTOFF=5.0
CONTACT_CUTOFF=1.65
HYDRO_R=1.38
HYDRO_WIDTH=0.43
CA_R=1.43
CA_WIDTH=0.45
SCREEN_LENGTH=2.2
ELECTRO_COUPLING=0.055
# in kBT (phenomenological bridge proxy) rather than molar [Ca2+]
BRIDGE_LEVELS={'off':0.0,'low':0.7,'high':1.4}


def sequence_features(sequence=MATURE_BETA_CASEIN_A2):
    if len(sequence)!=209 or any(sequence[i-1]!='S' for i in PHOSPHOSERINES):
        raise ValueError('Sequence or phosphorylation sites inconsistent with mature A2 source')
    fragments=[sequence[i:i+AA_PER_BEAD] for i in range(0,len(sequence),AA_PER_BEAD)]
    # hydrophobic residue enrichment only, NOT a calibrated free-energy scale.
    hydro=np.array([sum(aa in 'AVILMFWYP' for aa in segment)/len(segment) for segment in fragments],dtype=float)
    # partial effective charge proxy per coarse-grained bead at ~neutral pH; simplistic.
    charge=np.array([sum((1 if aa in 'KR' else -1 if aa in 'DE' else 0) for aa in segment) for segment in fragments],dtype=float)
    phospho_count=np.array([sum((start < p <= start+len(segment)) for p in PHOSPHOSERINES) for start,segment in zip(range(0,len(sequence),AA_PER_BEAD),fragments)],dtype=float)
    # additional phosphoserine negative charges; charge rescaled for a stable screening toy model
    charge=(charge-1.0*phospho_count)/3.0
    # contacts in each 10-residue bin: multiple pSer contribute more Ca affinity, capped.
    phos=np.minimum(1.0,phospho_count/2.0)
    electron_weight=np.array([len(seg)/AA_PER_BEAD for seg in fragments],float)
    return fragments,hydro,charge,phos,electron_weight,phospho_count


@njit(cache=True)
def total_forces(x, h, z, phos, chains, bridge, box=BOX):
    """Gradient-consistent bond + WCA + Gaussian hydrophobics/bridge + Yukawa.

    h,z,phos have one element per bead (whole system). Interaction functions
    are phenomenological in reduced units. No explicit calcium chemistry.
    """
    n=len(x)
    beadn=n//chains
    f=np.zeros_like(x)
    potential=0.0
    for i in range(n-1):
        for j in range(i+1,n):
            dx=x[i,0]-x[j,0];dx-=box*np.rint(dx/box)
            dy=x[i,1]-x[j,1];dy-=box*np.rint(dy/box)
            dz=x[i,2]-x[j,2];dz-=box*np.rint(dz/box)
            r2=dx*dx+dy*dy+dz*dz
            if r2<1e-10:r2=1e-10
            same=i//beadn==j//beadn
            if same and j==i+1:
                r=np.sqrt(r2)
                delta=r-R0
                potential+=.5*BOND_K*delta*delta
                factor=-BOND_K*delta/r
            else:
                if r2>=CUTOFF*CUTOFF:continue
                r=np.sqrt(r2)
                factor=0.0
                if r<RWCA:
                    inv2=1.0/r2
                    inv6=inv2**3
                    inv12=inv6*inv6
                    potential+=4*(inv12-inv6)+1
                    factor+=24*(2*inv12-inv6)*inv2
                # Intrachain repulsion + screened electrostatics; interchain also
                # sequence-dependent hydrophobic attraction and Ca-site bridging.
                if not same:
                    epsh=0.22+1.35*h[i]*h[j]
                    v=(r-HYDRO_R)/HYDRO_WIDTH
                    gh=np.exp(-.5*v*v)
                    potential-=epsh*gh
                    factor-=epsh*gh*(r-HYDRO_R)/(HYDRO_WIDTH**2*r)
                    if bridge>0 and phos[i]>0 and phos[j]>0:
                        v=(r-CA_R)/CA_WIDTH
                        gc=np.exp(-.5*v*v)
                        cstrength=bridge*phos[i]*phos[j]
                        potential-=cstrength*gc
                        factor-=cstrength*gc*(r-CA_R)/(CA_WIDTH**2*r)
                if z[i]!=0 and z[j]!=0:
                    coup=ELECTRO_COUPLING*z[i]*z[j]
                    exponent=np.exp(-r/SCREEN_LENGTH)
                    potential+=coup*exponent/r
                    factor+=coup*exponent*(1/r2+1/(SCREEN_LENGTH*r))/r
            f[i,0]+=factor*dx;f[j,0]-=factor*dx
            f[i,1]+=factor*dy;f[j,1]-=factor*dy
            f[i,2]+=factor*dz;f[j,2]-=factor*dz
    return f,potential


def init_bundle(seed, nchains=NCHAINS, box=BOX):
    """Same initially assembled elongated bundle for all Ca scenarios."""
    rng=np.random.default_rng(seed)
    nbeads=len(sequence_features()[0])
    centers=[(a,b) for a in (-1,1) for b in (-3,-1,1,3)]
    x=np.zeros((nchains*nbeads,3))
    for c in range(nchains):
        a,b=centers[c]
        for i in range(nbeads):
            x[c*nbeads+i]=[a*.81,b*.81,(i-(nbeads-1)/2)*R0]
    x+=rng.normal(0,.019,x.shape)
    x+=box/2
    return x


@njit(cache=True)
def integrate(x0,h,z,p,nchains,bridge,seed,steps,stride,dt=.003,gamma=1.0,temperature=1.0,box=BOX):
    np.random.seed(seed)
    x=x0.copy()
    v=np.random.normal(0,np.sqrt(temperature),x.shape)
    for j in range(3):
        v[:,j]-=np.mean(v[:,j])
    f,epot=total_forces(x,h,z,p,nchains,bridge,box)
    count=steps//stride+1
    traj=np.empty((count,len(x),3))
    times=np.empty(count)
    potentials=np.empty(count)
    traj[0]=x;times[0]=0.;potentials[0]=epot
    a=np.exp(-gamma*dt)
    b=np.sqrt(temperature*(1-a*a))
    k=1
    for step in range(1,steps+1):
        v+=.5*dt*f
        x+=.5*dt*v
        v=a*v+b*np.random.normal(0,1.,x.shape)
        x+=.5*dt*v
        f,epot=total_forces(x,h,z,p,nchains,bridge,box)
        v+=.5*dt*f
        if step%stride==0:
            traj[k]=x;times[k]=step*dt;potentials[k]=epot;k+=1
    return times,traj,potentials


def cluster_groups(x,nchains=NCHAINS,box=BOX,cutoff=CONTACT_CUTOFF):
    beadn=len(x)//nchains
    adj=[[] for _ in range(nchains)]
    for a in range(nchains):
        for b in range(a+1,nchains):
            p=x[a*beadn:(a+1)*beadn]
            t=x[b*beadn:(b+1)*beadn]
            dr=p[:,None,:]-t[None,:,:]
            dr-=box*np.rint(dr/box)
            r2=(dr*dr).sum(axis=2)
            if r2.min()<cutoff*cutoff:
                # coarse group alignment by minimum-distance contact translation
                ia,ib=np.unravel_index(np.argmin(r2),r2.shape)
                shift=box*np.rint((p[ia]-t[ib])/box)
                adj[a].append((b,shift));adj[b].append((a,-shift))
    comps=[];shifts={};seen=set()
    for seed in range(nchains):
        if seed in seen:continue
        seen.add(seed);shifts[seed]=np.zeros(3)
        stack=[seed];comp=[]
        while stack:
            i=stack.pop();comp.append(i)
            for j,shift in adj[i]:
                if j not in seen:
                    shifts[j]=shifts[i]+shift
                    seen.add(j);stack.append(j)
        comps.append(sorted(comp))
    comps.sort(key=lambda c:(-len(c),c[0]))
    return comps,shifts


def unwrapped_component(x, comp, shifts,nbeads):
    return np.concatenate([x[c*nbeads:(c+1)*nbeads]+shifts[c] for c in comp],axis=0)


def rg(x):
    return float(np.sqrt(np.mean(np.sum((x-x.mean(axis=0))**2,axis=1))))


def intensity_debye(x,weights,qs):
    """Shape-only point bead contrast; no physical electron density fitting."""
    ww=np.asarray(weights,float)
    if len(ww)!=len(x):raise ValueError('weights must match positions')
    ii,jj=np.triu_indices(len(x),1)
    distance=np.linalg.norm(x[ii]-x[jj],axis=1)
    product=ww[ii]*ww[jj]
    y=np.full(len(qs),float(np.dot(ww,ww)))
    for k,q in enumerate(qs):
        y[k]+=2*float(np.dot(product,np.sinc(q*distance/np.pi)))
    return y


def independent_cluster_saxs(x,w,qs,nchains=NCHAINS,box=BOX):
    """SAS of disconnected objects (dilute mixture approximation)."""
    comps,shifts=cluster_groups(x,nchains,box)
    nbeads=len(x)//nchains
    y=np.zeros(len(qs))
    for c in comps:
        xyz=unwrapped_component(x,c,shifts,nbeads)
        ws=np.concatenate([w[chain*nbeads:(chain+1)*nbeads] for chain in c])
        y+=intensity_debye(xyz,ws,qs)
    return y/w.sum(),comps,shifts


def feature_dump(out):
    chunks,hydro,charge,phos,weights,nphos=sequence_features()
    with (out/'sequence_features.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['bead_index_1based','residues_1based','sequence','hydrophobic_fraction_proxy','screened_charge_proxy','phosphoserine_count','calcium_site_weight','scattering_bead_weight_proxy'])
        for i,seg in enumerate(chunks):
            writer.writerow([i+1,f'{i*AA_PER_BEAD+1}-{i*AA_PER_BEAD+len(seg)}',seg,hydro[i],charge[i],int(nphos[i]),phos[i],weights[i]])
    return hydro,charge,phos,weights


def run(out,steps=3500,replicas=3,stride=100,seed=20261029,conditions=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    h,z,p,w=feature_dump(out)
    nchains=NCHAINS
    H=np.tile(h,nchains);Z=np.tile(z,nchains);P=np.tile(p,nchains);W=np.tile(w,nchains)
    qs=np.geomspace(.025,3.5,105)
    rows=[];curve_rows=[];summary=[];trajectories={};curves={}
    levels=BRIDGE_LEVELS if conditions is None else conditions
    for cond,bridge in levels.items():
        replica_curves=[]
        print(f'Running {cond}: bridge={bridge} kBT ...',flush=True)
        for rep in range(replicas):
            rseed=seed+rep*397
            initial=init_bundle(rseed)
            times,traj,potential=integrate(initial,H,Z,P,nchains,bridge,rseed+17,steps,stride)
            if not(np.isfinite(traj).all() and np.isfinite(potential).all()):raise RuntimeError('nonfinite positions or energies')
            if np.max(np.linalg.norm(np.diff(traj,axis=2),axis=2))>1e7:raise RuntimeError('unstable integrator')
            np.savez_compressed(out/f'{cond}_rep{rep+1}_trajectory.npz',time_reduced=times,xyz_reduced=traj,potential_reduced=potential,calcium_bridge_proxy_kbt=bridge,box_reduced=BOX)
            start=max(1,int(len(times)*.55))
            late_curves=[]
            for k in range(len(times)):
                comps,shifts=cluster_groups(traj[k],nchains)
                largest=unwrapped_component(traj[k],comps[0],shifts,len(h))
                metric={'condition':cond,'bridge_kbt':bridge,'replica':rep+1,'t_reduced':times[k],
                        'potential_reduced':potential[k],'largest_cluster_chains':len(comps[0]),
                        'cluster_count':len(comps),'largest_cluster_rg_reduced':rg(largest)}
                rows.append(metric)
                if k>=start and (k-start)%2==0:
                    curve,_,_=independent_cluster_saxs(traj[k],W,qs,nchains)
                    late_curves.append(curve)
            replica_curves.append(np.mean(late_curves,axis=0))
            trajectories[(cond,rep)]=traj
        curves[cond]=np.array(replica_curves)
        reps=[[] for _ in range(replicas)]
        for row in rows:
            if row['condition']==cond and row['t_reduced']>=steps*.003*.55:
                reps[row['replica']-1].append(row)
        repmeans={key:np.array([np.mean([x[key] for x in r]) for r in reps]) for key in ('largest_cluster_chains','cluster_count','largest_cluster_rg_reduced')}
        summary.append({'condition':cond,'bridge_strength_proxy_kbt':bridge,'replicas':replicas,'steps':steps,'late_window_fraction':.45,
                        **{f'{k}_mean':float(np.mean(v)) for k,v in repmeans.items()},
                        **{f'{k}_replica_sd':float(np.std(v,ddof=1)) if replicas>1 else None for k,v in repmeans.items()},
                        'mean_scattering_I0_per_scattering_weight':float(curves[cond][:,0].mean())})
    # CSV output and plots
    with (out/'time_metrics.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    with (out/'summary.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(summary[0]));wr.writeheader();wr.writerows(summary)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (out/'saxs_curve_replicas.csv').open('w',newline='') as f:
        wrc=csv.writer(f);wrc.writerow(['q_reduced','condition','replica','I_q_per_bead_scattering_weight','I_q_over_I0_model'])
        for cond,matrix in curves.items():
            for rep,y in enumerate(matrix):
                for q,yq in zip(qs,y):
                    wrc.writerow([q,cond,rep+1,yq,yq/matrix[rep,0]])
    meta={'kind':'Sequence-informed phenomenological Langevin bead-spring MD; NOT calibrated casein MD',
          'protein':'Bovine beta-casein A2 mature chain 209 amino acids (BIOPEP-UWM ID 1843; UniProt P02666)',
          'phosphoserines_mature_1based':PHOSPHOSERINES,'aa_per_bead':AA_PER_BEAD,'chains':nchains,'beads_per_chain':len(h),'box_reduced':BOX,
          'conditions_bridge_kBT':levels,'bridge_potential':'attractive gaussian between interchain phosphoserine-rich beads',
          'hydrophobicity':'sequence enriched, heuristic gaussian attraction','electrostatics':'charge proxy screened Yukawa',
          'thermostat':'Langevin BAOAB','initial_conditions':'all conditions start as a common preassembled bundle (retention study)',
          'dt_reduced':.003,'gamma':1.0,'kBT':1.0,'steps':steps,'replicas':replicas,'seed_base':seed,
          'SAXS':'unitless solvent contrast proxy, intracluster Debye forward calculation, no structure factor',
          'experimental_data_used':False,'experimental_validation':False,'physical_nm_or_ns_mapping':False,
          'limitations':'Not Martini/GROMACS, no explicit calcium ions, no ion binding thermodynamics, no milk micelle, no physical time/length units. Ca settings are effective attraction strengths, not molar concentrations.'}
    (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    draw_plots(out,rows,curves,qs,summary,trajectories,h)
    return summary


def draw_plots(out,rows,curves,qs,summary,trajs,h):
    plt.rcParams.update({'figure.dpi':125,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    names=list(curves)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for cond in names:
        r=[row for row in rows if row['condition']==cond and row['replica']==1]
        ax.plot([z['t_reduced'] for z in r],[z['largest_cluster_chains'] for z in r],label=cond)
    ax.set(xlabel='Reduced time (not ns)',ylabel='Chains in largest connected cluster',ylim=(0,NCHAINS+.5),title='Bundle retention — independent replica 1');ax.legend();fig.tight_layout();fig.savefig(out/'01_cluster_retention.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for cond in names:
        matrix=curves[cond]
        y=matrix.mean(axis=0)
        ax.loglog(qs,y,label=cond)
        if len(matrix)>1:
            sd=matrix.std(axis=0,ddof=1)
            ax.fill_between(qs,np.maximum(y-sd,1e-5),y+sd,alpha=.12)
    ax.set(xlabel='q (reduced inverse length; not nm^-1)',ylabel='I(q) / total bead weight (arbitrary)',title='Synthetic Debye scattering of simulated aggregates');ax.legend();fig.tight_layout();fig.savefig(out/'02_saxs_synthetic.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for cond in names:
        matrix=curves[cond]
        norm=matrix/matrix[:,0,None]
        ax.loglog(qs,norm.mean(axis=0),label=cond)
    ax.set(xlabel='q (reduced inverse length)',ylabel='I(q)/I(q_min), shape-normalized',title='Scattering shape (relative only)');ax.legend();fig.tight_layout();fig.savefig(out/'03_saxs_normalized.png');plt.close(fig)
    fig,axs=plt.subplots(1,len(names),figsize=(11,4),sharex=True,sharey=True)
    for ax,cond in zip(axs,names):
        xyz=trajs[(cond,0)][-1]
        for chain in range(NCHAINS):
            block=xyz[chain*len(h):(chain+1)*len(h)]
            ax.plot(block[:,0],block[:,2],alpha=.6,lw=1.1)
        ax.set_title(cond);ax.set_xlabel('x (reduced)')
    axs[0].set_ylabel('z (reduced)');fig.suptitle('Final chain projections (not calibrated structures)');fig.tight_layout();fig.savefig(out/'04_snapshots.png');plt.close(fig)
    return


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps',type=int,default=3500)
    parser.add_argument('--replicas',type=int,default=3)
    parser.add_argument('--stride',type=int,default=100)
    parser.add_argument('--seed',type=int,default=20261029)
    parser.add_argument('--out',type=Path,default=Path('results'))
    args=parser.parse_args()
    if args.steps<args.stride or args.steps%args.stride:parser.error('steps must be divisible by stride')
    if args.replicas<1:parser.error('replicas must be >=1')
    summary=run(args.out,args.steps,args.replicas,args.stride,args.seed)
    print('\nCondition summary (not experiment-validated):')
    for row in summary:
        print(f"{row['condition']:>4}: largest cluster {row['largest_cluster_chains_mean']:.2f} +/- {row['largest_cluster_chains_replica_sd']} chains")


if __name__=='__main__':main()
