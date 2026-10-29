#!/usr/bin/env python3
"""V5: genuine executed numerical study of SAXS identifiability + V4 MD audit.

No casein validation, no Martini 3 run, no direct Shape2SAS benchmarking.
All synthetic geometries and model coordinates are in uncalibrated length units.
"""
import argparse
import csv
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import t as tdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from saxs_core import (point_debye,sphere_form_factor,random_sphere_points,
                       weighted_rg,guinier_rg,homometric_pair,pair_distribution,
                       posterior_from_saxs)

HERE=Path(__file__).resolve().parent
LABELS={
 'no_ca':'0 Ca',
 'ca8_inert_bulk':'8 inert Ca, bulk',
 'ca8_coord_bulk':'8 binding Ca, bulk',
 'ca8_coord_near':'8 binding Ca, near',
 'ca16_coord_bulk':'16 binding Ca, bulk'
}

def csv_save(path, rows):
    if not rows: return
    with open(path,'w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        wr.writeheader();wr.writerows(rows)


def csv_read(path):
    with open(path,newline='',encoding='utf-8') as f:return list(csv.DictReader(f))


def v4_audit(v4dir,out):
    """Re-analysis of existing V4 runs, NOT additional MD."""
    if not (v4dir/'results'/'replica_summary.csv').exists():
        return {'status':'skipped','reason':'V4 result CSVs not found'}
    rows=csv_read(v4dir/'results'/'replica_summary.csv')
    per_cond=defaultdict(list)
    for r in rows:
        per_cond[r['condition']].append(r)
    outrows=[]
    for cond,parts in per_cond.items():
        early=np.array([float(r['calcium_contact_count_early_mean']) for r in parts])
        late=np.array([float(r['calcium_contact_count_late_mean']) for r in parts])
        d=late-early
        n=len(d); sd=float(d.std(ddof=1)) if n>1 else float('nan')
        hw=float(tdist.ppf(.975,n-1)*sd/(n**.5)) if n>1 else float('nan')
        outrows.append(dict(condition=cond,n_independent_replicas=n,
                           early_contact_mean=float(early.mean()),
                           late_contact_mean=float(late.mean()),
                           late_minus_early_mean=float(d.mean()),
                           late_minus_early_95pct_t_halfwidth=hw,
                           late_replica_sd=float(late.std(ddof=1)) if n>1 else float('nan')))
    csv_save(out/'v4_replica_audit.csv',outrows)
    curves=csv_read(v4dir/'results'/'synthetic_saxs_replicas.csv')
    byrep=defaultdict(list)
    for r in curves:
        byrep[(r['condition'],int(r['replica']))].append((float(r['q_reduced']),float(r['I_normalized'])))
    for k in byrep:byrep[k].sort(key=lambda x:x[0])
    keys=set(k[1] for k in byrep if k[0]=='ca8_coord_bulk')
    paired=[]
    for rep in sorted(keys):
        k1=('ca8_coord_bulk',rep);k2=('ca8_coord_near',rep)
        if k1 not in byrep or k2 not in byrep:continue
        q1=np.array([e[0] for e in byrep[k1]]);q2=np.array([e[0] for e in byrep[k2]])
        if not np.allclose(q1,q2): raise ValueError('V4 grids differ')
        a=np.array([e[1] for e in byrep[k1]])
        b=np.array([e[1] for e in byrep[k2]])
        paired.append(dict(replica=rep,matched_seed_rms_log_normalized_saxs=float(np.sqrt(np.mean(np.log(a/b)**2)))))
    csv_save(out/'v4_saxs_matched_replicas.csv',paired)
    fig,axs=plt.subplots(1,2,figsize=(11,4.4))
    for cond in ['ca8_coord_bulk','ca8_coord_near']:
        data=[r for r in outrows if r['condition']==cond][0]
        v=data['late_contact_mean'];err=data['late_replica_sd']
        ix=0 if cond=='ca8_coord_bulk' else 1
        axs[0].bar(ix,v,yerr=err,capsize=5,color='C'+str(ix))
    axs[0].set_xticks([0,1],['Bulk-start','Near-site-start'])
    axs[0].set_ylabel('Late mean Ca contacts ± between-replica SD')
    axs[0].set_title('V4 initiation dependence')
    for i,cond in enumerate(['ca8_coord_bulk','ca8_coord_near']):
        allcurves=np.stack([np.array([e[1] for e in vals]) for (c,r),vals in byrep.items() if c==cond])
        q=np.array([e[0] for e in byrep[(cond,1)]])
        mu=allcurves.mean(axis=0);sd=allcurves.std(axis=0,ddof=1)
        axs[1].semilogy(q,mu,label=('Bulk' if i==0 else 'Near'),color='C'+str(i))
        axs[1].fill_between(q,np.maximum(mu-sd,1e-6),mu+sd,alpha=.14,color='C'+str(i))
    axs[1].set(xlabel='q [inverse reduced length]',ylabel='V4 synthetic P(q)/P(q_min)',title='Protein-only forward scattering')
    axs[1].legend()
    fig.tight_layout();fig.savefig(out/'04_v4_reanalysis.png',dpi=165);plt.close(fig)
    return {'status':'executed_from_existing_v4_trajectories_and_CSVs',
            'n_trajectories':len(rows),'bulk_vs_near_matched_saxs_log_rms_mean':float(np.mean([r['matched_seed_rms_log_normalized_saxs'] for r in paired])),
            'contact_summary':outrows,'matched_saxs':paired}


def run(out,v4dir,seed=20261029,noise_reps=60):
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'figure.dpi':145,'axes.spines.top':False,
                         'axes.spines.right':False,'font.size':10})
    # Benchmark Debye using mathematically INDEPENDENT continuous-density form factor.
    qsp=np.linspace(.004,.20,90)
    radius=12.
    analytic=sphere_form_factor(qsp,radius)
    samples=[];rg_vals=[];guinier=[]
    for i in range(3):
        points=random_sphere_points(850,radius,seed+100*i)
        calc=point_debye(points,qsp)
        samples.append(calc)
        rg_vals.append(weighted_rg(points))
        guinier.append(guinier_rg(qsp,calc, qmax=.075))
    samples=np.asarray(samples)
    fit=samples.mean(axis=0)
    rms=float(np.sqrt(np.mean((fit-analytic)**2)))
    max_abs=float(np.max(np.abs(fit-analytic)))
    abs_at_high=float(abs(fit[-1]-analytic[-1]))
    sphere_rows=[{'q_inverse_reduced':float(q),'P_analytic':float(a),
                  'P_debye_average':float(b),'P_debye_replica_sd':float(s)}
                 for q,a,b,s in zip(qsp,analytic,fit,samples.std(axis=0,ddof=1))]
    csv_save(out/'sphere_forward_benchmark.csv',sphere_rows)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    axs[0].plot(qsp,analytic,'k--',label='Continuous sphere analytic')
    axs[0].plot(qsp,fit,label='Point Debye (3 × 850 beads)')
    axs[0].fill_between(qsp,fit-samples.std(axis=0,ddof=1),fit+samples.std(axis=0,ddof=1),alpha=.18)
    axs[0].set(xlabel='q [inverse reduced length]',ylabel='P(q)',title='Independent sphere benchmark')
    axs[0].legend(fontsize=9)
    axs[1].plot(qsp,fit-analytic,label='Debye - analytic')
    axs[1].axhline(0,c='k',lw=.8,alpha=.5)
    axs[1].set(xlabel='q [inverse reduced length]',ylabel='Absolute P(q) error',title='Finite-point representation error')
    fig.tight_layout();fig.savefig(out/'01_sphere_benchmark.png',dpi=165);plt.close(fig)

    # Exact, mathematically proven structural non-identifiability (homometry).
    A,B=homometric_pair()
    C=np.zeros((6,3));C[:,0]=[0,2,4,8,12,17]
    dA=np.sort(pdist(A));dB=np.sort(pdist(B))
    if not np.array_equal(dA,dB): raise AssertionError('Homometry failed')
    q=np.geomspace(.008,2.7,180)
    curves={'A':point_debye(A,q),'B':point_debye(B,q),'C':point_debye(C,q)}
    equality=float(np.max(np.abs(curves['A']-curves['B'])))
    # Exact homometry also implies identical Rg and P(r), but not geometries.
    centers,hA,edges=pair_distribution(A,bins=18,rmax=18.0)
    _,hB,_=pair_distribution(B,bins=18,rmax=18.0)
    csv_save(out/'homometric_scattering.csv',[{'q_inverse_reduced':float(qq),
        **{f'P_{k}':float(v[i]) for k,v in curves.items()}}
        for i,qq in enumerate(q)])
    csv_save(out/'homometric_pairdist.csv',[{'r_center_reduced':float(x),'pair_count_A':float(y),'pair_count_B':float(z)} for x,y,z in zip(centers,hA,hB)])
    fig,axs=plt.subplots(1,3,figsize=(13.2,4.1))
    axs[0].scatter(A[:,0],np.full(len(A),.1),s=80,label='Geometry A')
    axs[0].scatter(B[:,0],np.full(len(B),-.1),s=80,label='Geometry B',marker='x')
    axs[0].set(xlabel='Position [reduced]',yticks=[],ylim=(-.4,.4),title='Different structures')
    axs[0].legend(fontsize=8)
    axs[1].step(centers,hA,where='mid',label='A')
    axs[1].plot(centers,hB,'x',label='B')
    axs[1].set(xlabel='Pair separation [reduced]',ylabel='Number of pairs',title='Exactly the same P(r)')
    axs[1].legend()
    for key in curves: axs[2].loglog(q,curves[key],label=key)
    axs[2].set(xlabel='q [inverse reduced]',ylabel='P(q)',title='A and B exactly identical')
    axs[2].legend()
    fig.tight_layout();fig.savefig(out/'02_exact_homometry.png',dpi=165);plt.close(fig)

    # Synthetic inverse experiment: Gauss errors, no nuisance calibration factors.
    # Evaluate 60 independently generated observation vectors using a known A truth.
    sigma=np.full(len(q),.025)
    rng=np.random.default_rng(seed+9981)
    model_keys=list(curves)
    rows=[]
    for rep in range(noise_reps):
        y=curves['A']+rng.normal(0,sigma)
        pos,chis=posterior_from_saxs(y,sigma,curves)
        pos_with_prior,_=posterior_from_saxs(y,sigma,curves,priors={'A':.8,'B':.1,'C':.1})
        rows.append({'replicate':rep+1,**{f'posterior_{k}':float(pos[k]) for k in model_keys},
                     **{f'chi2_{k}':float(chis[k]) for k in model_keys},
                     'posterior_A_with_external_prior':float(pos_with_prior['A']),
                     'posterior_B_with_external_prior':float(pos_with_prior['B'])})
        if rep==0:example=y.copy()
    csv_save(out/'inverse_ensemble.csv',rows)
    csv_save(out/'inverse_example.csv',[dict(q_inverse_reduced=float(qq),
           intensity_synthetic=float(yy),sigma_absolute=float(ss),
           model_A=float(a),model_B=float(b),model_C=float(c))
           for qq,yy,ss,a,b,c in zip(q,example,sigma,*[curves[k] for k in model_keys])])
    fig,axs=plt.subplots(1,2,figsize=(11.2,4.4))
    axs[0].errorbar(q[::5],example[::5],yerr=sigma[::5],fmt='.',markersize=3,
                    label='Noisy synthetic observations (A)',color='k',alpha=.7)
    for k in model_keys:axs[0].plot(q,curves[k],label='Model '+k,lw=1.6)
    axs[0].set(xscale='log',xlabel='q [inverse reduced]',ylabel='P(q)',title='Identical forward models A and B')
    axs[0].legend(fontsize=8)
    axs[1].hist([r['posterior_A'] for r in rows],bins=np.linspace(0,1,21),label='Posterior A')
    axs[1].hist([r['posterior_B'] for r in rows],bins=np.linspace(0,1,21),histtype='step',lw=2,label='Posterior B')
    axs[1].set(xlabel='Posterior probability across synthetic repeats',ylabel='Count',title='Uniform model priors; n='+str(noise_reps))
    axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'03_inverse_uncertainty.png',dpi=165);plt.close(fig)
    # No actual experimental input; audit whether real data can directly validate the toy model.
    v4=v4_audit(v4dir,out)
    results={
       'status':'executed CPU analysis on synthetic geometries and existing V4 trajectories',
       'date':'2026-10-29',
       'shape2sas_direct_run':False,
       'experimental_saxs_loaded':False,
       'new_md_simulations':0,
       'sphere':{'radius_reduced':radius,'number_point_replicas':3,'points_per_replica':850,
                 'P_rms_absolute_error_vs_analytic':rms,'P_max_absolute_error_vs_analytic':max_abs,
                 'P_high_q_absolute_error':abs_at_high,
                 'direct_Rg_point_means_reduced':rg_vals,
                 'guinier_Rg_reduced': [rr[0] for rr in guinier],
                 'theoretical_sphere_Rg_reduced':float(np.sqrt(3/5)*radius)},
       'homometric':{'sets_A':A[:,0].tolist(),'sets_B':B[:,0].tolist(),
                     'number_identical_pair_distances':len(dA),
                     'max_absolute_SAXS_difference_all_q':equality,
                     'Rg_A_reduced':weighted_rg(A),'Rg_B_reduced':weighted_rg(B),
                     'distance_histograms_identical':bool(np.array_equal(hA,hB))},
       'inverse':{'synthetic_repetitions':noise_reps,'pointwise_sigma_absolute':float(sigma[0]),
                  'uniform_prior_mean_probability_A':float(np.mean([r['posterior_A'] for r in rows])),
                  'uniform_prior_mean_probability_B':float(np.mean([r['posterior_B'] for r in rows])),
                  'uniform_prior_mean_probability_C':float(np.mean([r['posterior_C'] for r in rows])),
                  'example_chi2_A':rows[0]['chi2_A'],'example_chi2_B':rows[0]['chi2_B'],
                  'external_prior_A_0.8_B_0.1_C_0.1_mean_posterior_A':float(np.mean([r['posterior_A_with_external_prior'] for r in rows]))},
       'v4_reanalysis':v4,
       'units':'all simulated q and distances are uncalibrated reduced units, not nm-1 or nm',
       'critical_limits':['No Martini 3','No explicit-solvent validated calcium force field',
                          'No directly matched experimental SAXS','No direct Shape2SAS installation in offline environment',
                          'Exact homometric sets are 1D mathematical examples, not casein conformations',
                          'V4 frame correlations and under-sampling make confidence intervals exploratory']
    }
    (out/'results_summary.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in results.items() if k in ['sphere','homometric','inverse']},indent=2))
    print('Saved',out)
    return results

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=HERE/'results')
    parser.add_argument('--v4',type=Path,default=HERE.parent/'v4_validation')
    parser.add_argument('--seed',type=int,default=20261029)
    parser.add_argument('--noise-reps',type=int,default=60)
    args=parser.parse_args()
    run(args.out,args.v4,args.seed,args.noise_reps)
