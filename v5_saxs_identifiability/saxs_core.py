"""Independent numerical SAXS forward models and diagnostics.

This module deliberately keeps all distances dimensionless unless the caller
supplies a *calibrated* length scale. Shape2SAS-inspired Debye calculation;
not an import of (or replacement for) the Shape2SAS implementation.
"""
from __future__ import annotations
import numpy as np
from scipy.spatial.distance import pdist


def point_debye(xyz: np.ndarray, q: np.ndarray, weights=None) -> np.ndarray:
    """P(q)=sum_ij w_i w_j sinc(q*r_ij/pi) / (sum_i w_i)^2.

    Identical point scatterers; solution dilute and orientationally averaged.
    Diagonal self-scattering and both off-diagonal directions included.
    No atom form factors, hydration shell, contrast, structure factor or resolution.
    """
    xyz = np.asarray(xyz, dtype=float)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 1:
        raise ValueError('xyz must be an (N,3) array')
    if np.any(q < 0) or not np.all(np.isfinite(q)):
        raise ValueError('q must be nonnegative and finite')
    if not np.isfinite(xyz).all():
        raise ValueError('xyz must be finite')
    n = len(xyz)
    w = np.ones(n) if weights is None else np.asarray(weights, float)
    if w.shape != (n,) or np.any(~np.isfinite(w)) or abs(w.sum()) < 1e-12:
        raise ValueError('invalid scattering weights')
    r = pdist(xyz)
    iu, ju = np.triu_indices(n, k=1)
    pair_weights = w[iu]*w[ju]
    pref = np.sum(w*w)
    out = np.empty(len(q), dtype=float)
    # Blocked to keep memory bounded and usable on a 6-8 GiB machine.
    for start in range(0,len(q),16):
        qq = q[start:start+16]
        out[start:start+16] = (pref + 2.0*np.sum(pair_weights[:,None]*np.sinc((r[:,None]*qq[None,:])/np.pi),axis=0)) / (w.sum()**2)
    return out


def pair_distribution(xyz, bins=140, weights=None, rmax=None):
    """Weighted unordered-pair histogram for plotting; excludes r=0 self terms."""
    xyz=np.asarray(xyz,float)
    w=np.ones(len(xyz)) if weights is None else np.asarray(weights,float)
    iu,ju=np.triu_indices(len(xyz),1)
    r=pdist(xyz)
    high=float(r.max()*1.001) if rmax is None else float(rmax)
    hist,edges=np.histogram(r,bins=bins,range=(0,high),weights=w[iu]*w[ju])
    return .5*(edges[:-1]+edges[1:]), hist, edges


def weighted_rg(xyz,weights=None):
    xyz=np.asarray(xyz,float)
    w=np.ones(len(xyz)) if weights is None else np.asarray(weights,float)
    if np.any(w<0) or w.sum()<=0: raise ValueError('Rg requires nonnegative positive-mass weights')
    center=np.average(xyz,axis=0,weights=w)
    return float(np.sqrt(np.sum(w*np.sum((xyz-center)**2,axis=1))/w.sum()))


def guinier_rg(q, p, qmax, intensity_unc=None):
    """Linear fit log[P(q)] = a - (Rg^2/3)q^2 at low q. For positive weights."""
    q=np.asarray(q,float);p=np.asarray(p,float)
    select=(q>0)&(q<=qmax)&(p>0)
    if select.sum()<4: raise ValueError('Need >=4 positive q values')
    x=q[select]**2;y=np.log(p[select])
    if intensity_unc is None:
        slope,intercept=np.polyfit(x,y,1)
    else:
        sigma=np.asarray(intensity_unc,float)[select]
        slope,intercept=np.polyfit(x,y,1,w=p[select]/sigma)
    return float(np.sqrt(max(0.,-3*slope))),float(np.exp(intercept)),int(select.sum())


def sphere_form_factor(q, radius):
    """Analytical dilute solid-sphere P(q), normalized to P(0)=1."""
    x=np.asarray(q,float)*float(radius)
    amplitude=np.ones_like(x)
    mask=np.abs(x)>1e-5
    xx=x[mask]
    amplitude[mask]=3*(np.sin(xx)-xx*np.cos(xx))/xx**3
    if np.any(~mask):
        ss=x[~mask];amplitude[~mask]=1-ss**2/10+ss**4/280
    return amplitude**2


def random_sphere_points(n, radius, seed):
    rng=np.random.default_rng(seed)
    v=rng.normal(size=(n,3));v/=np.linalg.norm(v,axis=1)[:,None]
    return v*(rng.random(n)**(1/3)*radius)[:,None]


def homometric_pair():
    """Classic 1D noncongruent homometric point sets (same pair distances)."""
    a=np.array([0,1,4,10,12,17],float)
    b=np.array([0,1,8,11,13,17],float)
    x=np.zeros((6,3));y=np.zeros((6,3))
    x[:,0]=a;y[:,0]=b
    return x,y


def posterior_from_saxs(y, sigma, curves, priors=None):
    """Discrete Gaussian-noise evidence. No parameter fitting or look-elsewhere correction."""
    names=list(curves)
    y=np.asarray(y,float);sig=np.asarray(sigma,float)
    if np.any(sig<=0) or y.shape!=sig.shape: raise ValueError('bad Gaussian observations')
    prior=np.full(len(names),1/len(names)) if priors is None else np.array([priors[k] for k in names],float)
    if np.any(prior<=0) or not np.isclose(prior.sum(),1):raise ValueError('Priors must be positive and sum to 1')
    chi=np.array([np.sum(((y-np.asarray(curves[k]))/sig)**2) for k in names])
    logp=np.log(prior)-.5*chi
    prob=np.exp(logp-logp.max());prob/=prob.sum()
    return dict(zip(names,prob)),dict(zip(names,chi))
