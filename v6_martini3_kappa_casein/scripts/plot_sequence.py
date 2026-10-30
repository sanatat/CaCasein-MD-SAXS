from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reference import MATURE

def main(out):
    x=np.arange(1,len(MATURE)+1)
    charge=np.array([-1 if aa in 'DE' else 1 if aa in 'KR' else 0 for aa in MATURE])
    hydro=np.array([1 if aa in 'AILMFWVYP' else 0 for aa in MATURE])
    # Descriptive local sequence proxies, not a prediction of micelle morphology.
    kernel=np.ones(15)/15
    smoothcharge=np.convolve(charge,kernel,'same')
    fig,ax=plt.subplots(figsize=(10,3.8))
    ax.axvspan(106,169,color='0.9',label='C-terminal tail (106-169)')
    ax.plot(x,smoothcharge,label='15-residue local charge proxy',lw=2)
    ax.scatter(x[charge>0],charge[charge>0]*0.11,s=7,label='Lys/Arg')
    ax.scatter(x[charge<0],charge[charge<0]*0.11,s=7,label='Asp/Glu')
    ax.axvline(105.5,color='0.3',ls=':',label='Chymosin F105|M106')
    ax.axvline(127,color='0.6',ls='--',label='Reference phosphoserine S127 (not simulated)')
    ax.set(xlabel='Position in mature bovine kappa-casein',ylabel='Approximate local side-chain charge / residue',xlim=(1,169))
    ax.legend(fontsize=7,loc='lower left',ncol=2)
    fig.tight_layout();Path(out).parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=180);plt.close(fig)

if __name__=='__main__':
    import sys
    main(sys.argv[1] if len(sys.argv)>1 else 'results/sequence_annotation.png')
