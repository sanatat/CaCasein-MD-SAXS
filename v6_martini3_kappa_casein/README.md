# V6: κ-casein Martini 3 and Martini3-IDP pilot

**Research status (30 October 2026): actual CPU GROMACS runs completed by Sana Tat on Linux Mint; preliminary technical study, not an equilibrated or experimentally validated casein-micelle model.** Earlier V1–V5 are different synthetic/coarse-grained studies and are unchanged.

## Scientific question

Can a protein-specific coarse-grained model be used to examine the sensitivity of the flexible C-terminal **κ-casein** segment to conformational sampling and parameterization? The study tests the unmodified **64-residue C-terminal segment** (mature κ-casein 106–169; UniProt P02668, precursor 127–190), *not* a whole micelle. Starting coordinates are from a low-confidence AlphaFold prediction, not an experimental conformational ensemble. Phosphorylation, glycosylation, explicit calcium-binding chemistry and protein aggregation are **not** modeled.

## Runs verified from author's GROMACS terminal logs

| Model | What was actually executed | Outcome |
| --- | --- | --- |
| Standard Martini 3 | 10 ns GROMACS 2023.3 production (10 fs step), following minimization and equilibration; 2001 frames | Backbone Rg averaged **2.0060 nm** over full 10 ns; **1.7311 ± 0.1528 nm** over 8–10 ns (framewise SD) |
| Martini3-IDP | Independent Martinize2 topology (`-ff martini3IDP`, `-id-regions 1:64`, `-idr-tune`); minimization | Minimization converged in 2236 steps; transient LINCS warnings at steps 33–35 |
| Martini3-IDP | 500 ps NVT equilibration (2 fs) and 100 ps NVT production probe (5 fs) | No matching production/equilibration LINCS warnings; mean temperature **299.911 K** during the 100 ps probe |

**These numbers are transcribed from the author's actual local terminal outputs, not recalculated from trajectories in this repository.** The 10 ns standard run, IDP minimization, equilibration and production were successfully executed locally; **raw trajectories/checkpoints have not been archived here**.

### Standard Martini 3 time-block Rg (backbone beads)

| Time (ns) | Mean Rg (nm) | Framewise SD (nm) |
| --- | ---: | ---: |
| 0–2 | 2.3386 | 0.0904 |
| 2–4 | 2.2082 | 0.1045 |
| 4–6 | 2.0204 | 0.0927 |
| 6–8 | 1.7320 | 0.1635 |
| 8–10 | 1.7311 | 0.1528 |

The apparent contraction **is not evidence for the native κ-casein brush size**. The standard and IDP simulations are very different lengths; **no comparative conclusion about either force field can currently be supported**. The final time-block similarity alone is not proof of thermodynamic convergence. Framewise standard deviations do not account for temporal autocorrelation.

## What is in this directory

- `scripts/reference.py` and `scripts/prepare_structure.py`: documented P02668 sequence verification and AlphaFold-derived 64-residue input extraction.
- `run_martini3_pilot.sh`: executable *short* standard Martini 3 smoke-test pipeline using Martinize2, Insane and GROMACS. Modes: `preflight`, `run`, `resume`.
- `run_idp_pilot.sh`: independent Martini3-IDP preparation and guarded 500 ps equilibration/100 ps probe. It deliberately requires manual review of any minimization LINCS warnings before proceeding.
- `scripts/wire_topology.py`: assembles proper protein/solvent/ion includes and molecule counts without overwriting the protein topology.
- `scripts/analyze_martini.py`: optional backbone geometry and **synthetic point-scatterer proxy**, not SAXS measurements or calibrated SAXS prediction.
- `mdp/*.mdp`: example inputs. Reaction-field setting `epsilon-rf = 0` corrects the earlier erroneous `inf` value. Standard smoke test is 200 ps; the 10 ns local run was separately extended by `gmx convert-tpr`.
- `tests/`: Python tests for reference sequence, structure extraction and topology assembly; **these do not test physics, validate force fields or substitute for a GROMACS run**.

## Reproduce the technical tests (Linux / CPU)

Required: Python 3.10+, `gmx` (local author run: GROMACS 2023.3), Martinize2 (local: Vermouth 0.15.0), Insane, `curl`, `unzip`, network for official Martini and AlphaFold input downloads. The **exact local Conda/Python/pip dependency lockfile was not archived**, so this is a runnable reconstruction, not a byte-identical reproducibility capsule.

```bash
cd v6_martini3_kappa_casein
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
# Install Martinize2/Vermouth and Insane separately from their official upstream projects.
python -m pytest -q tests
bash run_martini3_pilot.sh preflight
bash run_martini3_pilot.sh run
```

To prepare the IDP counterpart in another folder (only after a successful standard `run` downloaded the required force-field files):

```bash
bash run_idp_pilot.sh prepare
# REVIEW work_idp_v3/min.log and constrained geometry, particularly any LINCS warnings.
bash run_idp_pilot.sh continue-after-review
```

The IDP script pauses after minimization on purpose because the observed local minimization showed three transient LINCS warnings; use the second command only after review. This script does *not* execute a 10 ns Martini3-IDP run. Do **not** use unreviewed outputs as protein stability evidence.

### How the author extended the standard smoke test to 10 ns

From within `work/` after the 200 ps short run, the author continued from checkpoint:

```bash
gmx convert-tpr -s prod.tpr -extend 1800 -o prod_2ns.tpr
gmx mdrun -s prod_2ns.tpr -cpi prod.cpt -deffnm prod -nt 2 -append
gmx convert-tpr -s prod_2ns.tpr -extend 8000 -o prod_10ns.tpr
gmx mdrun -s prod_10ns.tpr -cpi prod.cpt -deffnm prod -nt 2 -append
```

The 10 ns data were analyzed with `gmx gyrate` on the 64 `BB` beads after PBC centering. For quantitative interpretation, a full molecule-unwrapping / sensitivity audit and independent seeds remain necessary.

## Limitations and integrity

- Single 10 ns standard trajectory and much shorter Martini3-IDP probe; **no ensemble convergence, uncertainty calibration, experimental SAXS/SANS agreement or force-field ranking**.
- No micelle environment, κ-casein post-translational modifications, explicit calcium-phosphoserine interaction, or scattering contrast/hydration model.
- Model starts from low-confidence predicted coordinates, which are **not a validated unfolded conformational ensemble**.
- Some minimization LINCS warnings occurred for the IDP model; their absence in later MD stages does not erase them.
- Topology differences and integration timesteps complicate comparisons; numerical stability is not validation of protein physics.

References: [UniProt P02668](https://www.uniprot.org/uniprotkb/P02668/entry); [Martini protein tutorials](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/); [Martini3-IDP tutorial](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut3.html).
