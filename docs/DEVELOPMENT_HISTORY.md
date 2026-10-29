# Current update: Version 5 (29 October 2026)

**Start with [`v5_saxs_identifiability/README.md`](../v5_saxs_identifiability/README.md).**
Version 5 adds independent Debye-vs-sphere validation, exact SAXS structural
non-identifiability, 60 synthetic-noise inversion trials, and a retrospective audit
of saved Version 4 MD trajectories. Version 5 does **not** run new MD, execute the
original Shape2SAS software, or fit a matched real experimental SAXS dataset.

---

# UPDATE: v2 sequence-informed β-casein + implicit-calcium proxy (executed 2026-10-29)

**New project stage:** `v2_sequence_calcium/` — includes 9 executed trajectories (8 chains, 21 beads per chain), real bovine β-casein A2 sequence and annotated phosphoserine sites, localized effective bridge interactions, 5 passing numerical/physics tests, SAXS forward calculations, plots, and experiment-download instructions. See `v2_sequence_calcium/README.md` and `v2_sequence_calcium/START_HERE_FA.md`. This v2 is **not an experimentally validated casein model** and has no explicit calcium ions. Original v1 results below are preserved.

---

# CaCasein-MD-SAXS | casein-inspired, coarse-grained MD to synthetic scattering

**Stage: reproducible generic-polymer pilot, NOT a molecular model of casein.**

This small, CPU-only project explores whether different inter-chain cohesive interactions lead to distinguishable clusters, and whether the resulting structures can be distinguished via a simulated small-angle-scattering signal. It connects a **real stochastic molecular-dynamics trajectory** (Langevin BAOAB integrator) to Debye forward scattering on the simulated coordinates, without claiming experimentally verified SAXS.

## Research question / testable hypothesis

When all simulations start from the same loose bundle of 12 six-bead chains, will weaker nonbonded cohesion cause the bundle to dissociate more than stronger cohesion? Can *in-silico* scattering distinguish the resulting cluster populations?

The pilot tests a **generic interaction-strength proxy** (`epsilon_interchain` in units of kBT), **not** explicit calcium concentration, pH, protein sequence or intermolecular biochemistry. The initialization is intentionally clustered. It measures short-time *bundle retention/dissociation*, not calibrated equilibrium association free energies.

## Reproducible quick start (Linux, CPU)

```bash
unzip CaCasein-MD-SAXS.zip
cd CaCasein-MD-SAXS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pytest -q tests
python cacasein_md_saxs.py --steps 12000 --replicas 3 --stride 200 --out results_reproduced
```

To run a quick smoke simulation:

```bash
python cacasein_md_saxs.py --steps 2000 --replicas 1 --stride 100 --out quickcheck
```

`results/` already contains the executed reference run (fixed base seed 20261029); rerunning should produce the same output for the same software versions and hardware-dependent floating-point behavior. No GPU and no OpenMM/GROMACS installation are needed.

## Model: 12 chains × 6 beads

- Mass = 1, Lennard-Jones distance sigma = 1, kBT = 1, periodic box length L = 30 (all **reduced units**, not physical nm or ns).
- Bonded neighbors: harmonic springs `U = 0.5 * k * (r - 1)^2`, `k = 100`.
- Intrachain non-neighbors: repulsive WCA potential (cutoff `2^(1/6) sigma`).
- Interchain pairs: *truncated-and-potential-shifted* Lennard-Jones, cutoff `2.5 sigma`, at three values `epsilon/kBT = 0.10, 0.65, 1.50`.
- Stochastic molecular dynamics with **BAOAB Langevin** integration, `dt = 0.004` reduced time units, friction `gamma = 1`, T = 1, 12,000 steps = 48 reduced time units. Three independent seeds per condition.
- Contact-based aggregation: an edge between two chains if any between-chain bead pair is within `1.5 sigma` under minimum-image periodic boundaries. Connected components define clusters. The cluster definition is arbitrary but explicit.
- Late-snapshot *analysis window* is 55–100% of the total trajectory. No equilibrium or decorrelation proof is claimed.

## Forward scattering

The isotropic Debye expression is evaluated on bead coordinates for each connected cluster separately:

\[ I(q)=\sum_i b_i^2 + 2\sum_{i<j} b_i b_j \frac{\sin(qr_{ij})}{qr_{ij}}. \]

All bead scattering weights are set to 1; intensities are averaged over late saved snapshots and replicas and normalized by the total number of beads. To mimic *dilute solution* measurements, between-cluster interference is omitted. Thus at exactly `q=0`, `I(0)/bead = sum_k n_k^2 / N`, where `n_k` is the number of beads in cluster k. This is a **synthetic SAS forward calculation in arbitrary units**, not an instrument-calibrated SAXS intensity curve. A full physical SAXS calculation would need solvent-excess electron contrast, realistic molecular form factors, sample concentration, resolution/background handling, and appropriate interparticle structure factors.

The within-cluster pair-distance histogram `P(r)` is a discretized **pair-count distribution**, not the inversion of real experimental SAXS data.

## Baseline results (executed, seed 20261029)

Three independent trajectories per condition; reported means use the last 45% of frames. The uncertainty figures are **standard deviations of independent-replica means** (not independent-frame errors):

| Interchain attraction | Mean largest cluster, chains | Replica SD | Mean cluster count | I(0)/bead, model units |
| --- | ---: | ---: | ---: | ---: |
| Weak: 0.10 kBT | 2.19 | 1.03 | 10.19 | 8.74 |
| Medium: 0.65 kBT | 8.77 | 3.21 | 3.68 | 45.40 |
| Strong: 1.50 kBT | 12.00 | 0.00 | 1.00 | 72.00 |

**Interpretation:** under these particular short-time conditions, stronger LJ cohesion preserves larger connected aggregates, generating greater *model* low-q scattering. The medium regime varies substantially across only three seeds. These trajectories do **not** demonstrate chemical calcium mediation, real casein self-assembly, SAXS validation, or statistical convergence to equilibrium.

## Files

- `cacasein_md_saxs.py` — complete dynamics, contact graph, Debye forward SAXS, diagnostics, figures.
- `analyze_replicas.py` — replicate-level summary, including between-replica SD.
- `tests/test_physics.py` — four tests: Debye q=0 + translation, finite-difference force/energy, PBC contacts + dilute sum, integrator smoke check.
- `results/*_trajectory.npz` — nine raw executed trajectories (time, XYZ, potential, box, epsilon).
- `results/frame_metrics.csv` — metrics across all saved frames.
- `results/summary.csv`, `summary.json` — group-level metrics.
- `results/per_replica_summary.csv` and `independent_replica_uncertainty.csv` — independent replicate analysis.
- `results/saxs_curves.csv`, `pair_distance.csv` — scattering curves and P(r) histograms.
- `results/01_rg_time.png` — radius of gyration of largest cluster vs. time (replica 1).
- `results/02_clusters.png` — chain count in largest cluster vs. time (replica 1).
- `results/03_debye_saxs.png` — synthetic Debye SAS over q.
- `results/04_pair_distance.png` — pair-distance distribution.
- `results/05_snapshot.png` — final coordinate snapshots (replica 1).
- `results/metadata.json` — parameter and model limitations.

## Checks performed

Four tests passed in the build environment. The fully run source produced all 9 finite trajectories and the figures/CSV files. The baseline ran with low memory on CPU. For an actual scientific study, add step-size sensitivity, long-window decorrelation analysis, convergence across >3 replicates, system-size sensitivity, cluster-cutoff sensitivity, initial-condition sensitivity (assembled vs separated), SAXS finite-size/PBC treatment, calibrated contrast and experimental validation.

## Relationship to the Copenhagen group and planned extensions

The conceptual link is the group's research on structural scattering and simulation-to-SAS analysis. Andreas Haahr Larsen and Jacob Judas Kain Kirkensgaard coauthored **Shape2SAS** (Journal of Applied Crystallography, 2023), which describes forward SAS and pair-distance distributions: https://doi.org/10.1107/S1600576723005848 . A general review of the physical basis of coarse-grained Langevin MD is https://doi.org/10.3390/biom11091347 . Neither publication is claimed as a validation of *this particular model*.

**Stage 2 (not done yet):** develop sequence-/structure-dependent alpha-, beta-, or kappa-casein coarse-graining with validated protein interactions, appropriate calcium-binding/charge model, realistic salt/pH effects, realistic solvent contrast; compare with published measured SAXS; use held-out experimental datasets and uncertainty propagation.

**Suggested truthful CV line:** "Built and tested a reproducible CPU-only coarse-grained Langevin MD pilot (12 bead-spring chains, 3 interaction strengths × 3 seeds) and computed Debye small-angle-scattering profiles from simulated cluster trajectories; model is not yet chemically calibrated to casein."

## V4 robustness audit (executed 2026-10-29)

See [`v4_validation/RESULTS_V4_EN.md`](../v4_validation/RESULTS_V4_EN.md) and [`v4_validation/START_HERE_FA.md`](../v4_validation/START_HERE_FA.md). This version runs 15 short independent trajectories and 6 longer stress-test trajectories, exposing a major ion-initialization bias in the earlier version. No Martini 3 simulation or experimental validation is claimed.
