# CaCasein-MD-SAXS V4 — ion initialization and sampling robustness

**Status (2026-10-29): Executed computational model audit, NOT an experimentally validated casein study.**

## Scientific question
How strongly are apparent calcium–phosphoserine contacts in a reduced-unit, sequence-informed, explicit-Ca *toy model* affected by where mobile ions are placed initially, and how much of the effect survives longer simulations?

This audit focuses on a methodological failure mode revealed by V3: mobile calcium-like beads were deliberately placed near phosphate-susceptible regions at the start, which can produce apparently strong occupancy even if equilibrium binding has not been demonstrated.

## Computational design and what was actually executed

- 4 × mature bovine β-casein A2 model chains, each represented by 21 reduced-unit beads of about 10 residues apiece (84 protein beads); 0, 8, or 16 Ca-like mobile beads.
- Simple harmonic bonds, repulsion, screened charge proxy, hydrophobic sequence-dependent attraction, and a phenomenological Gaussian Ca–phospho attraction, integrated with Langevin dynamics; **not Martini 3**.
- All five conditions are defined in `casein_v4.py`. Three independent seeds/replicas each; 4,000 integration steps per run and one frame every 100 steps: **15 simulation trajectories**, saved as `.npz` in `results/`.
- Additional stress test of the 8-ion, specific-attraction conditions with either *bulk* or *near-phosphate* ion starting positions: 3 replicas × 2 initializations, 12,000 integration steps each: **six additional long-check simulations**; time traces and split-third diagnostics in `long_time_results/`.
- Initial *protein* geometry was held constant across the initialization test, so this audit isolates calcium placement under fixed initial protein organization. The prepared protein bundle itself is artificial.

## Results: 4,000-step simulation cohort

| Condition | Bound Ca-like beads (replica mean ± SD) | Two-chain Ca bridges (replica mean ± SD) |
|---|---:|---:|
| No calcium | 0.000 ± 0.000 | 0.000 ± 0.000 |
| 8 ions, no specific attraction, uniform bulk start | 0.143 ± 0.247 | 0.000 ± 0.000 |
| 8 ions, specific attraction, uniform bulk start | 0.000 ± 0.000 | 0.000 ± 0.000 |
| 8 ions, specific attraction, **near-site start** | **3.302 ± 1.271** | **0.619 ± 1.031** |
| 16 ions, specific attraction, uniform bulk start | 0.286 ± 0.495 | 0.000 ± 0.000 |

Means are computed within the last half of the respective run, then averaged over **three independent replicas**. SD describes between-replica variation, not confidence that equilibrium has been attained. The occupancy metric counts *beads* near annotated phospho-containing coarse bins; it is not a binding constant.

## Long-run check: 12,000 steps, each condition in three independent replicas

**Mean Ca-bead contacts per frame by one-third of trajectory (across 3 replicas)**:

| Ion start | First third | Middle third | Final third |
|---|---:|---:|---:|
| Uniform bulk | 0.091 | 0.020 | **0.000** |
| Near phosphosites | 3.687 | 3.490 | **3.118** |

The outcomes remain markedly different even after tripling the simulation time. This suggests strong sensitivity to starting conditions and likely inadequate mixing in this low-density model, *not* quantitative biochemical affinity, successful calibration, or reliable equilibrium binding.

## SAXS result and interpretation

Synthetic Debye curves were calculated for **protein beads only**, treating separated clusters as independently scattering objects. They are displayed in reduced reciprocal length units, normalized by the smallest-q intensity, without solvent subtraction, polydispersity, instrumental resolution, calcium contrast, or a measured inter-aggregate structure factor. `04_saxs_distances.png` quantifies shape differences but **cannot** validate a casein micelle or compare to a real experiment without a justified physical mapping.

The scattering calculation passed zero-angle intensity, translation/rotation invariance, and bead-reordering tests. These check *mathematical implementation*, not scientific adequacy of the molecular model.

## Critical limitations

1. Unknown conversion of reduced lengths, time and ion numbers to nm, ns and molar Ca²⁺ concentration. Simulations cannot be compared to experimental q scales without an independently established scale.
2. Low effective encounter rates; calcium started near protein vs drawn uniformly within a much larger volume. The simulations did not establish ergodicity or adequate equilibration.
3. No explicit water/counterions; net electrostatic conditions cannot be mapped to aqueous buffered casein solutions.
4. Coarse binning merges several phosphoserines, so there is no resolution of calcium coordination geometries or binding thermodynamics.
5. Casein is intrinsically disordered and exists in an interacting multicomponent micellar environment; this is only four β-casein chains initialized in an idealized bundle.
6. SAXS is purely synthetic. No experimental scattering data have been fitted or used for validation.

## Why we did not claim a Martini 3 simulation

The present runtime has no `gmx` executable and does not contain a validated Martini 3 phosphorylated β-casein topology. Martinize2 documentation has explicitly cautioned that Martini 3 parameters for relevant post-translational modifications were not available in the version discussed. Arbitrarily substituting unphosphorylated serine or inventing phosphate–Ca²⁺ parameters would invalidate the result. See `MARTINI3_PATH.md`.

## Useful next experiment (after parameters and hardware setup)

- Use a validated Martini 3 topology with properly modeled phosphorylation and an explicit-solvent, neutralized ion solution.
- Equilibrate multiple, independently prepared beta-casein conformations and ion placements, demonstrate convergence of binding events/occupancy and test box size and calcium chemical potential effects.
- Compute SAXS from physically calibrated structures including solvent density contrast, compare to **composition-, pH-, temperature-, and sample-state-matched** scattering data; do not claim the enzymatic gelation time series is a matched calcium titration.

## Literature anchors

- Luo et al., *Food & Function* (2020): calcium–casein phosphopeptide binding and atomistic simulation, https://doi.org/10.1039/D0FO00844C . The paper studies the 1–25 peptide, **not** this four-chain model.
- Mekmene and Gaucheron, *Food Chemistry* (2011): experimental binding capacities and constants for casein, phosphoserine and other species, https://doi.org/10.1016/j.foodchem.2010.12.121 . Its measured constants **cannot** directly parametrize our dimensionless toy model.
- Martinize2/Vermouth, *eLife* (2025): PTM topology handling limitations, https://elifesciences.org/articles/90627v1 .
- Official Martini 3 protein tutorial: https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html .
- Unmatched public casein USAXS/SAXS dataset: https://zenodo.org/records/22044394 . It studies time-resolved **enzymatic gelation**, not a controlled Ca²⁺ titration.
