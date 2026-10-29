# CaCasein-MD-SAXS v3: explicit-Ca exploratory model

**Executed on 29 October 2026. Not a Martini simulation. Not validated against real SAXS.**

## Main change relative to v2

v2 treated calcium by changing inter-chain phosphate attraction strengths without calcium ions. v3 introduces **mobile divalent ion beads with charge +2** into a small sequence-informed bead/spring model. It retains simplified implicit solvent and phenomenological local attraction between Ca beads and phosphoserine-rich protein beads. This is a deliberately modest transition to explicit-ion *toy* MD, **not a calibrated Martini 3 force field, a physical molar Ca titration, or an all-atom simulation**.

The model uses 4 independent chains of bovine beta-casein A2. The 209-aa sequence of each chain is mapped to 21 segments; five phosphoserines are projected to two effective sites per chain. Approximate hydrophobicity and screened charge proxies are inherited from v2. Structural interpretation of a bin covering ~10 amino acids is limited.

## Four conditions

1. `no_ca`: protein only.
2. `ca8_no_coord`: 8 explicit charged Ca beads, **without** specific Ca–phosphate Gaussian binding (electrostatics and WCA still present).
3. `ca8_coord`: 8 explicit charged Ca beads, plus toy calcium–phosphate coordination potential (3 reduced kBT per effective site weight).
4. `ca16_coord`: 16 explicit charged Ca beads, same toy coordination potential.

These are **numbers of beads in a model box**, not concentrations; there are no physically justified volume units or counterions. Short runs start from a shared aligned protein bundle with calcium initially sampled near phosphoserine sites. This initial condition can bias binding and is explicitly *not* a demonstration of equilibrium.

## Simulation and analysis

- Pair potentials: harmonic protein backbone; WCA excluded volume; Gaussian interchain hydrophobic contacts; screened Yukawa charges; optional phenomenological local Ca–phosphoserine attraction.
- Integrator: stochastic BAOAB-like Langevin with velocity updates and force recomputation, reduced `dt=0.002`.
- Default outputs: 2,000 steps, 3 independent replicas per condition. Actual included run: 2,400 steps, 3 replicas, 4 conditions = **12 trajectories**.
- Protein-only SAXS: dilute/intracluster Debye scattering from saved snapshots; **shape only**, not absolute SAXS intensity and no nm scale or ionic/solvent contrast.
- Geometric measures: ion–phosphosite proximity, Ca ion adjacent to sites on >=2 different chains (a "bridge" proxy), interchain phosphate-site contacts, cluster retention, protein-bead pair-distance distribution (not a rigorous experimental P(r) inversion).
- Late-trajectory averages: last 50% of saved frames, not proven equilibrated. SD is across *independent replica means* (n=3), not across correlated frames.

## Required packages and running

Use Python 3.10+:

```bash
cd CaCasein-MD-SAXS/v3_explicit_ca
python3 -m venv .venv
source .venv/bin/activate
pip install -r ../requirements.txt
python -m pytest -q tests
python explicit_ca_md.py --steps 2400 --replicas 3 --stride 80 --out results_new
```

The parent `v2_sequence_calcium/casein_v2.py` is required for shared sequence and SAXS functions. Run from the full ZIP. No CUDA/GPU required.

## Files in results/

- `summary.csv`, `replica_summary.csv`, `metadata.json`: evaluated parameters and results.
- `*_trajectory.npz`: saved Cartesian coordinates of all protein and calcium beads, energy and reduced timestamps.
- `time_metrics.csv`: frame-level geometric metrics.
- `saxs_curve_replicas.csv`: forward synthetic normalized SAXS.
- `pair_distance_P_r_replicas.csv`: normalized intracluster protein-bead distance histograms.
- `01_synthetic_SAXS.png` to `05_snapshots.png`: figures.

## Experimental SAXS: found but not validated

Zenodo record **10.5281/zenodo.22044394** provides an August 2026 7z archive of real USAXS/SAXS measured during enzymatic gelation of reconstituted skim milk at pH 6.6 and 38°C. This is not beta-casein-only and not an isolated calcium-titration experiment. The archive could not be downloaded in the execution environment and **is not bundled**. Existing v2 helper can try downloading it on a locally internet-connected system:

```bash
python ../v2_sequence_calcium/experimental_data_helper.py fetch --out external/zenodo_casein_saxs_2026.7z
```

If 7z and corresponding Python library are installed, inspect the archive and select a numeric two-column profile (`q`, `I`). Be careful with units, blank/negative intensities, transmission background correction and condition matching.

An optional exploratory visual overlay exists in `compare_experiment.py`, but is **not an experiment-model fit**. It requires a physically justified length calibration (there is no default one in the code):

```bash
python compare_experiment.py --experiment /path/to/two_column_profile.txt \
  --model results/saxs_curve_replicas.csv --condition ca8_coord --q-unit nm^-1 \
  --length-nm-per-reduced LENGTH_VALUE_FROM_INDEPENDENT_CALIBRATION \
  --out EXPLORATORY_NOT_VALIDATION.png
```

Do **not** copy an arbitrary length value just to make curves align or report the result as validation.

## Scientific next steps before quantitative inference

1. Calibrate a protein structural ensemble, ideally multiple beta-casein conformers and measured solution data. Casein is largely intrinsically disordered: a single fixed protein fold is not adequate.
2. Build true Martini 3 topology (GROMACS + martinize2), explicitly **check availability and accuracy of phosphoserine and calcium coordination parameters**; do not mix Martini versions or improvise PTM beads.
3. Include realistic solvent, phosphate species, pH-related protonation, counterions, and ionic strength, calibrating protein-Ca binding against experiments.
4. Compare matched measurements across Ca concentrations: same casein composition, temperature, buffer, pH and background subtraction. Improve forward scattering with solvent electron-density and sample composition.
5. Demonstrate convergence with longer equilibration, multiple distinct initial configurations, finite-size tests and replicas.

For the correct method boundaries see `MARTINI3_UPGRADE_PLAN.md`.

## Bibliography

- Martini 3 force field: https://doi.org/10.1038/s41592-021-01098-3
- Calcium binding to beta-casein phosphopeptides: https://doi.org/10.1039/D0FO00844C
- Ca- and beta-casein enriched casein micelles (SAXS): https://www.sciencedirect.com/science/article/pii/S0268005X24000857
- Calcium caseinate micelle reassembly study: https://pubmed.ncbi.nlm.nih.gov/38395543/
- Experimental gelation SAXS archive: https://zenodo.org/records/22044394
- Martini 3 official: https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html
