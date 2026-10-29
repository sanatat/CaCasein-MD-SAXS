# CaCasein-MD-SAXS — V5: SAXS forward validation, structural non-identifiability, and sampling audit

**Status: EXECUTED computational methods pilot (29 October 2026).**

**Warning:** This is a *sequence-informed casein-related methods portfolio project*, not a validated quantitative model of calcium–casein affinity, a micelle reconstruction, or an experimental SAXS analysis. V5 reuses saved V4 MD data and independently executes new SAXS / inverse-problem numerical experiments. V5 does **not** run additional MD trajectories, Martini 3, the actual Shape2SAS program, or an experimental dataset fit.

## Scientific question

Can small-angle X-ray scattering, considered alone, uniquely identify a nanoscale arrangement, and how much could imperfect sampling and the use of external structural information affect an inference?

This question relates methodologically to Kirkensgaard and Larsen's Shape2SAS paper (Larsen et al., *J. Appl. Crystallogr.* 56, 1287–1294, 2023), but the exact homometric examples here are *not* casein structures. Our Debye implementation is independent, not copied from their software.

## Experiments that were actually performed

### 1. Analytic sphere benchmark

- Implemented point-scatterer Debye intensity and compared it to the independent analytical form factor of a uniform sphere.
- Radius = 12 uncalibrated length units; three independent uniformly sampled spheres of 850 points each.
- Evaluated 90 q points from 0.004 to 0.20 inverse reduced length.
- **Root-mean-square absolute P(q) error: 0.002945.** Finite random discretization is a main contribution.
- The analytical sphere predicts Rg = 9.29516 reduced units; direct sampled-point Rg values are logged in `results_summary.json`. Low-q Guinier estimates are also logged.
- This is a validation of one point-scatterer numerical implementation, *not* validation of a protein force field or Shape2SAS itself.

### 2. Mathematical counterexample: exactly homometric geometries

- A = [0, 1, 4, 10, 12, 17]; B = [0, 1, 8, 11, 13, 17] along the x-axis, identical positive point contrasts.
- They are not related by rigid translation/reflection, yet their fifteen unordered pair separations are identical.
- Orientationally averaged Debye SAXS depends only on these distances (and the point contrasts), so intensity is the same for *all q*, regardless of measurement precision.
- Across 180 sampled q values, max absolute numerical discrepancy: **3.33 × 10⁻¹⁶**.
- Same Rg = **6.155395** reduced units, and exactly identical pair-distance histograms.
- This is *proof of possible structural nonuniqueness*; it is **not** evidence that actual beta-casein conformers are exactly homometric.

### 3. An explicit inverse problem with measurement noise

- Synthetic truth: A. Models: A, B and distinguishable C. Known normalized scale and Gaussian independent noise σ = 0.025 (absolute P(q) units).
- 60 independent synthetic observation replicates, 180 q values per replica; likelihood from a Gaussian-noise model, initially equal candidate priors.
- Average posterior probabilities: A ≈ **0.50**, B ≈ **0.50**, C ≈ **0.00**. A and B have exactly equal likelihood for every possible observation.
- An additional *hypothetical external prior* of (0.8,0.1,0.1) on (A,B,C) shifts A's posterior mean to ≈ **0.8889**. The choice is illustrative, not derived from MD or actual data. **Priors do not create SAXS information.**
- Normalization, point-scatterer assumptions and independent-noise model are intentional simplifications. Real scattering errors may be correlated or heteroscedastic.

### 4. Retrospective V4 molecular-dynamics audit (no new MD)

- Analyzed existing 15 saved V4 trajectories (5 conditions × 3 independent replicas) and associated scattering CSVs.
- For eight ion-like beads with the same toy attractive potential, late-phase mean Ca-like bead contacts differ strongly between randomly distributed ions (**0.00**) and near-site initialization (**3.30**) across 3 replicas, from V4. Those are simulation contacts, **not physical binding affinities**.
- Paired-seed RMS difference of **log-normalized SAXS** between bulk and near-site start averaged ≈ **0.0547** across the three replicas.
- Reports early/late per-replica contact changes and conservative Student-t intervals (n=3, unstable estimates). No equilibrium or convergence claim.
- V4 curves are protein-only synthetic scattering and their q axis is an *inverse reduced length unit*. Ion contrast, hydration, interparticle interference and experimental beam corrections are absent.

## Figures

1. `results/01_sphere_benchmark.png`: point-Debye versus analytical sphere, residuals.
2. `results/02_exact_homometry.png`: noncongruent geometry, exact pair-separation match, same SAXS.
3. `results/03_inverse_uncertainty.png`: example noisy data, posterior across 60 repeats.
4. `results/04_v4_reanalysis.png`: V4 initialization effect and scattering curves.

## Reproduction (Linux Mint or other Python 3.10+)

```bash
cd CaCasein-MD-SAXS/v5_saxs_identifiability
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements_v5.txt
python -m pytest -q tests
python run_v5.py --v4 ../v4_validation --out results_recomputed
```

This produces all V5 plots, CSV files, numerical diagnostics, and `results_summary.json`. To run V1–V4 tests too:

```bash
cd ..
python -m pytest -q tests v2_sequence_calcium/tests v3_explicit_ca/tests v4_validation/tests v5_saxs_identifiability/tests
```

## Optional comparison with Shape2SAS (not part of the executed results)

Shape2SAS by Larsen, Brookes, Pedersen and Kirkensgaard is publicly available at:
https://github.com/andreashlarsen/Shape2SAS

Its documentation gives a sphere run, from the downloaded repo directory containing `shape2sas.py`:

```bash
python shape2sas.py --subunit sphere --dimension 50
```

The sphere dimension is in **Å** and q is **Å⁻¹** in Shape2SAS. Provide the generated `Iq_<model_name>.dat` to our local independent comparator:

```bash
python check_shape2sas_reference.py \
  --shape2sas-iq /path/to/Iq_Model_0.dat \
  --sphere-radius-angstrom 50
```

This compares both curves against the independent analytical sphere P(q) using one global amplitude scale. Ensure the exact geometry, contrast, polydispersity, q range and default configuration are understood before interpreting the result. **The official Shape2SAS software was not downloaded or executed in the current environment**, so no direct cross-package match is claimed.

## Optional actual experimental dataset analysis (not yet executed)

A public data source, https://zenodo.org/records/22044394, contains time-resolved enzymatic gelation SAXS/USAXS. It is NOT a calcium titration with conditions matched to the toy MD; it should not be called experimental validation of V3/V4. An external real experimental file, after responsible sample selection, background processing and export, can be independently inspected as a CSV with columns `q,intensity,sigma`:

```bash
python analyze_experimental_saxs.py \
  --csv /path/to/prepared_experimental.csv \
  --q-unit 'nm^-1' \
  --guinier-qmax 0.08 \
  --source 'DOI and exact sample/timepoint/measurement conditions'
```

This only reports a simple low-q Guinier fit and flags `qmax*Rg <= 1.3`. It does **not** perform an MD model fit. The dataset archive was not accessible in our restricted Python run; no empirical data points have been made up.

## Project integrity and scientific limits

- V1–V4 use reduced-unit phenomenological potentials; V5 has **zero new MD runs**.
- The geometries A and B are idealized one-dimensional, six-point examples; not conformations of casein.
- The calculated intensity assumes isolated, orientationally averaged particles with matched point weights, no explicit solvent, negligible interparticle structure factors and no instrument smearing.
- The absolute uncertainty model is synthetic; its Bayesian model selection probabilities are conditioned on three preselected candidate shapes and priors.
- Length-to-nm, time-to-ns and ion-count-to-molarity calibration does not exist for this model.
- The results are a demonstration of **model criticism, validation, uncertainty and identifiability**, NOT a test of calcium binding or a peer-reviewed biological result.

## References

1. Larsen AH, Brookes E, Pedersen MC, Kirkensgaard JJK (2023). Shape2SAS. *Journal of Applied Crystallography*, 56:1287–1294. https://doi.org/10.1107/S1600576723005848
2. Official Shape2SAS code: https://github.com/andreashlarsen/Shape2SAS
3. Bauland et al., background-subtracted USAXS/SAXS casein-micelle enzymatic gelation dataset (2026). https://doi.org/10.5281/zenodo.22044394

## Authoring

Prepared as a reproducible computational research portfolio example. Scripts, tests and synthetic numerical outputs are fully included.
