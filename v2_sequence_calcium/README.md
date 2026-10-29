# CaCasein-MD-SAXS v2 — sequence-informed calcium-bridge *proxy* + synthetic SAXS

**Status: executed, sequence-informed coarse-grained Langevin simulation; preliminary, uncalibrated. NOT a validated β-casein/milk micelle MD model. No claim of actual calcium concentrations, Martini/GROMACS implementation, or experimental SAXS validation.**

## Scientific design

This v2 addresses a more specific question than the old 72-bead toy model: **Can a localized bridge attraction between phosphorylated β-casein segments alter the organization of an initially bundled assembly, and would forward-calculated SAXS reflect those changes?**

- **Primary sequence:** mature bovine β-casein A2 from BIOPEP-UWM database ID 1843 (UniProt P02666); 209 amino acids. Five annotated mature-sequence phosphoserines at residues **15, 17, 18, 19, 35**. https://biochemia.uwm.edu.pl/biopep/proteinsdatacard.php?zm_ID=1843
- **Resolution:** 21 blocks / chain (mostly 10 amino acids each); 8 chains = 168 coarse-grained beads total. Protein compositional information is preserved at a coarse level but **neither backbone conformations nor actual structures/side-chain packing are resolved**.
- **Dynamics:** real stochastic Langevin MD using BAOAB with bonded harmonic restraints, WCA excluded-volume, weak sequence-dependent hydrophobic attraction, heuristic screened charge repulsion/attraction, **interchain** Gaussian bridge interactions localized to phosphorylated blocks. Arbitrary reduced length, energy/time; 3 independent stochastic runs for each condition.
- **Scenarios:** `off`, `low`, `high` vary the *bridge interaction strength*: 0.0, 0.7, 1.4 **kBT**. These are *not* Ca²⁺ concentrations, and the model lacks explicit ions, calcium-phosphate mineral, pH titration, ion-binding kinetics, ionic competition, and force-field validation.
- **Initialization:** identical (apart from seed-level perturbations) preassembled 8-chain elongated bundles. The experiment is **short-time bundle retention and rearrangement**, not equilibrium micelle nucleation or self-assembly from fully dispersed chains.
- **Observables:** contact clusters, radius of gyration of largest cluster, **localized interchain phosphosite contacts**, and *synthetic* Debye scattering from simulated snapshots. Independent clusters scatter incoherently (dilute mixture approximation); all beads have simple contrast *proxy* roughly proportional to represented residue count. **Neither an actual SAXS solvent contrast model nor a structure factor is calculated.**

### Executed runs and results

The bundled model retained all eight chains in its largest connected cluster across the three short runs for all conditions; thus **cluster size did not discriminate conditions**. The local phosphosite contact counts were more responsive. Each table row is the mean across **3 independent trajectory-level averages** from the final ~45% of saved frames. ± is **between-replica SD**, not a standard error. A contact is any pair of phosphate-containing coarse beads on separate chains within **1.75 reduced length units**.

| Bridge strength (proxy) | Localized phosphosite contact pairs, mean ± between-replica SD |
|---|---:|
| 0.0 kBT | 12.07 ± 1.88 |
| 0.7 kBT | 17.69 ± 0.68 |
| 1.4 kBT | 22.51 ± 2.07 |

These outcomes reflect **the behavior of the chosen phenomenological potential**. They are not predictive measurements of calcium-driven β-casein aggregation. The runs are short (3200 integration steps at dt 0.003 = 9.6 reduced time units); convergence, calibration, initial-state sensitivity and experimental benchmarking remain open.

A descriptive contact-cutoff sensitivity check in `results/contact_cutoff_sensitivity.csv` showed off < low < high at cutoffs 1.5, 1.75 and 2.0 (reduced units). This is still **not** an experimental test or an equilibrium result. The simulated SAXS curves for off/low/high overlap substantially; at q < 1 (reduced inverse length), the maximum relative difference between high/off mean normalized shapes is about **2%** in these trajectories. This does not establish measurement detectability because simulated experimental noise and instrument effects were not included.

The SAXS curves are forward calculations on the actual simulated positions, but the q axis is **inverse reduced length (NOT nm⁻¹)**. Without calibrated bead length and electron contrast, compare model *relative shapes only*, not absolute structural scales.

## Run (Linux Mint, CPU only)

```bash
cd CaCasein-MD-SAXS/v2_sequence_calcium
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r ../requirements.txt
python -m pytest -q tests
python casein_v2.py --steps 3200 --replicas 3 --stride 100 --out reproduced
python analyze_local_contacts.py --results reproduced
```

On an older CPU use `OPENBLAS_NUM_THREADS=1` before Python, and lower `--steps 300 --replicas 1` to smoke-test first. The included `results/` directory contains the independently **executed reference 3×3 trajectories**, full CSVs, figures, metadata, per-replica contact summaries, and the sequence feature table. Dependencies are NumPy, Numba, Matplotlib and pytest; no OpenMM or GROMACS is needed. This is an explicitly custom bead-spring MD integrator, **not** a force field that can be transferred to atomistic beta-casein.

## Experimental casein SAXS (available online, **NOT matched or validated here**)

A relevant open dataset was published on Aug 21, 2026 (DOI **10.5281/zenodo.22044394**), containing **experimental** background-subtracted USAXS/SAXS during **chymosin-mediated gelation of reconstituted skim milk at pH 6.6 and 38 °C**. The dataset archive is `substracted_merged_USAXS_SAXS.7z` (2.1 MB; MD5 `652981bb9f72efbf7949127b77393826`).

https://zenodo.org/records/22044394

**Important mismatch:** these are whole-milk *enzymatic gelation* data; our model is only 8 beta-casein A2 chains with a parameterized calcium bridge proxy. Therefore their absolute SAXS curves cannot validate this model, even if overlaid. No raw experimental data were fetched in this environment due to external binary-download restrictions.

From your Linux laptop (with internet) you can fetch the source data independently:

```bash
python experimental_data_helper.py fetch --out external/zenodo_casein_saxs_2026.7z
sudo apt install p7zip-full
7z x external/zenodo_casein_saxs_2026.7z -oexternal/extracted
```

The archive's internal file layout is not verified here. After finding an individual two-column `(q, I)` experimental text file, an *exploratory* (not validated) comparison can be plotted with a **justified, explicitly supplied** length-scale conversion:

```bash
python experimental_data_helper.py overlay \
 --experiment path/to/single_2_column_experimental_saxs.txt \
 --q-unit nm^-1 --length-nm-per-reduced YOUR_JUSTIFIED_LENGTH_NM \
 --condition off --model results/saxs_curve_replicas.csv
```

Do **not** tune the nm-per-bead conversion to improve visual agreement and call the result experimental validation. For actual experimental validation, obtain SAXS data from a matched isolated β-casein titration, including concentration, ionic strength, temperature, sample conditions, q-resolution, absolute scaling, and uncertainty.

## Explicit reproducibility/limitations

1. The phosphate modifications are sequence-informed; no experimentally established calcium equilibrium constants, bonded structure, local ion coordination, or electrostatic calibration is implemented.
2. Three short runs do not establish thermal equilibration, independent effective samples, or predictive uncertainty. No statistical-significance claims.
3. The same initially bundled geometry across conditions can bias results toward retention. Need dissolved starting conditions and finite-size studies.
4. Changing the user-defined bridge attraction unsurprisingly increases phosphosite contacts in this formulation; it is **an internal sensitivity test**, not discovery of a new biochemical mechanism.
5. The scattering uses an isotropic point-bead Debye formula. No experimentally calibrated form factor, electron-density solvent contrast, hydration shell, resolution convolution, or polydisperse real milk micelles.
6. Because there are only 8 chains and 21 beads per chain, no inference can be made about a full 100–200 nm casein micelle.

## Citable references

- BIOPEP-UWM modified mature bovine β-casein A2 primary sequence: https://biochemia.uwm.edu.pl/biopep/proteinsdatacard.php?zm_ID=1843
- Kaur & Bhattacharya (2025), on calcium-associated phosphorylated β-casein aggregation: https://doi.org/10.1016/j.jcis.2025.137916
- Larsen et al. (2023), Shape2SAS Debye forward-scattering method: https://doi.org/10.1107/S1600576723005848
- Bolinsson et al. (2025), AF4-UV-SAXS on casein micelles, data by request: https://www.sciencedirect.com/science/article/pii/S0268005X25003376
- Bauland et al. (2026), public real USAXS/SAXS data (enzymatic milk gelation, not matched to model): https://doi.org/10.5281/zenodo.22044394

**Accurate CV summary:** "Developed and executed a CPU-only, sequence-informed coarse-grained β-casein pilot (8 chains × 21 beads, 3 interaction conditions × 3 independent replicas), with localized phosphorylation-site bridging proxies, replica-level contact analysis, and synthetic Debye SAXS profiles. The phenomenological model is not experimentally calibrated."
