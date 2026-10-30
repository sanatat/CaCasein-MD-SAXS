# V6 — Bovine κ-casein C-terminal Martini 3 / Martini3-IDP pilot

**Status (30 October 2026): executed on the author's Linux Mint laptop; preliminary technical simulation, not a converged or experimentally validated protein model.**

This is a new, protein-specific branch of the project; the earlier V1–V5 work remains independent and unchanged. V6 models the **unmodified bovine κ-casein mature residues 106–169 (64 residues)**, corresponding to the C-terminal segment of UniProt P02668. It is **not a whole casein micelle**.

## Verified local execution reported by the researcher

| System | Procedure completed | Outcome |
|---|---|---|
| Standard Martini 3 | Structure preparation; Martinize2, Insane, minimization, equilibration; 10 ns GROMACS 2023.3 production, dt = 10 fs | 1,000,000 integration steps; complete trajectory of 2,001 frames (0–10 ns) |
| Standard Martini 3 | Backbone-bead Rg from the 10 ns trajectory | Full-trajectory mean **2.0060 nm**; 8–10 ns mean **1.7311 nm**, SD **0.1528 nm** (framewise); final frame **1.8917 nm** |
| Martini3-IDP | Independent topology using Martinize2/Vermouth 0.15.0, force field identifier `martini3IDP`, and IDR residues 1–64; energy minimization | Converged at step 2236, Fmax ~97.97 kJ mol⁻¹ nm⁻¹; transient LINCS warnings at minimization steps 33–35 |
| Martini3-IDP | NVT equilibration to 500 ps at dt = 2 fs | 501 frames, no matching LINCS warnings during MD; mean temperature **299.492 K** during final 250 ps |
| Martini3-IDP | Separate 100 ps production test at dt = 5 fs | 101 frames; mean temperature **299.911 K**, potential energy average ~−760,076 kJ/mol; no matching LINCS warnings during MD |

Numbers are transcribed from the author's GROMACS terminal outputs and are **not independent reanalysis of uploaded trajectories**. Raw `.xtc`, `.tpr`, `.edr` and `.cpt` files are held on the researcher's local computer, **not published here**. Numerical stability is not thermodynamic equilibration.

### Standard Martini 3: backbone Rg by 2 ns interval

| Production window (ns) | Mean Rg (nm) | Framewise SD (nm) |
|---|---:|---:|
| 0–2 | 2.3386 | 0.0904 |
| 2–4 | 2.2082 | 0.1045 |
| 4–6 | 2.0204 | 0.0927 |
| 6–8 | 1.7320 | 0.1635 |
| 8–10 | 1.7311 | 0.1528 |

The observed decrease in Rg is **not** evidence that this represents the true dimensions of the κ-casein brush: simulation length is short for an intrinsically flexible region, the input used a single low-confidence predicted conformation, and no independent ensemble validation or experimental matching has been completed. The Martini3-IDP pilot is much shorter than the standard Martini 3 trajectory, so **no force-field superiority or Rg comparison is claimed**.

## Reproducibility / scope

The source pilot workflow is based on the local `CaCasein-MD-SAXS-v6-Martini3-pilot` folder: `run_martini3_pilot.sh`, `scripts/`, `mdp/`, `tests/` and `results/`. The repaired topology assembler maintains separate Martinize2 and Insane prototype topologies. **Do not use `epsilon-rf = inf` in the MDP files; the corrected local runs use `epsilon-rf = 0`.** The original pilot automation reproduces only a short *standard Martini 3* smoke test; the Martini3-IDP comparative run was prepared and executed in a separate local `work_idp_v3/` directory.

**Original run configuration and derived results must not be implied reproducible from this README alone.** Publish scripts/configurations and a bounded log excerpt (without raw trajectories) before presenting this as a fully archived workflow.

### Scientific limitations

- Unmodified isolated 64-residue tail; **no phosphorylation, glycosylation, native casein micelle, or quantitative calcium association**.
- Coarse-grained geometry-based synthetic scattering **is not experimental or quantitatively calibrated SAXS**.
- The 10 ns standard Martini 3 trajectory is a single initialization and not a converged ensemble. Rg trends may reflect force-field bias, initial coordinates, sampling, or periodic-boundary handling.
- LINCS warnings occurred during **IDP minimization**, though no matching MD LINCS warnings were found in the follow-up equilibration / 100 ps production tests.
- Do not compare absolute potential energies across different box compositions, bonded parameters, or force fields.

## Source context

- [Martini 3 protein tutorials](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/)
- [Martini IDR tutorial](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut3.html)
- [UniProt P02668](https://www.uniprot.org/uniprotkb/P02668/entry)
- [Martini3-IDP research](https://doi.org/10.1038/s41467-025-58199-2)
