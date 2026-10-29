# Path to a physically grounded Martini 3 casein–calcium simulation

**Not executed yet.** This is a concrete plan and prerequisites, not a generated Martini 3 system.

### 1. Prerequisites
- Install **GROMACS** from its official installation instructions or distribution packages; check with `gmx --version`.
- Install **Martinize2** as documented by the official Martini initiative; verify its available Martini 3 force field definitions (`martinize2 -h`).
- Validate the **mature β-casein A2 sequence**, **all phosphorylation sites**, charge/protonation states, and an ensemble of initial conformations. This is an intrinsically disordered protein: one folded PDB and an elastic network should not be assumed representative.
- Confirm a **validated, version-specific Martini 3 phosphoserine** representation and verify Ca²⁺ interaction parameters using independent literature benchmarks. The 2025 Martinize2 paper says these PTM parameters were missing for Martini 3 at the time; verify any newer updates rather than blindly assuming availability.

### 2. First official tutorial before the casein model
Complete the official Martini 3 protein tutorial on its provided benchmark protein first:
https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html

It demonstrates a `martinize2` command with `-ff martini3001`, `-dssp`, and protein preparation. **That example does not produce a valid phosphorylated casein model.** Do not transfer its elastic restraints to disordered casein without physical justification.

### 3. Model and setup acceptance checklist
- Parameter provenance and versions in a machine-readable manifest.
- Correct topology counts/masses/charges; appropriate counterions; solvent and ionic strength matched to experimental conditions.
- Minimization, thermalization, then independent initializations and replicas.
- Timestep convergence, energy conservation sanity, density, box-size dependence, convergence diagnostics, and phosphate/Ca controls.
- Retain realistic ion concentrations in physical units and use appropriate reference for how free calcium is measured.
- Rg, site occupancy, ion residence, cluster populations, and time autocorrelation with block analyses.
- Predict SAXS with calibrated length, scattering contrast, hydration-shell and solvent effects as appropriate, and account for solution concentration and instrumental resolution.

### 4. Experimental comparison
Public record https://zenodo.org/records/22044394 contains casein scattering under **enzymatic gelation** at pH 6.6, 38°C, with experimental q between ~0.002 and 3 nm^-1. Those data may be inspected as an illustrative out-of-domain comparison, but they cannot validate calcium titration or phosphoserine affinities. Direct dataset download was **not possible in the current runtime**. No fit, R², chi², or claim of agreement is reported.

### 5. Why a dry-run-only `martinize2` command was not generated
Without a validated phosphorylated starting structure and force-field parameters, generating an apparently runnable topology encourages an unjustified biochemical claim. The official tutorial is a safe first working reference.

Sources:
- https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html
- https://elifesciences.org/articles/90627v1
- https://doi.org/10.1039/D0FO00844C
- https://doi.org/10.1016/j.foodchem.2010.12.121
