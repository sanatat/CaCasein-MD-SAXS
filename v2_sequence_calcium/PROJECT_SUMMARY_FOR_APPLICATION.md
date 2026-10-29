# Application-ready project explanation (no exaggerated claims)

**Project:** Sequence-Informed β-Casein Assembly: CPU-Only Langevin MD and Synthetic Small-Angle Scattering

**Summary (for CV):** Implemented and executed a Python/Numba coarse-grained β-casein pilot using the published 209-residue bovine β-casein A2 primary sequence, localized phosphorylated-segment attraction proxies, 8 protein chains (168 beads), and nine independent Langevin trajectories. Quantified site-specific interchain contacts and calculated synthetic SAXS-like Debye profiles directly from simulated coordinates, including independent-replica variability and finite-difference force tests. The model is phenomenological, not yet calibrated against measured SAXS or explicit calcium-ion chemistry.

**What makes it interesting:** Changing the effective interchain bridge attraction strengthened local phosphosite contacts, whereas whole-aggregate synthetic SAXS shape changed relatively little across three short runs. This raises a realistic and testable inverse-problem question: which site-specific molecular rearrangements can small-angle scattering actually resolve?

**What it does not show:** No experimental validation, no fully parameterized Martini molecular dynamics, no actual concentrations of calcium, no predictive casein micelle formation, no equilibrium demonstration, no biological efficacy claims.

**Ideas for next stage:** Validate against SAXS under genuinely matched solution conditions; obtain and document force-field/chemical parameterization; quantify identifiability with multi-contrast SAXS/SANS; compare isolated beta-casein with beta/kappa-casein mixtures; run longer convergence and alternate initial-state experiments; use Bayesian model comparison with realistic noise.

**Sources:**
- Mature beta-casein A2 phosphorylated primary sequence, BIOPEP-UWM: https://biochemia.uwm.edu.pl/biopep/proteinsdatacard.php?zm_ID=1843
- Calcium-induced beta-casein aggregation report (2025): https://doi.org/10.1016/j.jcis.2025.137916
- Larsen et al. Shape2SAS (2023): https://doi.org/10.1107/S1600576723005848
- Public 2026 SAXS dataset with different experimental conditions: https://doi.org/10.5281/zenodo.22044394
