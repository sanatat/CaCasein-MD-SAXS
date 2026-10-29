# Public repository audit — 29 October 2026

## Checks

- Python modules parse; local full-project suite: **36 passed**.
- V5 computations checked independently against an analytic sphere form factor and exact pair-distance invariance.
- Saved V5 results include CSV, figures and `results_summary.json`.
- Cached bytecode, JIT compiled artifacts, temporary test outputs, and virtual environments excluded from public distribution.
- Basic static scan found no obvious secrets or email addresses. This is not a formal security audit.

## Scientific claims that are NOT established

- No quantitative casein–calcium binding affinity or converged equilibrium ensemble.
- No realistic protein micelle structure, explicit-solvent validated force field, or Martini 3 run.
- No direct comparison against Shape2SAS executable.
- No experimental SAXS validation, q calibration or physically mapped time-scale.
- An identifiability proof for six-point mathematical examples must not be generalized to exact casein conformations.
- V5 conducts **no new MD**; its V4 component is a retrospective analysis of existing saved MD results.

## Suggested next work

1. Obtain and document suitable force-field parameters for phosphorylated casein and Ca²⁺.
2. Demonstrate sampling and box-size convergence with independent equilibrated starting conditions.
3. Run Shape2SAS independently on matched synthetic geometries for cross-package checks.
4. Analyze measurement-condition-matched SAXS/USAXS data with uncertainty propagation.
5. Decide on a code license before encouraging third-party reuse.
