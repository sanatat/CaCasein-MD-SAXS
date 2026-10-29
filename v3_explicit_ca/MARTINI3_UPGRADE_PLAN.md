# Moving to physical Martini 3 — NOT completed

**There is NO real Martini 3 simulation in this repository.** The computer used for v3 lacks GROMACS and a Martini 3 protein topology, and calibrated phosphoserine / Ca–phosphate coordination parameters have not been established. Passing this boundary without appropriate tools, data, and independent validation would produce misleading scientific claims.

## What is required

1. Obtain a structural **ensemble** of mature 209-residue bovine beta-casein A2 (verify variant and phosphorylation/protonation). Casein is largely disordered, so one arbitrary PDB starting conformation would be unrepresentative.
2. Use the official GROMACS + martinize2 workflow to create Martini-3-compatible peptide topology. Example **for an already validated input protein PDB**:

   ```bash
   martinize2 -f protein_input.pdb -x protein_cg.pdb -o protein_only.top \
      -ff martini3001 -p backbone -dssp
   ```

   This does not by itself produce a valid **phosphorylated beta-casein topology**. Check residue/PTM handling, whether martinize2 recognizes every modified site, and all parameter compatibility.
3. Check up-to-date phosphoserine PTM parameter support. The 2025 eLife description of martinize2 explicitly notes Martini 3 PTM parameters were not implemented at the time of writing. Do not silently assign ordinary serine beads to phosphorylated residues.
4. Parameterize and benchmark calcium binding and phosphate chemistry from published calcium–CPP experiments; explicit Ca2+ ions alone do not ensure realistic coordination.
5. Add Martini 3 water, ion species, box composition, appropriate ionic strength and charge neutrality, balance timescale, and check trajectory convergence.
6. Compute SAXS after validating physical length scale, excess-electron-density contrast, solvent contribution and structure factor (if concentrated). Compare to **matched** experimental data, not the unrelated enzymatic gelation Zenodo archive as if it were a fit.

## Source links

- Martini 3 validated methods: https://www.nature.com/articles/s41592-021-01098-3
- Official protein tutorial: https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html
- Martini 3 parameters: https://cgmartini.nl/docs/downloads/example-input-files/md-parameters.html
- Martinize2 PTM limitations: https://elifesciences.org/articles/90627v1
- Calcium–casein phosphopeptide atomistic evidence: https://doi.org/10.1039/D0FO00844C
