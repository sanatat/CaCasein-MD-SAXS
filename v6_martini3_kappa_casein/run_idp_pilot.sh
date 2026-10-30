#!/usr/bin/env bash
# Martini3-IDP continuation; deliberately requires review between minimization and MD.
set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-help}"
THREADS="${THREADS:-2}"
for command_name in gmx martinize2 insane python3; do
  command -v "$command_name" >/dev/null || { echo "Missing executable: $command_name" >&2; exit 2; }
done
if [[ "$MODE" == 'prepare' ]]; then
  test -s inputs/tail_64.pdb || { echo 'First run the standard workflow to obtain and verify inputs/tail_64.pdb' >&2; exit 2; }
  mkdir -p work_idp_v3
  cp work/martini_v3.0.0*.itp work_idp_v3/
  (
    cd work_idp_v3
    # Case-sensitive force field ID, as reported by Vermouth 0.15.0 -list-ff
    martinize2 -f ../inputs/tail_64.pdb -x kappa_idp_cg.pdb -o topol.top \
      -ff martini3IDP -p backbone -ss C -id-regions 1:64 -idr-tune \
      > martinize.log 2>&1 || { tail -n 60 martinize.log; exit 1; }
    test -s molecule_0.itp
    insane -f kappa_idp_cg.pdb -o system.gro -p insane.top \
      -pbc cubic -d 4 -sol W -salt 0.15 > insane.log 2>&1 || { tail -n 60 insane.log; exit 1; }
    python3 ../scripts/wire_topology.py assemble topol.top insane.top system.top system.gro
    gmx grompp -f ../mdp/min.mdp -c system.gro -p system.top -o min.tpr -maxwarn 0
    gmx mdrun -deffnm min -nt "$THREADS"
    echo 'Minimization complete; examine work_idp_v3/min.log, especially for LINCS warnings.'
    echo 'Do not mistake minimized coordinates for a validated IDP ensemble.'
  )
elif [[ "$MODE" == 'continue-after-review' ]]; then
  (
    cd work_idp_v3
    for f in min.gro min.log system.top; do test -s "$f" || { echo "Missing $f" >&2; exit 2; }; done
    # Caller is responsible for reviewing any minimization LINCS warnings before continuing.
    gmx grompp -f ../mdp/eq_idp_safety.mdp -c min.gro -p system.top -o eq_safety.tpr -maxwarn 0
    gmx mdrun -deffnm eq_safety -nt "$THREADS"
    test -s eq_safety.cpt
    gmx convert-tpr -s eq_safety.tpr -extend 490 -o eq500.tpr
    gmx mdrun -s eq500.tpr -cpi eq_safety.cpt -deffnm eq_safety -nt "$THREADS" -append
    gmx grompp -f ../mdp/prod_idp_probe.mdp -c eq_safety.gro -t eq_safety.cpt \
      -p system.top -o idp_probe100.tpr -maxwarn 0
    gmx mdrun -deffnm idp_probe100 -nt "$THREADS"
    gmx check -f idp_probe100.xtc
    if grep -iE '^Step [0-9]+.*LINCS WARNING|^LINCS WARNING|not finite|Fatal error' \
      eq_safety.log idp_probe100.log; then
      echo 'Review numerical warnings before interpreting outputs' >&2; exit 3
    fi
  )
  echo 'IDP 500 ps equilibration and 100 ps production probe completed; no convergence claim.'
else
  echo 'Usage: bash run_idp_pilot.sh prepare|continue-after-review' >&2
  exit 2
fi
