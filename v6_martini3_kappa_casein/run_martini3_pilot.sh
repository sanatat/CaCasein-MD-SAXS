#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
MODE="${1:-preflight}"
mkdir -p results inputs work
python3 scripts/reference.py --out results >/dev/null
if [[ "$MODE" == 'preflight' ]]; then
  echo 'Checking tools: GROMACS, martinize2, insane, curl, unzip'
  for x in gmx martinize2 insane curl unzip; do
    if command -v "$x" >/dev/null 2>&1; then printf 'OK   %-12s %s\n' "$x" "$(command -v "$x")"; else printf 'MISS %-12s\n' "$x"; fi
  done
  echo 'The repaired script supports run (fresh) and resume (reuse existing coordinates/topology).'
  exit 0
fi
if [[ "$MODE" != 'run' && "$MODE" != 'resume' ]]; then echo 'Usage: bash run_martini3_pilot.sh [preflight|run|resume]' >&2; exit 2; fi
command -v gmx >/dev/null || { echo 'Missing gmx' >&2; exit 3; }
run_md() {
  (
    cd work
    echo 'Running GROMACS minimization ...'
    gmx grompp -f ../mdp/min.mdp -c system.gro -p system.top -o min.tpr -maxwarn 0
    gmx mdrun -deffnm min -nt 2
    echo 'Running GROMACS equilibration ...'
    gmx grompp -f ../mdp/eq.mdp -c min.gro -p system.top -o eq.tpr -maxwarn 0
    gmx mdrun -deffnm eq -nt 2
    echo 'Running GROMACS short technical smoke test ...'
    gmx grompp -f ../mdp/prod.mdp -c eq.gro -t eq.cpt -p system.top -o prod.tpr -maxwarn 0
    gmx mdrun -deffnm prod -nt 2
  )
  python3 scripts/analyze_martini.py --tpr work/prod.tpr --xtc work/prod.xtc --out results
  echo 'Technical Martini 3 pilot completed. NOT an experimental SAXS or calcium-binding validation.'
}

if [[ "$MODE" == 'resume' ]]; then
  for f in work/topol.top work/system.top work/system.gro work/molecule_0.itp; do
    [[ -s "$f" ]] || { echo "Missing required existing file $f; use run instead" >&2; exit 4; }
  done
  # Legacy run wrote Insane's prototype to work/system.top. Read that before overwrite.
  INPUT_TOP=work/insane.top
  [[ -s "$INPUT_TOP" ]] || INPUT_TOP=work/system.top
  python3 scripts/wire_topology.py assemble work/topol.top "$INPUT_TOP" work/system.top work/system.gro
  run_md
  exit 0
fi
for x in martinize2 insane curl unzip; do command -v "$x" >/dev/null || { echo "Missing $x" >&2; exit 3; }; done
if [[ ! -s inputs/tail_64.pdb ]]; then
  echo 'Downloading AlphaFold P02668 prediction and extracting tail.'
  python3 scripts/prepare_structure.py --out inputs
fi
if [[ ! -s work/martini_v3.0.0.itp || ! -s work/martini_v3.0.0_solvents_v1.itp || ! -s work/martini_v3.0.0_ions_v1.itp ]]; then
  mkdir -p work/forcefield_download
  curl -fLSs --retry 3 --max-time 120 \
    https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip \
    -o work/forcefield_download/martini_v300.zip
  (cd work/forcefield_download && unzip -qo martini_v300.zip)
  for ff in martini_v3.0.0.itp martini_v3.0.0_solvents_v1.itp martini_v3.0.0_ions_v1.itp; do
    candidate="$(find work/forcefield_download -type f -name "$ff" | head -1)"
    [[ -n "$candidate" ]] || { echo "Missing $ff in official archive" >&2;exit 4; }
    cp "$candidate" work/"$ff"
  done
fi
# Unmodified kappa-casein C-terminal peptide. No phosphosites, glycans or explicit Ca2+.
SS="$(python3 -c "import sys;sys.path.insert(0,'scripts');from reference import TAIL;print('C'*len(TAIL))")"
(
  cd work
  martinize2 -f ../inputs/tail_64.pdb -x kappa_c_tail_cg.pdb -o topol.top \
    -ff martini3001 -p backbone -ss "$SS" \
    > martinize2.stdout 2> martinize2.stderr
  # Generate coordinates and molecule counts to a separate prototype topology.
  # NEVER let Insane overwrite the final combined system.top.
  insane -f kappa_c_tail_cg.pdb -o system.gro -p insane.top \
    -pbc cubic -d 4 -sol W -salt 0.15 > insane.stdout 2> insane.stderr
)
python3 scripts/wire_topology.py assemble work/topol.top work/insane.top work/system.top work/system.gro
run_md
