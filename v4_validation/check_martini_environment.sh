#!/usr/bin/env bash
# Read-only checks; does not install software or claim successful Martini modelling.
set -u
printf '=== GROMACS / Martini 3 readiness check ===\n'
if command -v gmx >/dev/null 2>&1; then
  echo '[FOUND] gmx:' "$(command -v gmx)"
  gmx --version | head -18
else
  echo '[MISSING] GROMACS command gmx'
  echo 'On Linux Mint/Ubuntu: sudo apt update && sudo apt install gromacs'
  echo 'Or follow official instructions for a newer supported version.'
fi
if command -v martinize2 >/dev/null 2>&1; then
  echo '[FOUND] martinize2:' "$(command -v martinize2)"
  martinize2 -h 2>&1 | head -6
else
  echo '[MISSING] martinize2 (check official Martini 3 protein tutorial)'
fi
if command -v python3 >/dev/null 2>&1; then
  python3 - <<'PY'
try:
    import numba, numpy, scipy
    print('[READY] Toy simulation Python deps: NumPy=%s Numba=%s SciPy=%s' % (numpy.__version__,numba.__version__,scipy.__version__))
except ImportError as e:
    print('[MISSING] Python package:',e)
PY
fi
printf '\nCritical: phosphorylated beta-casein topology and Ca2+ parameter validation are separate scientific prerequisites.\n'
