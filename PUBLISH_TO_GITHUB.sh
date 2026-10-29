#!/usr/bin/env bash
# Publish the audited CaCasein-MD-SAXS portfolio to its existing GitHub repository.
# This script never force-pushes and never creates a new GitHub repository.
set -euo pipefail

REPO='https://github.com/sanatat/CaCasein-MD-SAXS.git'
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

if ! command -v git >/dev/null 2>&1; then
  echo 'ERROR: git is missing. Install using: sudo apt install git' >&2
  exit 1
fi
for f in README.md requirements.txt v5_saxs_identifiability/run_v5.py .gitignore; do
  test -f "$f" || { echo "ERROR: missing $f" >&2; exit 1; }
done
if ! git config user.name >/dev/null || ! git config user.email >/dev/null; then
  cat >&2 <<'EOT'
ERROR: Git author identity is not configured.
Configure your own Git identity (NOT your GitHub password):
  git config --global user.name "Your Name"
  git config --global user.email "you@example.com"
Then rerun this script.
EOT
  exit 1
fi

if [ ! -d .git ]; then
  git init -b main
fi
current_branch="$(git symbolic-ref --quiet --short HEAD || true)"
if [ "$current_branch" != 'main' ]; then
  echo "ERROR: current branch '$current_branch' is not main; stopping safely." >&2
  exit 1
fi
if git remote get-url origin >/dev/null 2>&1; then
  existing="$(git remote get-url origin)"
  case "$existing" in
    "$REPO"|git@github.com:sanatat/CaCasein-MD-SAXS.git) ;;
    *) echo "ERROR: origin points to unexpected remote: $existing" >&2; exit 1 ;;
  esac
else
  git remote add origin "$REPO"
fi

# Commit only the reviewed project folder. Honor .gitignore.
git add -A
if ! git diff --cached --quiet; then
  git commit -m 'Initial public computational methods portfolio release (v5)'
else
  echo 'No local changes to commit.'
fi

echo 'Pushing to GitHub. GitHub may ask you to log in on your own computer.'
git push -u origin main
printf '\nSUCCESS: https://github.com/sanatat/CaCasein-MD-SAXS\n'
