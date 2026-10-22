#!/bin/sh
# Build the separate GPL-licensed PROVEAN component outside its source tree.
set -eu
mutabind2_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mutabind2_prefix=${1:-"$mutabind2_root/external/provean-1.1.5-modern"}
mutabind2_psiblast=${MUTABIND2_PSIBLAST:-$(command -v psiblast)}
mutabind2_blastdbcmd=${MUTABIND2_BLASTDBCMD:-$(command -v blastdbcmd)}
mutabind2_cdhit=${MUTABIND2_CDHIT:-$(command -v cd-hit)}
mutabind2_build=$(mktemp -d /tmp/mutabind2_provean_build_XXXXXX)
trap 'rm -rf "$mutabind2_build"' EXIT HUP INT TERM
cd "$mutabind2_build"
"$mutabind2_root/vendor/provean-1.1.5/configure" \
  --prefix="$mutabind2_prefix" \
  PSIBLAST="$mutabind2_psiblast" \
  BLASTDBCMD="$mutabind2_blastdbcmd" \
  CDHIT="$mutabind2_cdhit"
make -j2
make install
