#!/usr/bin/env bash
# Compatibility wrapper. New automation should use reproduce_verification.sh.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mapped=()
for arg in "$@"; do
    case "$arg" in
        mono_hex) mapped+=(monodomain_cartesian) ;;
        mono_tet) mapped+=(monodomain_tet_generic) ;;
        eikonal_hex) mapped+=(eikonal_cartesian) ;;
        eikonal_tet) mapped+=(eikonal_tet_generic) ;;
        bidomain_hex) mapped+=(bidomain_cartesian) ;;
        bidomain_tet) mapped+=(bidomain_tet_generic) ;;
        bath_hex) mapped+=(bath_bidomain_cartesian) ;;
        bath_tet) mapped+=(bath_bidomain_tet_conformal) ;;
        coupling1D3D_hex) mapped+=(purkinje_monodomain_coupled) ;;
        niederer_hex) mapped+=(niederer_cartesian) ;;
        *) mapped+=("$arg") ;;
    esac
done
exec bash "$REPO_ROOT/reproduce_verification.sh" "${mapped[@]}"
