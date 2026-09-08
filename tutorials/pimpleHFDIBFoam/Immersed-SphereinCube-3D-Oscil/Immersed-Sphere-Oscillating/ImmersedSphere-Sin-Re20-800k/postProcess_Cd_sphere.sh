#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Extract Fx from bodiesInfo/<time>/body<ID>.info (FCoupling_F)
# and compute Cd for 3D oscillating sphere:
#
#   Cd = Fx / (0.5*rho*Uref^2*Aref)
#   Aref = pi*D^2/4
#
# Plots:
#   Present immersed-boundary Cd = solid line
#   Gilmanov & Sotiropoulos reference = points from ./Cd.dat
# ============================================================

# ----------------------------
# User settings: 3D sphere case
# ----------------------------
BODY_ID=0
RHO=1.0

# Sphere diameter used in the benchmark/reference scaling
D=1.0

# Reference velocity.
# For sinusoidal motion, usually Uref = 2*pi*f*A.
# For the Gilmanov & Sotiropoulos style benchmark, this is commonly pi/4.
UREF=0.7853981633974483

# Reference data file: red-circle Gilmanov & Sotiropoulos data
REF_CD_FILE="./Cd.dat"

# Optional y-range
YRANGE_MIN=${YRANGE_MIN:--30}
YRANGE_MAX=${YRANGE_MAX:-30}

# ----------------------------
# Paths
# ----------------------------
BODIES_DIR="bodiesInfo"
OUTDIR="postProcessing/forces"

FORCES_DAT="$OUTDIR/forces_body${BODY_ID}.dat"
DRAG_DAT="$OUTDIR/dragFx_body${BODY_ID}.dat"
CD_DAT="$OUTDIR/Cd_body${BODY_ID}.dat"
REF_SORTED_DAT="$OUTDIR/reference_Gilmanov_Sotiropoulos_Cd_sorted.dat"

PLOT_CD_GP="$OUTDIR/plot_Cd_body${BODY_ID}.gp"
CD_PDF="$OUTDIR/Cd_body${BODY_ID}.pdf"

# ----------------------------
# Checks
# ----------------------------
if [[ ! -d "$BODIES_DIR" ]]; then
  echo "ERROR: Cannot find '$BODIES_DIR' in case root."
  echo "Run this script from the case root where bodiesInfo/ exists."
  exit 1
fi

if [[ ! -f "$REF_CD_FILE" ]]; then
  echo "ERROR: Cannot find reference file '$REF_CD_FILE'"
  echo "Expected Gilmanov & Sotiropoulos reference file: ./Cd.dat"
  exit 1
fi

mkdir -p "$OUTDIR"

# ----------------------------
# Derived values
# ----------------------------
AREF=$(awk -v D="$D" 'BEGIN{pi=atan2(0,-1); printf "%.15g", pi*D*D/4.0}')
DEN=$(awk -v rho="$RHO" -v u="$UREF" -v a="$AREF" 'BEGIN{printf "%.15g", 0.5*rho*u*u*a}')

echo "Using:"
echo "  BODY_ID = $BODY_ID"
echo "  rho     = $RHO"
echo "  D       = $D"
echo "  Aref    = $AREF"
echo "  Uref    = $UREF"
echo "  den     = $DEN"
echo "  refCd   = $REF_CD_FILE"
echo

# ----------------------------
# 1) Extract force: time Fx Fy Fz
# ----------------------------
{
  echo "# time Fx Fy Fz"
  for d in $(find "$BODIES_DIR" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | sort -g); do
    f="$BODIES_DIR/$d/body${BODY_ID}.info"
    [[ -f "$f" ]] || continue

    vec=$(grep -m1 -E "^[[:space:]]*FCoupling_F" "$f" | sed -E 's/.*\(\s*([^)]*)\s*\).*/\1/')

    if [[ -n "${vec:-}" ]]; then
      echo "$d $vec"
    fi
  done
} > "$FORCES_DAT"

echo "Wrote: $FORCES_DAT"

# ----------------------------
# 2) Drag force Fx
# ----------------------------
awk '
  BEGIN{print "# time Fx"}
  $1 ~ /^#/ {next}
  NF>=2 {print $1, $2}
' "$FORCES_DAT" > "$DRAG_DAT"

echo "Wrote: $DRAG_DAT"

# ----------------------------
# 3) Cd = Fx / denominator
# ----------------------------
awk -v den="$DEN" '
  BEGIN{print "# time Cd"}
  $1 ~ /^#/ {next}
  NF>=2 {print $1, $2/den}
' "$DRAG_DAT" > "$CD_DAT"

echo "Wrote: $CD_DAT"

# ----------------------------
# 4) Sort Gilmanov & Sotiropoulos reference data
# ----------------------------
awk '
  $1 ~ /^#/ {next}
  NF>=2 {print $1, $2}
' "$REF_CD_FILE" | sort -g -k1,1 > "$REF_SORTED_DAT"

echo "Wrote: $REF_SORTED_DAT"
echo

# ----------------------------
# 5) Basic stats
# ----------------------------
awk '
  $1 ~ /^#/ {next}
  NF>=2 {
    Fx=$2
    if(n==0){min=Fx; max=Fx}
    if(Fx<min) min=Fx
    if(Fx>max) max=Fx
    sum+=Fx
    n++
  }
  END{
    if(n>0){
      mean=sum/n
      amp=0.5*(max-min)
      print "---- Fx stats ----"
      print "N       =", n
      print "Fx_mean =", mean
      print "Fx_min  =", min
      print "Fx_max  =", max
      print "Fx_amp  =", amp
    } else {
      print "No valid force rows found."
    }
  }
' "$DRAG_DAT"

awk '
  $1 ~ /^#/ {next}
  NF>=2 {t=$1; cd=$2}
  END{
    if(t!=""){
      print "---- Last Cd ----"
      print "Last time =", t
      print "Last Cd   =", cd
    }
  }
' "$CD_DAT"

echo

# ----------------------------
# 6) Plot Cd PDF
# ----------------------------
if command -v gnuplot >/dev/null 2>&1; then

  cat > "$PLOT_CD_GP" <<EOF
set terminal pdfcairo size 12in,7in font ",14"
set output '$(basename "$CD_PDF")'

set grid
set xlabel 't [s]'
set ylabel 'C_d [-]'
set title 'Oscillating sphere: drag coefficient vs time'
set key at graph 0.02, 0.98 left top
set border lw 1.5
set datafile commentschars "#"
set yrange [${YRANGE_MIN}:${YRANGE_MAX}]

plot \
  '$(basename "$CD_DAT")' using 1:2 with lines lw 2 title 'Present immersed-boundary', \
  '$(basename "$REF_SORTED_DAT")' using 1:2 with points pt 6 ps 0.6 title 'Gilmanov & Sotiropoulos'

set output
EOF

  ( cd "$OUTDIR" && gnuplot "$(basename "$PLOT_CD_GP")" )

  echo "Plot written:"
  echo "  $CD_PDF"
else
  echo "NOTE: gnuplot not found. Data files were generated but plot was not."
fi

echo
echo "Done. Contents of $OUTDIR:"
ls -lh "$OUTDIR"
