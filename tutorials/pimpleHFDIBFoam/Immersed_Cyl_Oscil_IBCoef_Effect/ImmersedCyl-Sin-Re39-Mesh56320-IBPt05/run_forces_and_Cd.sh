#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Extract Fx from bodiesInfo/<time>/body<ID>.info (FCoupling_F)
# and compute Cd = Fx / (0.5*rho*Uref^2*Aref).
# Also plots local reference data from:
#   ./Cd.dat
# ============================================================

# ----------------------------
# User settings
# ----------------------------
BODY_ID=0
RHO=1.0
D=0.1
SPAN=0.1
UREF=0.3926990817

# Exact local reference file in this case folder
REF_CD_FILE="./Cd.dat"

# Optional y-range (leave empty to autoscale)
YRANGE_MIN=${YRANGE_MIN:--4.5}
YRANGE_MAX=${YRANGE_MAX:-4.5}

# ----------------------------
# Paths
# ----------------------------
BODIES_DIR="bodiesInfo"
OUTDIR="postProcessing/forces"

FORCES_DAT="$OUTDIR/forces_body${BODY_ID}.dat"
DRAG_DAT="$OUTDIR/dragFx_body${BODY_ID}.dat"
CD_DAT="$OUTDIR/Cd_body${BODY_ID}.dat"
REF_SORTED_DAT="$OUTDIR/reference_Cd_sorted.dat"

PLOT_DRAG_GP="$OUTDIR/plot_dragFx_body${BODY_ID}.gp"
PLOT_CD_GP="$OUTDIR/plot_Cd_body${BODY_ID}.gp"

DRAG_PDF="$OUTDIR/dragFx_body${BODY_ID}.pdf"
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
  exit 1
fi

mkdir -p "$OUTDIR"

# ----------------------------
# Derived values
# ----------------------------
AREF=$(awk -v D="$D" -v s="$SPAN" 'BEGIN{printf "%.15g", D*s}')
DEN=$(awk -v rho="$RHO" -v u="$UREF" -v a="$AREF" 'BEGIN{printf "%.15g", 0.5*rho*u*u*a}')

echo "Using:"
echo "  BODY_ID = $BODY_ID"
echo "  rho     = $RHO"
echo "  D       = $D"
echo "  span    = $SPAN"
echo "  Aref    = $AREF"
echo "  Uref    = $UREF"
echo "  den     = $DEN"
echo "  refCd   = $REF_CD_FILE"
echo

# ----------------------------
# 1) Build forces_body*.dat: time Fx Fy Fz
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
# 2) Drag-only: time Fx
# ----------------------------
awk '
  BEGIN{print "# time Fx"}
  $1 ~ /^#/ {next}
  NF>=2 {print $1, $2}
' "$FORCES_DAT" > "$DRAG_DAT"

echo "Wrote: $DRAG_DAT"

# ----------------------------
# 3) Cd: time Cd
# ----------------------------
awk -v den="$DEN" '
  BEGIN{print "# time Cd"}
  $1 ~ /^#/ {next}
  NF>=2 {print $1, $2/den}
' "$DRAG_DAT" > "$CD_DAT"

echo "Wrote: $CD_DAT"

# ----------------------------
# 4) Sort local reference Cd.dat
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
# 6) Plot PDFs with gnuplot
# ----------------------------
if command -v gnuplot >/dev/null 2>&1; then

  cat > "$PLOT_DRAG_GP" <<EOF
set terminal pdfcairo size 12in,7in font ",12"
set output '$(basename "$DRAG_PDF")'
set grid
set xlabel 'time [s]'
set ylabel 'Fx from bodiesInfo (FCoupling_F.x)'
set title 'Drag force Fx vs time (body${BODY_ID})'
set key left top
set datafile commentschars "#"
plot '$(basename "$DRAG_DAT")' using 1:2 with lines lw 2 title 'Fx'
set output
EOF

  cat > "$PLOT_CD_GP" <<EOF
den = $DEN
set terminal pdfcairo size 12in,7in font ",14"
set output '$(basename "$CD_PDF")'
set grid
set xlabel 't [s]'
set ylabel 'C_d [-]'
set title 'Drag Coefficient (C_d) vs. Time (t)'
set key at graph 0.02, 0.98 left top
set border lw 1.5
set datafile commentschars "#"
set yrange [${YRANGE_MIN}:${YRANGE_MAX}]
plot \
  '$(basename "$CD_DAT")' using 1:2 with lines lw 2 title 'Present', \
  '$(basename "$REF_SORTED_DAT")' using 1:2 with points pt 6 ps 0.5 title 'Wan&Turek'
set output
EOF

  ( cd "$OUTDIR" && gnuplot "$(basename "$PLOT_DRAG_GP")" )
  ( cd "$OUTDIR" && gnuplot "$(basename "$PLOT_CD_GP")" )

  echo "Plots written:"
  echo "  $DRAG_PDF"
  echo "  $CD_PDF"
else
  echo "NOTE: gnuplot not found. Data files were generated but plots were not."
fi

echo
echo "Done. Contents of $OUTDIR:"
ls -lh "$OUTDIR"
