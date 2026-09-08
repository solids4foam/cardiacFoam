#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# User settings
# ----------------------------
BODY_ID=0
RHO=1000.0
UINF=1.0
AREF=0.20295

# ----------------------------
# Paths
# ----------------------------
BODIES_DIR="bodiesInfo"
OUTDIR="postProcessing/forces"

FORCES_DAT="$OUTDIR/forces_body${BODY_ID}.dat"
DRAG_DAT="$OUTDIR/dragFx_body${BODY_ID}.dat"
CD_DAT="$OUTDIR/Cd_body${BODY_ID}.dat"

PLOT_DRAG_GP="$OUTDIR/plot_dragFx_body${BODY_ID}.gp"
PLOT_CD_GP="$OUTDIR/plot_Cd_body${BODY_ID}.gp"

DRAG_PDF="$OUTDIR/dragFx_body${BODY_ID}.pdf"
CD_PDF="$OUTDIR/Cd_body${BODY_ID}.pdf"

# Denominator for Cd
DEN=$(awk -v rho="$RHO" -v u="$UINF" -v a="$AREF" 'BEGIN{printf "%.15g", 0.5*rho*u*u*a}')
echo "Using: rho=$RHO, Uinf=$UINF, Aref=$AREF  => den=0.5*rho*U^2*Aref=$DEN"

# ----------------------------
# Checks
# ----------------------------
if [[ ! -d "$BODIES_DIR" ]]; then
  echo "ERROR: Cannot find '$BODIES_DIR' in case root."
  echo "Run this script from the case root where bodiesInfo/ exists."
  exit 1
fi

mkdir -p "$OUTDIR"

# ----------------------------
# 1) Build forces_body*.dat  (time Fx Fy Fz)
# ----------------------------
{
  echo "# time Fx Fy Fz"
  # sort numeric (handles 0.2, 1, 10, ...)
  for d in $(ls -1 "$BODIES_DIR" | sort -g); do
    f="$BODIES_DIR/$d/body${BODY_ID}.info"
    [[ -f "$f" ]] || continue

    # Extract vector from a line like:
    # FCoupling_F   (Fx Fy Fz)
    vec=$(awk '
      $1=="FCoupling_F"{
        gsub(/[()]/,"",$2);          # in case vector is split oddly (rare)
      }
    ' "$f" 2>/dev/null)

    # Robust extraction using grep+sed:
    # pulls whatever is inside parentheses on FCoupling_F line
    vec=$(grep -m1 "^FCoupling_F" "$f" | sed -E 's/.*\(\s*([^)]*)\s*\).*/\1/')

    # Expect 3 numbers
    if [[ -n "${vec:-}" ]]; then
      echo "$d $vec"
    fi
  done
} > "$FORCES_DAT"

# Quick sanity
echo "Wrote: $FORCES_DAT"
echo "Preview:"
head -n 5 "$FORCES_DAT" || true
echo "Tail:"
tail -n 5 "$FORCES_DAT" || true

# ----------------------------
# 2) Drag-only file (Fx): dragFx_body*.dat (time Fx)
# ----------------------------
awk 'NF==4 && $1!="#" {print $1, $2} END{}' "$FORCES_DAT" \
  | awk 'BEGIN{print "# time Fx"} {print $0}' > "$DRAG_DAT"

echo "Wrote: $DRAG_DAT"

# ----------------------------
# 3) Cd file: Cd_body*.dat (time Cd)
# ----------------------------
awk -v den="$DEN" '
  BEGIN{print "# time Cd"}
  NF==2 && $1!="#" {print $1, $2/den}
' "$DRAG_DAT" > "$CD_DAT"

echo "Wrote: $CD_DAT"

# ----------------------------
# 4) Stats + last Cd
# ----------------------------
awk '
  NF==4 && $1!="#" {
    Fx=$2
    if(n==0){min=Fx; max=Fx}
    if(Fx<min) min=Fx
    if(Fx>max) max=Fx
    sum+=Fx; n++
  }
  END{
    if(n>0){
      mean=sum/n
      amp=0.5*(max-min)
      print "---- Fx stats (from forces file) ----"
      print "N       =", n
      print "Fx_mean =", mean
      print "Fx_min  =", min
      print "Fx_max  =", max
      print "Fx_amp  =", amp
    } else {
      print "No valid force rows found in forces file."
    }
  }
' "$FORCES_DAT"

awk '
  NF==2 && $1!="#" {t=$1; cd=$2}
  END{
    if(t!=""){
      print "---- Last Cd ----"
      print "Last time =", t
      print "Last Cd   =", cd
    }
  }
' "$CD_DAT"

# ----------------------------
# 5) Plot PDFs with gnuplot (run inside OUTDIR so paths are simple)
# ----------------------------
if command -v gnuplot >/dev/null 2>&1; then
  cat > "$PLOT_DRAG_GP" <<EOF
set terminal pdfcairo size 12in,7in font ",12"
set output '$(basename "$DRAG_PDF")'
set grid
set xlabel 'time'
set ylabel 'Drag force Fx (FCoupling_F.x)'
set title 'Drag force Fx vs time (body${BODY_ID})'
set key left top
set datafile commentschars "#"
plot '$(basename "$DRAG_DAT")' using 1:2 with lines lw 2 title 'Fx'
set output
EOF

  cat > "$PLOT_CD_GP" <<EOF
den = $DEN
set terminal pdfcairo size 12in,7in font ",12"
set output '$(basename "$CD_PDF")'
set grid
set xlabel 'time'
set ylabel 'Cd = Fx / (0.5*rho*U^2*Aref)'
set title sprintf('Cd vs time (body${BODY_ID}), den=%.6g', den)
set key left top
set datafile commentschars "#"
plot '$(basename "$CD_DAT")' using 1:2 with lines lw 2 title 'Cd'
set output
EOF

  ( cd "$OUTDIR" && gnuplot "$(basename "$PLOT_DRAG_GP")" )
  ( cd "$OUTDIR" && gnuplot "$(basename "$PLOT_CD_GP")" )

  echo "Plots written:"
  echo "  $DRAG_PDF"
  echo "  $CD_PDF"
else
  echo "NOTE: gnuplot not found. Data files are ready, but plots were not generated."
fi

echo "Done. Contents of $OUTDIR:"
ls -lh "$OUTDIR"
