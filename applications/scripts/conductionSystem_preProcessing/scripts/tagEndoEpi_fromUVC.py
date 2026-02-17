#!/usr/bin/env python3
"""Tag endocardial/epicardial surfaces from UVC fields.

Inputs:
- Legacy ASCII VTK (volume or surface) with POINT_DATA:
  - uvc_transmural
  - uvc_intraventricular

Outputs:
- outputs/<input>_endo_epi_surface.vtk (surface tags, unstructured)
- outputs/<input>_complementary_endocardium_epicardium_label.vtk (optional)
- outputs/<input>_endo_epi_volume.vtk (optional volume-mapped tags)

Main:
- Thin wrapper around conduction_preproc.tagging.tag_endo_epi_surface.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from conduction_preproc.tagging.tag_endo_epi_surface import main  # noqa: E402


if __name__ == "__main__":
    main()
