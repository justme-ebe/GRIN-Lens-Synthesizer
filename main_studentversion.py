"""
main.py (Homogenized Version)
Builds a homogenized (effective-permittivity) GRIN lens in HFSS via PyAEDT.
Optimized to stay under the 64,000 volume mesh element limit.
"""

import numpy as np
import pandas as pd
from ansys.aedt.core import Hfss

# Configuration
CSV_PATH  = "lens_hole_table.csv"
TILE_MM   = 5.0
T_LENS_MM = 48.2
EPS_ABS   = 2.6
TAND_ABS  = 0.008
AEDT_VER  = "2025.2"
STUDENT   = True

def mg_forward(f, eps_h=EPS_ABS, eps_i=1.0):
    return eps_h * ((eps_i + eps_h) + f * (eps_i - eps_h)) / ((eps_i + eps_h) - f * (eps_i - eps_h))

# 1. Load Data
df = pd.read_csv(CSV_PATH)
pitch_mm = float(np.median(np.diff(np.sort(df.x_mm.unique()))))
f_fill = (np.pi / 4.0) * (df.d_mm / pitch_mm)**2
df["eps_eff"] = mg_forward(f_fill)

# 2. Bin into coarse tile grid
df["ix"] = np.round(df.x_mm / TILE_MM).astype(int)
df["iy"] = np.round(df.y_mm / TILE_MM).astype(int)

tiles = (df.groupby(["ix", "iy"])
           .agg(eps_eff=("eps_eff", "mean"))
           .reset_index())

# Rigid snap to prevent HFSS part intersections
tiles["x_mm"] = tiles["ix"] * TILE_MM
tiles["y_mm"] = tiles["iy"] * TILE_MM

# 3. Build in HFSS
hfss = Hfss(version=AEDT_VER, student_version=STUDENT, new_desktop=True,
            project="GRIN_Lens_Student", design="lens_homogenized")
hfss.modeler.model_units = "mm"

for _, row in tiles.iterrows():
    tag = f"{row.eps_eff:.3f}".replace(".", "p")
    mat_name = f"ABS_eps_{tag}"
    
    if mat_name not in hfss.materials.mat_names_aedt:
        mat = hfss.materials.add_material(mat_name)
        mat.permittivity = float(row.eps_eff)
        mat.dielectric_loss_tangent = TAND_ABS * (row.eps_eff - 1.0) / (EPS_ABS - 1.0)

    hfss.modeler.create_box(
        origin=[row.x_mm - TILE_MM / 2, row.y_mm - TILE_MM / 2, 0],
        sizes=[TILE_MM, TILE_MM, T_LENS_MM],
        name=f"tile_{int(row.ix)}_{int(row.iy)}",
        material=mat_name,
    )

hfss.save_project()