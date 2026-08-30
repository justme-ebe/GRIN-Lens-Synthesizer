"""
main.py (Explicit Version)
Builds the exact physical GRIN lens geometry in HFSS via PyAEDT.
Performs Boolean subtractions for all drilled holes. 
"""

import numpy as np
import pandas as pd
from ansys.aedt.core import Hfss

# Configuration
CSV_PATH  = "lens_hole_table.csv"
T_LENS_MM = 48.2
EPS_ABS   = 2.6
TAND_ABS  = 0.008
AEDT_VER  = "2025.2"
STUDENT   = False 

# 1. Load Data
df = pd.read_csv(CSV_PATH)
df = df[df["d_mm"] > 1e-4].reset_index(drop=True)

x_unique = np.sort(df["x_mm"].unique())
pitch_mm = float(np.median(np.diff(x_unique)))

x_min, x_max = df["x_mm"].min(), df["x_mm"].max()
y_min, y_max = df["y_mm"].min(), df["y_mm"].max()

lens_x_length = (x_max - x_min) + pitch_mm
lens_y_length = (y_max - y_min) + pitch_mm
lens_x_origin = x_min - (pitch_mm / 2.0)
lens_y_origin = y_min - (pitch_mm / 2.0)

# 2. Initialize HFSS
hfss = Hfss(version=AEDT_VER, student_version=STUDENT, new_desktop=True,
            project="GRIN_Lens_Commercial", design="lens_explicit")
hfss.modeler.model_units = "mm"

mat_name = "ABS_Solid"
if mat_name not in hfss.materials.mat_names_aedt:
    mat = hfss.materials.add_material(mat_name)
    mat.permittivity = EPS_ABS
    mat.dielectric_loss_tangent = TAND_ABS

# 3. Construct Geometry
lens_body = hfss.modeler.create_box(
    origin=[lens_x_origin, lens_y_origin, 0],
    sizes=[lens_x_length, lens_y_length, T_LENS_MM],
    name="Lens_Solid_Body",
    material=mat_name,
)

tool_list = []
for i, row in df.iterrows():
    radius = row["d_mm"] / 2.0
    hole_cyl = hfss.modeler.create_cylinder(
        cs_axis="Z",
        position=[row["x_mm"], row["y_mm"], 0],
        radius=radius,
        height=T_LENS_MM,
        name=f"hole_{i}",
        material="vacuum"
    )
    tool_list.append(hole_cyl.name)

hfss.modeler.subtract(
    blank_list=[lens_body.name], 
    tool_list=tool_list, 
    keep_originals=False
)

hfss.save_project()