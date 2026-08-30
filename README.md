# 3D-Printed Ku-Band GRIN Lens Synthesizer

This repository contains a workflow to design, validate, and simulate a 3D-printed Gradient-Index (GRIN) phase-correction lens optimized for a 120 mm x 80 mm horn antenna.

It calculates the required Fermat phase delays across the antenna aperture, converts those delays into a physical 2D hole-punch pattern using Maxwell-Garnett effective medium theory, and builds the 3D model in Ansys HFSS via PyAEDT.

## Workflow

### 1. Phase Map & Fill Fraction (MATLAB)
`phase_to_hole_diameter.m` calculates the physical print geometry.
* Rebuilds the spherical-wave phase error map based on the horn's E-plane and H-plane phase centers.
* Inverts the phase map to calculate the required physical delay.
* Uses the 2D Maxwell-Garnett mixing formula to determine the required air-fill fraction for a given base material (e.g., ABS, $\epsilon_r = 2.6$).
* Maps the fill fraction to exact hole diameters on a square lattice.
* **Output:** `lens_hole_table.csv` containing the $(x, y)$ coordinates and diameter of every hole.

### 2. HFSS Model Construction (Python)
Two versions of `main.py` are provided to accommodate different HFSS license constraints. 

* **Student License (Homogenized):** Bins the continuous hole data into coarse structural tiles to bypass the 64,000 mesh element limit. Calculates the averaged effective permittivity for each tile and draws the lens as a grid of solid blocks with dynamically generated materials.
* **Commercial License (Explicit):** Constructs the exact physical geometry by creating a solid dielectric slab and executing a Boolean subtraction for every individual cylindrical hole. Used for rigorous cross-polarization and near-field analysis prior to manufacturing.

## Prerequisites
* MATLAB
* Python 3.7+
* PyAEDT (`pip install pyaedt`)
* Ansys Electronics Desktop (AEDT)

## Usage
1. Run `phase_to_hole_diameter.m` to generate `lens_hole_table.csv`.
2. Place the CSV in the same directory as the selected Python script.
3. Execute `main.py`. Ansys HFSS will open automatically and build the geometry.

*Note: When using specialty RF materials (e.g., PREPERM, HIPS), update the base permittivity and loss tangent variables in both scripts to match measured coupon data.*