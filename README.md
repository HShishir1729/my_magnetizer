# Magnetic Field Profile Analysis Pipeline

## Overview

This repository contains scripts and configuration used to analyze magnetic field profiles from galaxy simulations. The workflow focuses on:

* Extracting radial magnetic field profiles (`B_\phi`)
* Matching profiles between raw `.dat` outputs and HDF5 snapshots
* Matching seed field profiles to evolved profiles
* Visualizing the evolution and seed placement across redshifts

This pipeline is particularly useful for identifying **where seed magnetic fields are injected between GALFORM snapshots**.

---

## 1. Global Configuration (Fortran)

File: `global_input_parameter.f90`

This file controls which physical quantities are written to `.dat` files during simulation runs.

### Relevant Flags

```fortran
logical :: write_bfloor   = .true.
logical :: write_afloor   = .false.
logical :: write_bp       = .true.
logical :: write_r        = .true.
logical :: write_seed     = .true.
logical :: write_seed_r   = .true.
logical :: write_seed_beq = .true.

character(len=100) :: output_suffix = '_fid_10times_random'
```

### Description

Each logical flag enables/disables output:

* `write_bp` → azimuthal magnetic field profiles (`Bp`)
* `write_r` → corresponding radial grid
* `write_seed` → seed magnetic field profiles
* `write_seed_r` → seed radii
* `write_seed_beq` → equipartition field at seed stage
* `write_bfloor`, `write_afloor` → floor-related quantities

> These outputs are written as `.dat` files and later used in Python scripts.

---

## 2. Profile Matching: Raw vs HDF5

Script: `profile_matching_corrected.py`

### Purpose

Matches profiles from:

* Raw `.dat` files (direct outputs)
* HDF5 simulation snapshots

This ensures consistency and helps identify which snapshot corresponds to which raw profile.

### Key Steps

1. **Input**

   ```bash
   python profile_matching_corrected.py galaxy_list.txt TASK_ID
   ```

2. **Read Data**

   * Variable-length `.dat` profiles (`Bp`, `r`)
   * HDF5 datasets:

     * `Output/Bp`
     * `Output/r`

3. **Preprocessing**

   * Select valid region between zero crossings
   * Sort profiles radially
   * Normalize radius using:

     ```
     r_disk = max(r) / 2.7
     ```

4. **Sampling**
   Profiles are interpolated at fixed disk fractions:

   ```
   r = r_disk × [1/4, 2/4, 3/4, 4/4, 5/4]
   ```

5. **Matching Algorithm**

   * Adaptive tolerance matching (`np.allclose`)
   * Fallback: partial matching (≥3 points)
   * Best match chosen via RMSE

6. **Outputs**
   Saved in:

   ```
   results/files_corrected/
   ```

   * Sampled profiles (raw & HDF5)
   * Row-to-row mapping:

     ```
     row_matches_reduced_<gal>.txt
     ```

---

## 3. Seed Matching: Seed vs Raw Profiles

Script: `seed_matching_corrected2.py`

### Purpose

Matches **seed magnetic field profiles** to raw profiles using a **global optimization approach**.

This is crucial to determine **where the seed enters the galaxy evolution sequence**.

### Key Steps

1. **Read Inputs**

   * Seed profiles (`seed_output_*`)
   * Seed radii (`radius_seed_*`)
   * Raw profiles (`Bp_profiles_*`, `r_profiles_*`)

2. **Sampling**
   Same disk-fraction interpolation as before:

   ```
   r = r_disk × [1/4, 2/4, 3/4, 4/4, 5/4]
   ```

3. **Stability Handling**

   * Values rounded to 7 decimal places
   * Removes numerical noise

4. **Filtering**

   * Removes invalid (all-zero) rows

5. **Cost Function**
   Combines:

   * Shape similarity (normalized Bp)
   * Radial consistency

   ```
   cost = bp_error + 0.3 × r_error
   ```

6. **Global Matching**
   Uses Hungarian Algorithm:

   ```python
   linear_sum_assignment()
   ```

7. **Outputs**

   * Sampled seed and raw profiles
   * Matching index map:

     ```
     seed_row_matches_<gal>.txt
     ```

---

## 4. Visualization: Overlay Plot

Script: `overlay_plot3.py`

### Purpose

Generates a **multi-page PDF** showing:

* Evolution of magnetic field profiles across redshift
* Comparison between:

  * Fiducial run
  * 100× random seed run
* Seed field and equipartition field overlay

---

### Features

#### 1. Profile Evolution

* Plots `B_\phi(r)` for all redshifts
* Radius normalized by:

  ```
  r / r_half  (r_half = max(r)/2.7)
  ```

#### 2. Two Simulation Comparisons

* Solid line → Fiducial
* Dotted line → 100× random

#### 3. Seed Field Overlay

* Plotted at **z = 6**
* Includes:

  * Seed field
  * Equipartition field (`Beq/100`)

#### 4. Zoom Panel

* Focuses on seed amplitude scale
* Dynamically rescales y-axis

#### 5. Color Coding

* Redshift encoded via colormaps:

  * Blue → Fiducial
  * Orange → 100× random

---

### Output

```
galaxy_profiles_overlay_with_seed_zoom_corrected.pdf
```

Each page corresponds to one galaxy.

---

## 5. Full Workflow

### Step-by-Step Pipeline

1. **Run Simulation**

   * Enable outputs via `global_input_parameter.f90`

2. **Generate `.dat` Files**

   * Profiles: `Bp`, `r`
   * Seed: `seed`, `r_seed`, `seed_beq`

3. **Match Profiles (Raw ↔ HDF5)**

   ```bash
   python profile_matching_corrected.py galaxy_list.txt TASK_ID
   ```

4. **Match Seed to Profiles**

   ```bash
   python seed_matching_corrected2.py
   ```

5. **Visualize Results**

   ```bash
   python overlay_plot3.py
   ```

---
