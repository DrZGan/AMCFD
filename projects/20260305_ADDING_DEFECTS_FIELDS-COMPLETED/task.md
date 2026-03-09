# Task: Adding Defects Fields

## Objective
Create a new `defect.f90` module that computes and prints a defect field array (similar to `temp`).
Minimize impact on other modules unless strictly necessary.
The module supports multiple defect computation methods, each with its own parameters defined at the top of `defect.f90`.

---

## Method 1: `maxtemp_determ` (Deterministic Max-Temperature Method)

### Parameters (defined at top of `defect.f90`, labeled as `maxtemp_determ` parameters)
- `k_lof` : calibration parameter for lack-of-fusion, default = 1.0, typical range 0–1
- `k_kep` : calibration parameter for keyhole porosity scaling

---

### Task 1 — Create `max_temp` Array
- Allocate `max_temp` and `defect` arrays covering the **full X-Y control volume**, but only from `kmax` down to `layerheight` in Z (i.e., Z extent = one layer height). Defect computation is limited to this single-layer region.
- At every time step, for each cell within this Z range: if the current temperature exceeds the recorded `max_temp`, update `max_temp` with the current temperature.

---

### Task 2 — Compute Defect Array (call after simulation completes)
- `defect` array is already allocated (see Task 1), initialized to `0`.
- Apply the following rules cell by cell using `max_temp`:
  - If `max_temp < T_solid` → `defect = -k_lof` (lack-of-fusion)
  - If `max_temp > T_boiling` → `defect = k_kep * (max_temp - T_boiling) / T_boiling` (keyhole porosity)
  - If `T_solid ≤ max_temp ≤ T_boiling` → `defect = 0`
- **Clean-up step**: Read the toolpath file to determine the maximum laser scanning range in the X-Y plane. Set `defect = 0` for all cells outside this range.

---

### Task 3 — Placeholder Method `maxtemp_stochas`
- Create an empty subroutine `maxtemp_stochas` in `defect.f90` as a future stochastic alternative.
- The simulation workflow should allow switching between `maxtemp_determ` and `maxtemp_stochas`.

---

### Task 4 — Timing (`mod_timing.f90`)
- Add a dedicated timer for the `max_temp` update step (called every time step).
- Record and accumulate the total time spent updating `max_temp`.
- Report this in the timing summary.

---

### Task 5 — Output (`mod_print.f90`)
After the defect array is computed:

**Arrays to output:**
- `max_temp` field
- `defect` field

**Defect metrics** (computed over the volume defined by max laser scanning range × layer height):

| Metric | Definition |
|--------|-----------|
| Defect fraction | `Σ abs(defect(i)) × V_cell(i) / V_total` |
| Lack-of-fusion fraction | `Σ abs(defect(i)) × V_cell(i) / V_total` where `defect(i) < 0` |
| Keyhole porosity fraction | `Σ defect(i) × V_cell(i) / V_total` where `0 < defect(i) ≤ 1` |

> Note: `defect fraction = lack-of-fusion fraction + keyhole porosity fraction`

**`defect_report.txt`** (written to `results/`), containing:
- Defect fraction
- Lack-of-fusion fraction
- Keyhole porosity fraction
- Max laser scanning range (node coordinates, X-Y plane boundary)
- Defect volume
- Lack-of-fusion volume
- Keyhole porosity volume
- Total reference volume (max scan range × layer height)

---

## Notes

