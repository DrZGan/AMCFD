# Task: Fix Temperature Field Accuracy Outside Local Region

## Problem

When `localnum > 0`, the simulation alternates between local and global enthalpy solves (e.g., `localnum=3` means 3 local steps then 1 global step). During local steps, only cells within `(ilo:ihi, jlo:jhi, klo:khi)` are solved. Cells **outside** the local region are completely frozen for `localnum` consecutive time steps — no conduction, no boundary heat loss, no enthalpy update. When the global step finally runs, those cells have effectively skipped `localnum * delt` of physical cooling time, producing inaccurate temperatures (too hot) in the far field.

## Root Cause

In `main.f90`, every solver call (`properties`, `bound_enthalpy`, `discretize_enthalpy`, `source_enthalpy`, `solution_enthalpy`, `enthalpy_to_temp`) receives `(ilo, ihi, jlo, jhi, klo, khi)`. During local steps, these bounds exclude the outer region entirely. No compensation is applied when any step (local or global) later solves those cells.

## Scope

All changes are in `fortran_new/`. Do NOT modify `fortran_benchmark/`. Do NOT modify any code in `fortran_new/` that is unrelated to this task.

---

## Solution: Skipped-Timestep Compensation via Effective delt

Instead of rewriting an explicit enthalpy update, leverage the existing implicit solver. The idea:

1. Track how many timesteps each cell has been skipped (frozen during local steps).
2. When a cell is **next solved** (whether in a local or global step), use an effective `delt_eff = (n_skipped + 1) * delt` so the implicit solver naturally compensates for the accumulated time.

The key discretization line is in `mod_discret.f90`:
```fortran
apnot(i,j,k) = den(i,j,k) / delt * volume(i,j,k)
su(i,j,k) = apnot(i,j,k) * hnot(i,j,k)
```

By replacing `delt` with `delt_eff(i,j,k)`, the transient term correctly represents the actual elapsed time since the cell was last solved. The implicit solver handles all the physics (conduction, boundary loss) naturally — no new equations needed.

### Key insight: the local region moves every timestep

The local region follows the laser beam, so its bounds `(ilo:ihi, jlo:jhi, klo:khi)` change each step. This means:
- A cell that was **outside** the local region in step N may enter the local region in step N+1 (because the laser moved toward it). That cell has `n_skipped > 0` and must be compensated when it is solved in step N+1's local solve.
- A cell that was **inside** the local region in step N may fall outside in step N+1. It starts accumulating `n_skipped` from that point.
- Therefore, `delt_eff` compensation must apply to **every** solve step (both local and global), not just global steps.

---

### Task 1 — Create `n_skipped` array and management routines

**In `mod_local_enthalpy.f90`:**

1. Declare a 3D integer array `n_skipped(ni, nj, nk)`, initialized to 0.
2. Add `subroutine allocate_skipped(ni, nj, nk)` — allocates and zeros the array.
3. Add `subroutine update_skipped(ilo, ihi, jlo, jhi, klo, khi, is_local)`:
   - If `is_local`:
     - For all interior cells `(2:nim1, 2:njm1, 2:nkm1)` **outside** `(ilo:ihi, jlo:jhi, klo:khi)`: increment `n_skipped(i,j,k)` by 1.
     - For cells **inside** `(ilo:ihi, jlo:jhi, klo:khi)`: reset `n_skipped(i,j,k) = 0` (they were just solved with compensation).
   - If global (`is_local = .false.`):
     - Reset the entire `n_skipped` array to 0 (all cells were just solved).

---

### Task 2 — Create `delt_eff` array

**In `mod_local_enthalpy.f90`:**

1. Declare a 3D real array `delt_eff(ni, nj, nk)`.
2. Add `subroutine compute_delt_eff()`:
   - For each cell: `delt_eff(i,j,k) = delt * (n_skipped(i,j,k) + 1)`
   - Cells with `n_skipped = 0` get `delt_eff = delt` (normal behavior).
   - This must be called **before** the iteration loop each step, so that `discretize_enthalpy` and `source_enthalpy` use the correct effective timestep.

---

### Task 3 — Modify `discretize_enthalpy` and `source_enthalpy` to use `delt_eff`

**In `mod_discret.f90`:**

Change the enthalpy transient term from:
```fortran
apnot(i,j,k) = den(i,j,k) / delt * volume(i,j,k)
```
to:
```fortran
apnot(i,j,k) = den(i,j,k) / delt_eff(i,j,k) * volume(i,j,k)
```

**In `mod_sour.f90`:**

Change the latent heat source term from:
```fortran
volht = volume(i,j,k) * hlatnt * den(i,j,k) / delt
```
to:
```fortran
volht = volume(i,j,k) * hlatnt * den(i,j,k) / delt_eff(i,j,k)
```

**Important:** `delt_eff` must be accessible from these modules. Add `use local_enthalpy, only: delt_eff` or pass it through module scope.

---

### Task 4 — Integrate into `main.f90`

In the time loop:

1. **After** `call allocate_fields(...)`, add `call allocate_skipped(ni, nj, nk)`.

2. **After** `call get_enthalpy_region(...)` and the local/global decision, **before** the iteration loop:
   - Call `call compute_delt_eff()` — this sets `delt_eff` for the current step based on `n_skipped`.

3. **After** `end do iter_loop` (after convergence), **before** the `tnot/hnot` copy block:
   ```fortran
   call update_skipped(ilo, ihi, jlo, jhi, klo, khi, is_local)
   ```

No other changes to the main loop structure. The existing solver calls (`discretize_enthalpy`, `source_enthalpy`, etc.) automatically pick up `delt_eff`.

---

### Task 5 — Timing

Add a dedicated timer `t_skipped_mgmt` in `mod_timing.f90`:
- Accumulate time spent in `update_skipped` and `compute_delt_eff`.
- Report in the timing summary as "Skipped-step management".
- Add "Total time used" (wall-clock, as shown in output.txt: `hr m s` format) to `timing_report.txt`.

---

### Task 6 — Verification

All test scripts and analysis code go in this project folder (`projects/20260304_LOCAL_REGION/`), NOT in the source code directories.

**Test procedure:**

1. Create a bash script `projects/20260304_LOCAL_REGION/run_verification.sh` that:
   - Copies `fortran_new/inputfile/input_param.txt` to a backup.
   - Sets `localnum = 0` in `input_param.txt`, compiles and runs the simulation, copies `result/thermal_history.txt` and `result/timing_report.txt` to `projects/20260304_LOCAL_REGION/baseline_thermal_history.txt` and `baseline_timing_report.txt`. Also copies `result/thermal_history.png` to `projects/20260304_LOCAL_REGION/baseline_thermal_history.png` if it exists.
   - Sets `localnum = 3` in `input_param.txt`, compiles and runs the simulation, copies results to `projects/20260304_LOCAL_REGION/local3_thermal_history.txt` and `local3_timing_report.txt`. Also copies `result/thermal_history.png` to `projects/20260304_LOCAL_REGION/local3_thermal_history.png` if it exists.
   - Restores the original `input_param.txt` from backup.
   - **Important:** Both runs must use identical input parameters except for `localnum`. The input file from `fortran_new/` is the single source of truth — do NOT copy from `fortran_benchmark/`.

2. Create a Python script `projects/20260304_LOCAL_REGION/compare_results.py` that:
   - Reads `baseline_thermal_history.txt` and `local3_thermal_history.txt` from this project folder.
   - Computes the 4 quantitative metrics for each of the 10 monitoring points.
   - Reads `baseline_timing_report.txt` and `local3_timing_report.txt`, extracts key timing values.
   - Generates a comparison PNG (`thermal_history_comparison.png`) showing all 10 points: solid lines = baseline, dashed lines = local. Save in this project folder.
   - Writes the full comparison to `projects/20260304_LOCAL_REGION/local_region_verification.txt` (NOT to `result/`).

Write results in `projects/20260304_LOCAL_REGION/local_region_verification.txt` including:

**Thermal history quantitative comparison (`thermal_history.txt`):**

The file has columns: `time(s) T1 T2 T3 T4 T5 T6 T7 T8 T9 T10` corresponding to 10 monitoring points:
- P1  (1.0mm, 0.50mm, 0.695mm) — scan track 1, near start, surface
- P2  (2.0mm, 0.50mm, 0.695mm) — scan track 1, centre, surface
- P3  (3.0mm, 0.50mm, 0.695mm) — scan track 1, near end, surface
- P4  (1.0mm, 0.50mm, 0.655mm) — scan track 1, near start, 40um below surface
- P5  (2.0mm, 0.50mm, 0.655mm) — scan track 1, centre, 40um below surface
- P6  (3.0mm, 0.50mm, 0.655mm) — scan track 1, near end, 40um below surface
- P7  (2.0mm, 0.25mm, 0.695mm) — offset from scan, surface
- P8  (2.0mm, 0.75mm, 0.695mm) — offset from scan, surface
- P9  (2.0mm, 1.50mm, 0.695mm) — far from scan, surface (substrate ref)
- P10 (2.0mm, 0.50mm, 0.200mm) — deep substrate below track 1

For each of the 10 points, report:
1. **Max absolute temperature difference** `max|T_local(t) - T_global(t)|` across all timesteps
2. **Max relative difference** `max|T_local(t) - T_global(t)| / T_peak_global` where `T_peak_global` is the peak temperature at that point in the baseline run
3. **Relative RMSD** `sqrt(mean(((T_local(t) - T_global(t)) / T_global(t))^2))` — RMS of pointwise relative differences
4. **Peak temperature difference** `|T_peak_local - T_peak_global|` — the error in the maximum temperature reached

Present as a table with columns: Point | MaxAbsDiff(K) | MaxRelDiff(%) | RelRMSD(%) | PeakTempDiff(K)

Acceptance criteria:
- Max relative difference across all 10 points: `max(MaxRelDiff) < 5%`
- Relative RMSD across all 10 points: `max(RelRMSD) < 3%`
- If not met, investigate and document the cause.

**Important constraints:**
- Do NOT modify the local region parameters: `local_half_x=1.0e-3, local_half_y=2.0e-4, local_depth_z=2.0e-4`
- For faster testing, set `timax=0.005` first. After passing, re-run with `timax=0.01` to confirm.

**Timing comparison (`timing_report.txt`):**
- Total wall time for localnum=0 vs localnum=3
- Per-module breakdown comparison
- Speedup ratio. The point of local solver is performance — if there is no meaningful speedup, something is wrong.

**Conclusion:** whether compensation is working (accuracy preserved + speedup achieved)

---

### Task 7 — Parametric Study: localnum Sensitivity (run after Task 6 passes)

**Prerequisite:** Task 6 with `localnum=3` must pass all acceptance criteria before starting this task.

**Objective:** Characterize how `localnum` (0, 4, 5, 6, 7, 8, 9, 10) affects accuracy and performance, to guide users in choosing the optimal value.

**Important:** Use `timax=0.01` for all parametric runs (not 0.005). Include `localnum=0` as baseline in the parametric sweep — do NOT reuse Task 6 baseline (which used `timax=0.005`).

**Execution strategy:**

- Single input file (`fortran_new/inputfile/input_param.txt`) — do NOT generate multiple input files
- `run.sh` modifies `case_name` in the input file via sed, then launches `./cluster_main`
- `cluster_main` reads the input file immediately at startup and closes it — it does not access the file again
- Parallel execution: before each launch, sed modifies `localnum`, then `run.sh` is called (which seds `case_name` and starts the process). A `sleep 2` between launches ensures the previous case has finished reading the file
- 24 cores / 6 threads per case = 4 cases per batch, 2 batches for 8 cases
- Compile once with `bash compile.sh` before running, then only call `run.sh`

**Execution steps:**

```bash
# 1. Compile
cd fortran_new && bash compile.sh

# 2. Batch 1: localnum=0,4,5,6 (4 cases x 6 threads = 24 cores)
for N in 0 4 5 6; do
    CASE="param_local${N}"
    [ "$N" -eq 0 ] && CASE="param_baseline"
    sed -i "s/localnum=[0-9]*/localnum=$N/" inputfile/input_param.txt
    bash run.sh "$CASE" 6 &
    sleep 2
done
wait

# 3. Batch 2: localnum=7,8,9,10 (4 cases x 6 threads = 24 cores)
for N in 7 8 9 10; do
    CASE="param_local${N}"
    sed -i "s/localnum=[0-9]*/localnum=$N/" inputfile/input_param.txt
    bash run.sh "$CASE" 6 &
    sleep 2
done
wait

# 4. Restore localnum=0
sed -i "s/localnum=[0-9]*/localnum=0/" inputfile/input_param.txt
```

Results are saved automatically in `result/<case_name>/<case_name>_thermal_history.txt` and `result/<case_name>/<case_name>_timing_report.txt`.

**Scripts (all in `projects/20260304_LOCAL_REGION/`):**

1. `compare_parametric.py`:
   - Read baseline and all `local{N}` thermal_history / timing_report from `fortran_new/result/param_*/`.
   - For each `localnum` value, compute the same 4 metrics (MaxAbsDiff, MaxRelDiff, RMSD, PeakTempDiff) for all 10 points.
   - Extract total wall time and speedup ratio for each `localnum`.
   - Generate `parametric_report.txt` containing:
     - **Summary table:** localnum | Speedup | P1-MaxAbsDiff | P1-RMSD | ... | P10-MaxAbsDiff | P10-RMSD
     - **Trend analysis:** for each point, how MaxAbsDiff and RMSD grow with `localnum`
     - **Recommended localnum range:** the largest `localnum` that still meets the acceptance criteria (max MaxRelDiff < 5%, max relative RMSD < 3%), with the best speedup
     - **Conclusion:** tradeoff summary — accuracy vs performance

---

## Performance

- Optimize OpenMP parallelization in `mod_local_enthalpy.f90` (`update_skipped`, `compute_delt_eff`) and in the hnot update block in `main.f90` to improve speedup ratio.
- The hnot/tnot/fraclnot update and enthalpy restoration loops in `main.f90` should use `!$OMP PARALLEL DO` for better performance.

## Constraints

- Do NOT relax convergence criteria for local steps. The convergence tolerance (`conv_res_heat`, `conv_res_cool`) must remain identical for both local and global steps.

## Notes

- This approach reuses the existing implicit solver entirely — no new physics equations or explicit updates.
- The `delt_eff` compensation applies to **all** solve steps (local and global), because the local region moves with the laser. A cell entering the local region after being skipped must also be compensated.
- `hnot` already stores the enthalpy from the last time the cell was actually solved, so the `apnot * hnot` source term naturally represents the correct old-time value.
- For cells with `n_skipped = 0`, `delt_eff = delt` — behavior is identical to the current code.
- No new input parameters are needed. All existing parameters are reused.

## Permissions

- You may modify files in `fortran_new/` that are **directly related to this task** (e.g., `mod_local_enthalpy.f90`, `mod_discret.f90`, `mod_sour.f90`, `main.f90`, `mod_timing.f90`). Do NOT touch unrelated code.
- You may create new files in `fortran_new/` if needed for the implementation.
- You may temporarily modify `fortran_new/inputfile/input_param.txt` for testing (must restore afterward).
- You may create scripts, data files, and results in `projects/20260304_LOCAL_REGION/`.
- Do NOT modify anything in `fortran_benchmark/`.
- Do NOT write test/verification scripts or result files into `fortran_new/` or `result/` — all test artifacts go in `projects/20260304_LOCAL_REGION/`.
- Do NOT ask for manual confirmation — proceed with implementation directly.

## Allowed Bash Commands

The following commands may be executed without user confirmation:

**Build & Run (working directory: `fortran_new/`):**

`compile.sh` already handles clean, compile, and run (it calls `clean.sh`, compiles all modules, links, and executes `./cluster_main`). So a single command is sufficient:

- `cd /mnt/d/Fortran/AMCFD/fortran_new && bash compile.sh` — clean + compile + run (all-in-one)
- `cd /mnt/d/Fortran/AMCFD/fortran_new && bash clean.sh` — clean build artifacts only
- `cd /mnt/d/Fortran/AMCFD/fortran_new && ./cluster_main` — run simulation only (if already compiled)

**Test & Analysis:**
- `bash projects/20260304_LOCAL_REGION/run_verification.sh` — run verification tests
- `bash projects/20260304_LOCAL_REGION/run_parametric.sh` — run parametric study (localnum=4–10)
- `python3 projects/20260304_LOCAL_REGION/compare_results.py` — generate comparison report
- `python3 projects/20260304_LOCAL_REGION/compare_parametric.py` — generate parametric report
- `python3 <any script in projects/20260304_LOCAL_REGION/>` — run any analysis script in this project

**File Operations:**
- `ls`, `cat`, `head`, `tail`, `wc` — inspect files and results
- `diff` — compare files
- `grep`, `find` — search code
- `cp`, `mv` — copy/move files within `fortran_new/` or `projects/20260304_LOCAL_REGION/`
- `rm` — remove temporary files in `projects/20260304_LOCAL_REGION/`
- `mkdir -p` — create directories
- `chmod +x` — make scripts executable
- `sed` — for modifying `localnum` in `input_param.txt` during test runs
