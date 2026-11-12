# AMG experiments

These Fortran experiments compare three scores for selecting the sparsity pattern of AMG prolongation. For each selected pattern, Algorithm 3 minimizes the interpolation energy.

## Choose the case and score

Edit [src/main.f90](src/main.f90). Enable one case setup call and comment out the other:

- `setup_scalar_1d_case`: the weighted 1D diffusion example with 12 internal points.
- `setup_scalar_2d_case`: the rail371 example with $A = S + M_\Gamma$.

Then enable only one of the three score calls:

- `compute_coupling_strength`
- `compute_connection_score`
- `compute_energy_score`

Set `run_name` to match the case and score, for example `rail371_energy`. This name only controls the output filenames. The enabled calls select the case and score.

## Data for rail371

Download `fenics_rail_371.zip` from [Zenodo](https://doi.org/10.5281/zenodo.5113560). Put these two files from the archive in `experiments/data/fenics_rail_371/`, relative to the repository root:

- `ODE_371_S.mtx`
- `ODE_371_M_GAMMA.mtx`

The program reads both matrices and adds them to build $A$.

## Build and run

Use GNU Fortran and CMake 3.20 or later. Run these commands from the repository root:

```sh
cmake -S experiments -B build
cmake --build build
./build/amg_test
```

Rebuild after changing the case, score or settings in `main.f90`. Run the program from the repository root so that it can find the data and output directories.

## Output

Results are saved in `experiments/result/`, relative to the repository root. Keep this directory in place before running.

- `<run_name>_scores.txt` is written for every run. Its columns are `rank candidate_id row col score`. The rows follow the ranking; `candidate_id` is the original candidate number, `row` is the fine node and `col` is the coarse column.
- If the number of candidates is at most `exhaustive_limit` (currently 12), the program computes an exhaustive reference. It prints `failed solves` and the table `budget ranked_J best_J gap` in the terminal. This is the branch used by the 1D example; it does not write an energy file.
- For larger candidate sets, including rail371, exhaustive enumeration is skipped. The program writes `<run_name>_energy.txt` with the columns `budget ranked_J`.

The budget is the number of added candidate positions. Running again with the same `run_name` replaces the corresponding output files.

## Data citation

Saak, J., & Behr, M. (2021). *The Oberwolfach Steel-Profile benchmark revisited* (Version 1.0) [Data set]. [Zenodo](https://doi.org/10.5281/zenodo.5113560).
