# Master Thesis

This repository contains my notes and Fortran experiments on a priori pattern selection for energy-minimizing AMG prolongation.

My main reference is *Energy Optimization of Algebraic Multigrid Bases* by Mandel, Brezina and Vaněk. I currently compare coupling strength, connection score and energy score. Each method selects positions from the same candidate set. For each selected pattern, I use Algorithm 3 to minimize the interpolation energy.

## Notes and code

- [Paper notes](notes/notes_on_paper.md) contain my understanding of the paper and the score derivations.
- [Experiment notes](notes/result_of_scalar_experiments.md) contain the scalar results and my explanations.
- [Fortran source](experiments/src) contains the prototype, including the 1D example and the rail371 example.
- [Results](experiments/result) contains the saved scores, energies and comparison figure.

## Running the experiments

I use GNU Fortran and CMake 3.20 or later. Run these commands from the repository root.

```sh
cmake -S experiments -B build
cmake --build build
./build/amg_test
```

For the rail371 example, download `fenics_rail_371.zip` from [Zenodo](https://doi.org/10.5281/zenodo.5113560) and put the extracted matrix files in `experiments/data/fenics_rail_371/`. The data is downloaded separately.

The test case and score are selected by enabling the corresponding calls in [main.f90](experiments/src/main.f90). Update `run_name` to match the selection and rebuild after making changes.

The 1D case uses an exhaustive reference and prints the energy comparison in the terminal. The rail371 case writes its score and energy tables to `experiments/result/`, with filenames based on `run_name`.
