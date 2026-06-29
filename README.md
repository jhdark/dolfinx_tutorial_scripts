# DOLFINx tutorial scripts

A set of small, self-contained finite-element examples to help new fellows get started
with **DOLFINx** and **FESTIM**. Each example introduces one new concept on top of the
previous one.

## The examples

| Folder | Topic | New concept |
| --- | --- | --- |
| [`Example_1`](Example_1) | Steady state heat transfer | The anatomy of a DOLFINx script |
| [`Example_2`](Example_2) | Transient heat transfer | Time stepping |
| [`Example_3`](Example_3) | Two species | Multiple fields (mixed function space) |
| [`Example_4`](Example_4) | Hydrogen trapping | Coupling fields through a reaction |

Work through them in order; Example 4 is the basic FESTIM model.

## Building the environment

The examples need a conda environment with DOLFINx 0.10. From this folder:

```bash
conda env create -f environment.yml
conda activate dolfinx-tutorial-env
```

(`mamba env create -f environment.yml` works too and is faster.)

## How to work through an example

Each `Example_N` folder contains the same three things:

1. **`eg_N.ipynb`**: a notebook that explains the problem step by step. **Read this
   first.**
2. **`my_script.py`**: a near-empty scaffold for **you** to fill in. After reading the
   notebook, try to reproduce the solution here from scratch. This is where the
   learning happens.
3. **`eg_N.py`**: the complete reference script. Check your version against it once you
   have had a go, or use it if you get stuck.

Each folder also has its own `README.md` with the full problem statement.

## Viewing results

Example 1 plots its result with PyVista. Examples 2 to 4 write `.bp` files; open these
in [ParaView](https://www.paraview.org/) to scrub through time.
