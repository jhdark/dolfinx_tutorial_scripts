# Example 4: Hydrogen trapping in 1D

## Problem statement

This is the model at the heart of **FESTIM**: hydrogen transport with trapping. Two
fields live on the rod `Ω = (0, 1)` and now they are **coupled** by a reaction:

- `cm`: the **mobile** concentration, which diffuses through the material.
- `ct`: the **trapped** concentration, which is stuck in defects and does **not**
  diffuse.

Mobile hydrogen falls into empty traps and trapped hydrogen is released back, so the
two fields exchange material everywhere in the domain.

- The domain is the 1D interval `Ω = (0, 1)`.
- Diffusion coefficient `D = 0.1` (mobile field only).
- Trapping rate `k = 0.0001`, detrapping rate `p = 0.01`, trap density `n = 1`.
- Initial condition: `cm = ct = 0`.
- Mobile field `cm`: left end `= 100`, right end `= 0`. (`ct` has no boundary
  conditions; it only changes through the reaction.)
- Time step `dt = 0.05`, final time `t = 10`.

### Strong form

$$
\frac{\partial c_m}{\partial t} - \nabla \cdot (D \, \nabla c_m) = -R,
\qquad
\frac{\partial c_t}{\partial t} = R,
$$

with the reaction (trapping minus detrapping)

$$
R = k \, c_m \,(n - c_t) - p \, c_t .
$$

Mobile hydrogen lost to traps (`-R` in the first equation) reappears as trapped
hydrogen (`+R` in the second), so the coupling conserves material.

## What this example teaches (new compared to Example 3)

1. **Coupling fields through a reaction term.** The same expression `R` appears with
   opposite sign in the two equations, tying them together.
2. **Mixing element types.** `cm` uses a continuous (CG) element because it diffuses;
   `ct` uses a discontinuous (DG) element because it has no spatial derivatives. They
   are combined in one `mixed_element`.
3. **A field with no boundary condition** (`ct`), governed purely by the reaction.
4. **A nonlinear problem:** the `cm * ct` product makes the system genuinely nonlinear,
   which is exactly why we have been using the Newton (`NonlinearProblem`) interface
   all along.

## Files

- [`eg_4.ipynb`](eg_4.ipynb): the example explained step by step. **Start here.**
- [`my_script.py`](my_script.py): a scaffold (imports only) for **you** to fill in.
- [`eg_4.py`](eg_4.py): the complete reference script to check yourself against.

## Running

From the repository root, create and activate the environment (see the top-level
`environment.yml`), then:

```bash
conda activate dolfinx-tutorial-env
jupyter lab eg_4.ipynb   # read through the example
python my_script.py      # run your own version
python eg_4.py           # compare against the reference
```

Running the script writes `trapping_cm.bp` and `trapping_ct.bp`, which you can open in
[ParaView](https://www.paraview.org/).
