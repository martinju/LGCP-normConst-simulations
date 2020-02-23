
# <img src="logo.png" align="right" height="50px"/>

# LGCP normalization constant approximation - simulation source code

This repository contains the results and source code for the simulation
experiments in the paper

*Investigating mesh based approximation methods for the normalization
constant in the Cox process likelihood* by Martin Jullum.

A link to the pre-print version of the paper will be made available here
in due time.

#### Main files

  - `source_code/help_functions.R` Functions containing the
    implementation of all the differnet approximation methods for the
    integral in question, in addition to a series of help functions
    related to the simulation experiments.
  - `source_code/full_simulations_script.R` The script used to carry out
    the full simulation experiment in the paper. The script also
    contains the short test version that quickly runs trough some
    parameter combinations (to check that the script works as intended).
  - `source_code/permutation_tests.R` The script used to carry out the
    permutation tests mentioned in the paper.
  - `source_code/illustrating_mesh_approximated_field.R` A
    self-contained script for illustating how to approximate a given
    field on specific mesh using the finite element method (FEM).

#### Full result table for simulations experiments

A sort and searchable table with the full simulation results is
available
[here](https://martinju.github.io/LGCP-normConst-simulations/sim_res.html).

#### Permutation tests

A sort and searchable table with results for permutation tests of
pairwise differences between all combinations of the approximation
methods, for each parameter combination is available
[here](https://martinju.github.io/LGCP-normConst-simulations/permut_tests.html).

Any questions about the paper, simulations or scripts can be raised by
opening a GitHub issue, or sending an email to the author
[jullum@nr.no](mailto:jullum@nr.no?subject=LGCP-normConst-simulations).

<!-- ### Source code files -->

<!-- All source code for the simulations are available under /Source code -->

<!-- File  | Description -->

<!-- ------------- | ------------- -->

<!--   simulation_script.R                         | This is the main script executing the simulation experiment in the paper. Settings for full simulation is commented out, such that the script can be ran quickly.  -->

<!-- timing_deterministic_integration_methods.R  | Script for timing the deterministic integration methods. -->

<!-- help_function.R | All functions used in the above scripts. Ufortunatly, almost no documentation is available.  -->
