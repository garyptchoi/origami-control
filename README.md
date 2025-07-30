# Rigidity control of general origami structures

<img src = "https://github.com/garyptchoi/origami-control/blob/main/cover.jpg" height="300" />

This repository contain simulation codes for controlling the rigidity of general origami structures by incrementally adding facet planarity constraints.

At every step, k random candidate facets are chosen and one among them is selected based on some selection rule.

Changing the number of choices k and the selection rule allows us to effectively control the explosive rigidity percolation behavior in general origami.

Any comments and suggestions are welcome. 

If you use this code in your work, please cite the following paper:

R. Li and G. P. T. Choi,
"[Rigidity control of general origami structures.](https://arxiv.org/abs/2507.16934)"
Preprint, arXiv:2507.16934.

Copyright (c) 2025, Rongxuan Li and Gary P. T. Choi

https://github.com/garyptchoi/origami-control

===============================================================
# Preparation codes:

## Folder: `Preprocess`->`origami_source`

### Instructions

1. Download the `.obj` file of the desired pattern from [https://origamisimulator.org/](https://origamisimulator.org/).
2. Change the file extension from `.obj` to `.m` save in the `origami_source` folder.
3. Use the following naming convention:
origami_source/`<OrigamiType>`/`<OrigamiType>`_`<Size>`_`<FoldingPercentage>`%.m

- `<OrigamiType>`: type of the origami structure (e.g., `auxetic_triangle`)
- `<Size>`: resolution index (e.g., `2`)
- `<FoldingPercentage>`: folding percentage (e.g., `25%`)

---

## File: `Preprocess`->`read_file.m`

### Function
This script reads and processes origami or kirigami geometry exported from [https://origamisimulator.org/](https://origamisimulator.org/). It visualizes and saves the vertices, edges, and facet information of the pattern.


### Instructions

Modify the input file list/path **around line 30** in `read_file.m` to match the structure you want to load,  
   for example:

   ```matlab
   % Example path
   file_list  = 'origami_source/auxetic_triangle/auxetic_triangle_2_25%.m'';
```
---
## Folder: `Preprocess`->`origami_data`

### Function

The output of `Preprocess`->`read_file.m` will be saved in `.mat` format inside the `origami_data/` folder. The following variables are generated:

- `V` — vertex coordinates matrix (size: `n × 3`), where each row is `[x_i, y_i, z_i]`
- `F` — facet matrix (triangles), size `m × 3`, each row is `[v1, v2, v3]`
- `E` — edge list, size `k × 2`, each row is `[v1, v2]`
- `L` — edge label vector, size `k × 1`

---

### Edge Labeling Rules

| Label | Type                  | Color      |
|-------|-----------------------|------------|
| 0     | Marginal edge         | Light gray |
| 1     | Planar diagonal edge  | Red        |
| 2     | Mountain edge         | Green      |
| 3     | Valley edge           | Blue       |

---
# Simulation codes: 

## File: `Simulation.m`


### Function
This is the **main simulation scrip** for computing the Degrees of Freedom (DOF) during the constraint sampling process.

For each simulation, the script performs a complete DOF evolution by progressively adding constraints to the origami structure.


### Output
- The simulation results are saved as a matrix in the `simulation_results/` folder.
- The output matrix has size: [`n_sim` × (`num_facets` + 1)]

where:
- `n_sim` is the number of simulations 
- `num_facets` is the number of facets in the structure

Each row corresponds to one simulation and contains the DOF values recorded after each step of constraint addition.

---
# Simulation results:
- The simulation results are saved as `.mat` files in the `simulation_results/` folder.

- Each file is named using the following convention: `simulation_results`/`<OrigamiType>`/`<OrigamiType>`_`<Size>`_`<FoldingPercentage>`%.m

---
# Analysis codes:

- `analysis_dof_norm.m`  
  Main script for analyzing **normalized Degrees of Freedom (DOF)** versus constraint density **$\rho$** for general origami structures.
  
  **plotted results** will be saved in folder `plots_k_norm`, data loaded from `simulation_results`. 

- `analysis_P_vs_r.m`    
  Main script for analyzing the **probability of achieving a minimum DOF origami (P)** at each constraint density $\rho$, computed from `n_sim` independent simulations.

  **plotted results** will be saved in folder `plots_P_k`, data loaded from `simulation_results`. 

- `analysis_P_summary.m`  
  Computes the **triangular facet ratio and probability (P)** along with other parameters in folder.

  **numerical results**  will be saved in folder `P_summary_data`, data loaded from `simulation_results`. 

- `analysis_3d_fitting.m`  
  Main scripts for analyzing the **critical transition** using both structural and selection parameters across different types of general origami structures. Data loaded from `P_summary_data` folder.