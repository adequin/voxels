# Voxels: A minimal 3D voxel-beam lattice toolkit (MATLAB)

> Lightweight utilities to assemble voxel lattices, generate 3D beam element configurations, build a global stiffness matrix, solve with Dirichlet boundary conditions, and (optionally) infer top‑face displacements from a few strain sensors. This README was AI generated.

---

## Table of contents

* [Overview](#overview)
* [Coordinate system & conventions](#coordinate-system--conventions)
* [Quick start](#quick-start)
* [Repository layout](#repository-layout)
* [Core API](#core-api)

  * [`get_lattice`](#get_lattice)
  * [`define_nodes`](#define_nodes)
  * [`config_array`](#config_array)
  * [`compute_matrix`](#compute_matrix)
  * [`solve_with_dirichlet`](#solve_with_dirichlet)
  * [`estimate_with_sensors` & `estimate_topface`](#estimate_with_sensors--estimate_topface)
  * [`get_xyz`](#get_xyz)
  * [Plotting helpers (`beamCoords`, `plotBeamDeform`)](#plotting-helpers-beamcoords-plotbeamdeform)
* [Examples](#examples)

  * [Compression/Tension via sensors](#compressiontension-via-sensors)
  * [Three-voxel bending test](#three-voxel-bending-test)
* [Assumptions & limitations](#assumptions--limitations)
* [Troubleshooting](#troubleshooting)
* [Roadmap / TODO](#roadmap--todo)
* [License](#license)

---

## Overview

This library builds a voxelized beam lattice, produces element-level orientation/section data, assembles a global stiffness matrix for 3D Euler–Bernoulli/Timoshenko-style beams, and solves for nodal displacements using displacement (Dirichlet) boundary conditions. It also supports a simple sensor model that infers a rigid top‑face motion ((u_z,;\theta_x,;\theta_y)) from a sparse set of axial strain readings on vertical members.

**What you get**

* Voxel indexing and face/edge node bookkeeping.
* Connector edge discovery between neighboring voxels.
* Per‑element local bases and cross‑section properties.
* Sparse global stiffness assembly (12×12 per beam) in one pass.
* Displacement‑BC solver and reaction extraction.
* Small plotting utilities for undeformed/deformed shape.
* Sensor‑based top face fit for quick experiments.

Minimum dependency is MATLAB (no toolboxes assumed).

---

## Coordinate system & conventions

* **Global axes**: X right, Y "into the page", Z up.

* **Face numbering** and local UV per face (from `define_nodes.m`):

  | Face | Normal | Up | Right |
  | ---- | ------ | -- | ----- |
  | f1   | −Y     | +Z | +X    |
  | f2   | +X     | +Z | +Y    |
  | f3   | +Y     | +Z | −X    |
  | f4   | −X     | +Z | −Y    |
  | f5   | +Z     | +Y | +X    |
  | f6   | −Z     | −Y | −X    |

* **Node ordering**: each voxel contributes 42 node IDs:

  * 12 edges × 3 points each = 36, plus 6 face centers.
  * Node IDs are contiguous per voxel; face helpers return named nodes like `f5.tm`, `f5.c`, etc.

* **DOF ordering (per node)**: `[ux, uy, uz, rx, ry, rz]` → 6 DOF/node.

* **Element orientation**: local basis `[ex; ey; ez]` is built from element axis and the host face normal/up so the local stiffness can be rotated to global via `K_beam = R * K_local * R'` where `R = blkdiag(R3,R3,R3,R3)` and `R3 = [ex'; ey'; ez']`.

---

## Quick start

```matlab
% 1) Define a 3-voxel vertical stack (nx1x1 logical array)
A = [1;1;1];
lat = get_lattice(A);          % lattice + defaults
cfg = config_array(lat);       % element list & properties
K   = compute_matrix(cfg);     % global sparse stiffness

% 2) Build coordinates at each node id (global)
coords = beamCoords(lat, cfg);

% 3) Simple Dirichlet solve: fix bottom face, push top center down
n_top = define_nodes(lat,1,1,1);      % top voxel nodes
n_bot = define_nodes(lat,3,1,1);      % bottom voxel nodes

bc_idx = [];
bc_val = [];
% fix all DOF of bottom face nodes
bot = [];
for nm = n_bot.face_nodes, bot = [bot, n_bot.f6.(nm)]; end
for node = bot, bc_idx = [bc_idx, 6*node-5:6*node]; bc_val = [bc_val, zeros(1,6)]; end
% impose uz at top center
top_c = n_top.f5.c; bc_idx = [bc_idx, 6*top_c-3]; bc_val = [bc_val, -1e-3];

U = solve_with_dirichlet(K, bc_idx, bc_val);
plotBeamDeform(cfg, coords, U, 50);
```

---

## Repository layout

```
compute_matrix.m          % element rotation + sparse global assembly
config_array.m            % voxel + connector element generation
get_lattice.m             % lattice creation + defaults
solve_with_dirichlet.m    % partitioned solve w/ reactions
get_xyz.m                 % node id → (x,y,z), with optional lattice offset

% Sensor-based estimation
estimate_with_sensors.m   % wrapper + top-face fit (dz,thx,thy)

% Node bookkeeping
define_nodes.m            % face/edge maps, normals, up/right, ids

% Example / utility code (grouped in main.m)
main.m                    % example flows + plotting helpers
```

---

## Core API

### `get_lattice`

```matlab
lat = get_lattice(A, conn_l, conn_b, conn_h, conn_e, conn_g, voxel_size, chamfer_size)
```

* **A**: logical 3D array (Rows×Cols×Slices) indicating voxel presence.
* Returns a struct with `A`, `dims`, `id` (voxel ids), connector/voxel dimensions and material defaults.

### `define_nodes`

```matlab
n = define_nodes(lat, r, c, s)  % or define_nodes(vox_id)
```

* Returns a struct with face fields `f1..f6`, each exposing named nodes (e.g. `bm`, `c`, `tm`…), plus `normal/up/right` per face.

### `config_array`

```matlab
array_config = config_array(lat)
```

* Builds the **element table** by concatenating voxel internals and inter‑voxel connectors.
* Each row packs: `[node1, node2, ex, ey, ez, L, E, G, b, h]` (orientation as three row vectors, section and material).

### `compute_matrix`

```matlab
K = compute_matrix(array_config)
```

* Assembles a sparse global stiffness (size `6N × 6N`).
* Uses local 12×12 beam stiffness then rotates to global: `K_beam = R * K_local * R'`.
* Symmetrizes for numerical safety.

### `solve_with_dirichlet`

```matlab
[U, reactions] = solve_with_dirichlet(K, bc_idx, bc_val, f)
```

* **bc_idx**: global DOF indices to constrain; **bc_val**: prescribed values.
* **f**: optional global force vector (`6N×1`).
* Returns **U** as `N×6` (reshape is handled internally) and **reactions** at constrained DOFs.

### `estimate_with_sensors` & `estimate_topface`

```matlab
U = estimate_with_sensors(lat, K, beam_n1, beam_n2, readings)
[bc_idx, bc_val, fit] = estimate_topface(lat, beam_n1, beam_n2, readings, weights)
```

* Fits (u_z(x,y) = dz + \theta_x y - \theta_y x) over the **top face** from axial strains on several vertical members.
* Assumes columnar assembly (size `n×1×1` without gaps) to compute the global height used in the design matrix.
* Produces top‑face `uz` BCs and clamps the bottom face.

### `get_xyz`

```matlab
xyz = get_xyz(node_id)           % local voxel coordinates
xyz = get_xyz(node_id, lat)      % global coordinates (with lattice spacing)
```

### Plotting helpers (`beamCoords`, `plotBeamDeform`)

* `beamCoords(lat, cfg)` collects coordinates for all node ids appearing in `cfg`.
* `plotBeamDeform(cfg, coords, U, scale, sensor_beams)` overlays undeformed and deformed beams; optionally highlights sensor segments.

---

## Examples

### Compression/Tension via sensors

```matlab
A = [1;1;1];
lat = get_lattice(A);
cfg = config_array(lat);
K   = compute_matrix(cfg);

n = define_nodes(1);
sensor_beams = [n.f1.c  n.f1.bm;
                n.f2.c  n.f2.bm;
                n.f3.c  n.f3.bm];
readings = -[250e-6, 250e-6, 250e-6];   % axial strain
U = estimate_with_sensors(lat, K, sensor_beams(:,1), sensor_beams(:,2), readings);
coords = beamCoords(lat,cfg);
plotBeamDeform(cfg, coords, U, 50, sensor_beams);
```

### Three-voxel bending test

```matlab
A = [1,1,1];
lat = get_lattice(A);
cfg = config_array(lat);
K = compute_matrix(cfg);

n1 = define_nodes(1);
n2 = define_nodes(2);
n3 = define_nodes(3);

bc_idx = []; bc_val = [];
% support: pick nodes on bottom faces of voxels 1 & 3
bot_nodes = [n1.f6.tr, n1.f6.br, n3.f6.tl, n3.f6.bl];
for node = bot_nodes
    bc_idx = [bc_idx, 6*node-5:6*node-3];
    bc_val = [bc_val, 0,0,0];
end
% impose uz at top nodes of middle voxel
z_displ = -1e-2; % 1 cm
for node = [n2.f5.tm, n2.f5.c, n2.f5.bm]
    bc_idx = [bc_idx, 6*node-3];
    bc_val = [bc_val, z_displ];
end

U = solve_with_dirichlet(K, bc_idx, bc_val);
coords = beamCoords(lat,cfg);
plotBeamDeform(cfg, coords, U, 1);
```

---

## Assumptions & limitations

* Cross‑sections are rectangular (width **b**, height **h**), with Saint‑Venant torsion (J) approximation.
* Some connector/“stub” dimensions are approximations; parallel‑axis corrections for compound sections are TODO.
* Sensor fit assumes uniform axial strain through height (i.e., linear top‑face motion over (x,y)).
* Lattice in `estimate_topface` must be a continuous vertical column (`size(A,2)==1 & size(A,3)==1`).

---

## Troubleshooting

* **"Matrix is close to singular"**: over‑constrained or under‑constrained model; check `bc_idx` duplicates and ensure adequate supports.
* **Dimension errors in plotting**: ensure `U` is `N×6` (function reshapes if needed) and `coords` covers all node ids used by `cfg`.
* **Wrong element orientation**: verify `local_basis` inputs; element axis must be parallel to face directions or normal within tolerance.
* **Duplicate BC indices**: `estimate_with_sensors` asserts uniqueness; if you add BCs, keep them disjoint.

---

## Roadmap / TODO

* Parallel‑axis corrections for stub beams; revisit (J) for non‑square aspect ratios.
* Cleaner boundary condition helpers (e.g. `bc = bc.addNodes(face, dofs, values)`).
* Pack element rows in a typed struct to reduce magic column indices.
* Unit tests for `local_basis` and `define_nodes` maps.
* Diagram(s) for face numbering and node aliases.

---

## License

Specify a license (e.g., MIT) in `LICENSE` if you intend to share or open‑source.
