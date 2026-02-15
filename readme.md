# conCluster_search

`conCluster_search` finds the **strongest contiguous voxel cluster of fixed size** from an arbitrary set of 3D voxel coordinates and values. It is designed for neuroimaging and volumetric analyses where voxel locations are **sparse**, **irregular**, and **not defined on a full 3D grid**.

Unlike naive approaches that select the top-K voxels by value (which can break connectivity), this function **explicitly enforces spatial connectivity** throughout the search.

---

## Edit history
- 02/15/2026 by Alex He - created first working version

---

## Key Features

- ✅ Works on **arbitrary voxel lists** (no full image required)
- ✅ Enforces **3D connectivity** (6 / 18 / 26)
- ✅ Finds a **fixed-size connected subcluster**
- ✅ Supports customizable scoring (e.g. `mean`, `sum`)
- ✅ Two visualization modes:
  - `volshow` (volume rendering)
  - explicit voxel cubes (`patch`)
- ✅ Robust to sparse and irregular voxel distributions

---

## Function Signature

```matlab
[maxClusterValue, maxClusterIndices, maxClusterCoords] = ...
    conCluster_search(input_coords, input_voxels, ...
                      cluster_size, connectivity, ...
                      average_method, plot_mode)
```

---

## Inputs

| Argument | Description |
|--------|-------------|
| `input_coords` | **[N × 3]** integer voxel coordinates `(X, Y, Z)` |
| `input_voxels` | **[N × 1]** voxel values aligned with `input_coords` |
| `cluster_size` | Desired number of voxels in the cluster (default: `5`) |
| `connectivity` | 3D connectivity: `6`, `18`, or `26` (default: `6`) |
| `average_method` | Function handle for scoring clusters (e.g. `@mean`, `@sum`) |
| `plot_mode` | Visualization mode: `'3d'`, `'patch'`, or `false` |

---

## Outputs

### `maxClusterValue`

Scalar value representing the **score of the selected cluster**, computed as:

```matlab
average_method(values_of_selected_voxels)
```

Examples:
- If `average_method = @mean`: mean voxel value
- If `average_method = @sum`: total cluster strength

---

### `maxClusterIndices`

**[cluster_size × 1]** indices into the original input arrays  
(`input_coords`, `input_voxels`) identifying which voxels belong to the cluster.

```matlab
input_coords(maxClusterIndices, :) == maxClusterCoords
```

---

### `maxClusterCoords`

**[cluster_size × 3]** voxel coordinates of the selected cluster, in the **same coordinate system as `input_coords`**, ordered consistently with `maxClusterIndices`.

---

## Algorithm Overview

1. **Local bounding box construction**  
   Converts sparse voxel coordinates into a compact local volume for efficient connectivity analysis.

2. **Threshold sweep**  
   Progressively lowers the voxel-value threshold until at least one connected component contains ≥ `cluster_size` voxels.

3. **Connectivity-preserving shrinking**  
   For each candidate component larger than `cluster_size`, a greedy growth algorithm selects a **connected subcluster** of exactly `cluster_size` voxels that maximizes the chosen score.

4. **Cluster selection**  
   The highest-scoring candidate cluster is returned.

---

## Visualization Modes

### `plot_mode = '3d'`

Uses `viewer3d` + `volshow`:

- Grayscale semi-transparent background volume
- Red opaque cluster overlay
- Volume axes are permuted so `(X, Y, Z)` matches Cartesian coordinates

Best for:
- Large voxel sets
- Quick inspection
- Interactive slicing

---

### `plot_mode = 'patch'`

Uses explicit voxel cubes (`patch`):

- Background voxels: grayscale, highly transparent
- Cluster voxels: solid red cubes
- True Cartesian axes with lighting

Best for:
- Precise spatial interpretation
- Figures and presentations
- Small to moderate voxel counts

---

## Example Usage

```matlab
[maxVal, idx, coords] = conCluster_search( ...
    voxel_coords, voxel_values, ...
    10, 18, @mean, 'patch');
```

---

## Notes & Limitations

- `cluster_size` must be ≤ number of input voxels.
- Coordinates must be **integer voxel indices**.
- The cluster shrinking step is greedy but connectivity-safe; it is not guaranteed to find the global optimum in adversarial cases, but works well in practice.
- Visualization assumes `(X, Y, Z)` Cartesian interpretation.

---

## Typical Use Cases

- Identifying peak clusters in statistical maps
- EEG / MEG source-space voxel selection
- Sparse fMRI / PET voxel analyses
- Any 3D point cloud where spatial contiguity matters

---

## License

This project is licensed under the BSD 3-Clause License.  
See the [LICENSE](LICENSE) file for details.
