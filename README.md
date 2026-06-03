# Contacts Routine

A pipeline for computing minimal residue-residue distances (contacts) from MD trajectories. Supports single and multiple trajectories.

---

## Overview

The pipeline runs in four sequential steps:

1. **`estimate_contacts.py`** — Computes the fraction of frames in which each residue pair is in contact (minimum heavy-atom distance ≤ 4.5 Å).
2. **`extract_indices.py`** — Selects residue pairs whose contact frequency falls within a user-defined window (min–max threshold).
3. **`contacts.py`** — Computes all-atom pairwise distances for the selected residue pairs across the trajectory.
4. **`extract_contacts.py`** — Identifies the minimal distances for atom pairs forming contacts above the minimum threshold, and writes the final output.

The entry point for the full pipeline is `run_contacts_routine.sh`.

---

## Requirements

- Python with: `MDAnalysis`, `mdtraj`, `msmhelper`, `numpy`, `click`, `tqdm`, `prettypyplot`, `matplotlib`
- Trajectory files in `.xtc` format
- A topology file (`.pdb` or `.tpr`)
- An index file (`.ndx`) listing residue pairs to analyze (1-indexed, shape `(n, 2)`, where n is the number of residue pairs)
- PDB residue numbering must be positive and sequential with no gaps (required by `mdtraj`)

---

## Usage

```bash
./run_contacts_routine.sh PARAM1 PARAM2 PARAM3 PARAM4 PARAM5 PARAM6 PARAM7 PARAM8
```

| Parameter | Description |
|-----------|-------------|
| `PARAM1` | `.xtc` trajectory file **or** path to folder containing `.xtc` files |
| `PARAM2` | Topology file (`.pdb`) |
| `PARAM3` | Minimum contact frequency threshold (0–1), e.g. `0.1` |
| `PARAM4` | Maximum contact frequency threshold (0–1), e.g. `0.9` |
| `PARAM5` | Index file (`.ndx`) of residue pairs |
| `PARAM6` | Base name for the system (used for output file names) |
| `PARAM7` | Trajectory mode: `single` or `multi` |
| `PARAM8` | Threshold mode: `overall` or `per-trajectory` |

### Threshold window

Contacts are selected if their formation frequency falls between `MIN_THR` and `MAX_THR` (inclusive). This allows filtering out both transient contacts (too rare) and permanently formed contacts (too stable). To apply only a lower bound, set `MAX_THR` to `1.0`.

### Threshold modes

- **`overall`**: Contact frequency is averaged across all frames from all trajectories. A contact is selected if the total average falls within the threshold window.
- **`per-trajectory`**: Contact frequency is computed per trajectory. A contact is selected only if it meets the threshold in **every** trajectory individually.

> Note: `per-trajectory` mode requires `multi` trajectory mode.

---

## Examples

**Single trajectory:**
```bash
./run_contacts_routine.sh traj.xtc system.pdb 0.1 0.9 indices.ndx my_system single overall
```

**Multiple trajectories — overall threshold:**
```bash
./run_contacts_routine.sh /path/to/traj_folder system.pdb 0.1 0.9 indices.ndx my_system multi overall
```

**Multiple trajectories — per-trajectory threshold:**
```bash
./run_contacts_routine.sh /path/to/traj_folder system.pdb 0.1 0.9 indices.ndx my_system multi per-trajectory
```

**For our HP35 benchmark example:**
```bash
./run_contacts_routine.sh traj.xtc system.pdb 0.3 1.0 indices.ndx HP35 single overall
```

When using `multi` mode, the folder should contain `.xtc` files named e.g.:
```
/path/to/traj_folder/traj1.xtc
/path/to/traj_folder/traj2.xtc
/path/to/traj_folder/traj3.xtc
```

---

## Output Files

| File | Description |
|------|-------------|
| `<s>.is_mindist` | Per-residue-pair contact fractions (from step 1) |
| `<s>.is_mindist.thr<MIN>-<MAX>.ndx` | Residue pairs within the threshold window (from step 2) |
| `<s>.all_thr<MIN>-<MAX>_selected_atom_distances` | All-atom pairwise distances for selected pairs (from step 3) |
| `<s>.all_thr<MIN>-<MAX>_selected_atom_distances.atom_indices` | Atom index mapping `(res_i, res_j, atom_i, atom_j)` (from step 3) |
| `<s>.mindist` | **Final output**: minimal distances per frame for selected contacts |
| `<s>.mindist.ndx` | Residue pair indices of the final selected contacts |

The pipeline skips step 1 if `<s>.is_mindist` already exists. The final output files (.mindist and .mindist.ndx) are those used to compute the similarity matrix and correlation-based clusters with MoSAIC. The intermediate output files are usually not used in our workflow and can be deleted (unless the user needs them for some different analysis).

---

## Pipeline Details

### Step 1 — `estimate_contacts.py`

Computes whether each residue pair is in contact in each frame using `MDAnalysis`. Parallelized across CPU cores using Python's `multiprocessing`. A contact is defined as a minimum heavy-atom distance ≤ **4.5 Å**. Outputs the contact fraction for each pair.

Before processing, validates the PDB residue numbering: raises an error if any residue IDs are negative or if the index file references residues not present in the PDB, and warns if there are gaps in residue numbering.

### Step 2 — `extract_indices.py`

Reads the contact fractions from step 1 and writes an index file containing only the residue pairs whose contact frequency falls within `[MIN_THR, MAX_THR]`.

### Step 3 — `contacts.py`

Uses `mdtraj` to compute distances between all heavy-atom pairs within each selected residue pair. Outputs a large distance matrix (one row per frame) and an atom index file.

Validates that PDB residue numbering is positive and fully sequential (no gaps), as `mdtraj` indexes residues by position. If your PDB has non-sequential numbering, renumber it first, e.g.:
```bash
gmx editconf -f in.pdb -o out.pdb -resnr 1
```

### Step 4 — `extract_contacts.py`

For each residue pair, identifies atom pairs that individually form contacts (distance ≤ **0.45 nm**) above the minimum threshold. Writes the per-frame minimal distance for each such residue pair to the final output file.

---

## Notes

- Index files are **1-indexed** throughout.
- Hydrogen atoms are excluded from distance calculations by default.
- The cutoff used in step 1 (`estimate_contacts.py`) is **4.5 Å**; the cutoff in step 4 (`extract_contacts.py`) is **0.45 nm** — these are equivalent.
