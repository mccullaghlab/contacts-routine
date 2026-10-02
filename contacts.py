import itertools

import click
import numpy as np
import MDAnalysis
from MDAnalysis.lib.distances import calc_bonds
from tqdm import tqdm

from trajectory_files import trajectory_files


@click.command(
    no_args_is_help='-h',
    help='Mindist computes the minimal distance between pairs of residues.',
)
@click.option(
    '--trajectory',
    '-f',
    'trajfile',
    type=click.Path(exists=True),
    help='Path to any trajectory file supported by MDAnalysis',
)
@click.option(
    '--trajectory-list',
    '--traj-list',
    'trajlist',
    type=click.Path(exists=True),
    help='Path to a trajectory folder OR file containing trajectory paths',
)
@click.option(
    '--topology',
    '-s',
    'topfile',
    required=True,
    type=click.Path(exists=True),
    help='Path to topology file (.tpr or .pdb)',
)
@click.option(
    '--index',
    '-n',
    'ndxfile',
    required=True,
    type=click.Path(exists=True),
    help='Path to index file, 1 (!) indexed of shape (n, 2).',
)
@click.option(
    '--output',
    '-o',
    required=True,
    type=click.Path(),
    help='Path to output file',
)
def main(trajfile, trajlist, topfile, ndxfile, output):
    traj_files = trajectory_files(trajfile, trajlist)
    
    print(f"Processing {len(traj_files)} trajectory file(s)")
    if len(traj_files) <= 10:
        for f in traj_files:
            print(f"  - {f}")
    else:
        print(f"  - {traj_files[0]}")
        print(f"  - ...")
        print(f"  - {traj_files[-1]}")
    
    # load files
    index_pairs = np.loadtxt(ndxfile, dtype=int)
    # Validate residue numbering
    universe = MDAnalysis.Universe(topfile, traj_files[0])
    pdb_resids = sorted(set(universe.residues.resids))

    # Check for negative residues
    if any(r < 0 for r in pdb_resids):
        raise click.UsageError(
            f"PDB contains negative residue indices. Please renumber starting from 1."
        )

    # convert residue indices to heavy atoms
    atoms_per_res = {
        index: universe.select_atoms(
            f'resid {index} and not (type H or name H*)'
        ).indices.tolist()
        for index in np.unique(index_pairs)
    }

    missing_resids = [resid for resid, atoms in atoms_per_res.items() if not atoms]
    if missing_resids:
        raise click.UsageError(
            f"Index file contains residue IDs with no heavy atoms: {missing_resids}"
        )

    atom_pairs = [
        list(
            itertools.product(
                atoms_per_res[i], atoms_per_res[j],
            ),
        ) for i, j in index_pairs
    ]

    index_output = np.concatenate([
    [
        (res_i, res_j, atom_i + 1, atom_j + 1)  # only atoms get +1
        for atom_i, atom_j in atom_pair
    ]
    for (res_i, res_j), atom_pair in zip(index_pairs, atom_pairs)
    ])
    np.savetxt(
        f'{output}.atom_indices',
        index_output,  # no global +1
        fmt='%.0f',
        header='res_i res_j atom_i atom_j',
    )

    # flatten atom_pairs to loop over all
    atom_pairs = np.concatenate(atom_pairs)

    with open(output, 'w') as ostream:
        for traj_idx, traj_file in enumerate(traj_files):
            print(f"Processing trajectory {traj_idx + 1}/{len(traj_files)}: {traj_file}")
            for distances in compute_distances(
                traj_file, topfile, atom_pairs,
            ):
                ostream.write(
                    ' '.join([f'{d:.5f}' for d in distances]) + '\n',
                )


def compute_distances(trajfile, topfile, atom_pairs):
    universe = MDAnalysis.Universe(topfile, trajfile)
    first_atoms = universe.atoms[atom_pairs[:, 0]]
    second_atoms = universe.atoms[atom_pairs[:, 1]]
    for _ in tqdm(universe.trajectory, desc=f"Processing {trajfile}"):
        # MDAnalysis reports Angstrom; the pipeline's distance files use nm.
        yield calc_bonds(
            first_atoms.positions,
            second_atoms.positions,
            box=universe.dimensions,
        ) / 10.0


if __name__ == '__main__':
    main()
