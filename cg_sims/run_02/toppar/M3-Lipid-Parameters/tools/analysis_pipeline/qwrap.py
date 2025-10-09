from MDAnalysis import transformations
import MDAnalysis as mda
import numpy as np
import sys

# Diffusion analysis - added to analysis_lipids by Balazs Fabian <fbalazsf@gmail.com>

""" Create CoM trajectory and unwrap it correctly.
    The molecules must be whole. Otherwise, MDA has
    to make it so, which is slow.
"""

def qwrap(trp, xtc, out_gro, out_xtc, lipid_names_str):
    u = mda.Universe(trp, xtc)

    # Add transformation to make whole (slow!)
    if False:
        workflow = [transformations.unwrap(u.atoms)]
        u.trajectory.add_transformations(*workflow)

    # TODO: this must be general enough!
    membrane = u.select_atoms(f"resname {lipid_names_str}")

    # dummy CoM group for writing trajectory
    dummyCoMGroup = sum([ lipid.atoms[0] for lipid in membrane.residues ])

    # write a gro file
    dummyCoMGroup.write(out_gro)

    # Perform the "0-th" step
    # xw: wrapped
    # xu: unwrapped
    xw_prev = np.array([ lipid.atoms.center_of_mass() for lipid in membrane.residues ])
    xu_prev = xw_prev

    with mda.Writer(out_xtc, dummyCoMGroup.n_atoms) as W:
        for ts in u.trajectory:

            L = u.dimensions[:3]
            xw = np.array([ lipid.atoms.center_of_mass() for lipid in membrane.residues ])

            delta = xw - xw_prev
            shift = np.floor(delta/L + 0.5)*L
            xu = xu_prev + delta - shift

            dummyCoMGroup.positions = xu
            W.write(dummyCoMGroup)

            # Update
            xw_prev = xw
            xu_prev = xu

def main():
    # Needs sims tpr, xtc, out com.gro name, out xtc name, and string list of lipid names
    qwrap(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

if __name__ == "__main__":
    main()
