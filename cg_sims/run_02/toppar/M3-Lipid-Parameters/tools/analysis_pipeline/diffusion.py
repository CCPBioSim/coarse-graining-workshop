import Dfit
import numpy as np
import MDAnalysis as mda
import sys

# Diffusion analysis - added to analysis_lipids by Balazs Fabian <fbalazsf@gmail.com>

def diffusion(conf, traj, out_file, glsIntervalLower, glsIntervalUpper, lipid, leaflet_group):
    # Load trajectory
    #conf = sys.argv[1]
    #traj = sys.argv[2]
    u = mda.Universe(conf,traj)

    # Remove Center of Mass of the first frame
    #pos = u.atoms.positions - u.atoms.center_of_mass()
    pos = u.atoms.positions - u.atoms.center_of_geometry() # needed as MDAnalysis does not capture a mass in all cases current COM input conf .gro all have the same mass so this shoudl be ok
    u.atoms.positions = pos

    # Determine time step
    dt = u.trajectory.dt

    # Choose leaflet.
    # NOTE: after Center of Mass removal,
    # one leaflet must be at z<0 and the other at z>0
    if str(leaflet_group) == "0":
        symbol = "<"
    elif str(leaflet_group) == "1":
        symbol = ">"
    else:
        print(f"ERROR leaflet_group={leaflet_group} not supported")

    # Create selection for analysis
    # leaflet: the complete upper/lower leaflet irrespective of lipid type
    leaflet = u.select_atoms(f'prop z {symbol} 0')
    # sel: a single lipid type in the above leaflet
    sel = u.select_atoms(f'resname {lipid} and prop z {symbol} 0')

    # Remove CoM and convert to nm
    #pos = [ 0.1 * (sel.positions[:,:2] - leaflet.atoms.center_of_mass()[:2]) for ts in u.trajectory ]
    pos = [ 0.1 * (sel.positions[:,:2] - leaflet.atoms.center_of_geometry()[:2]) for ts in u.trajectory ] # see same change above
    pos = np.array(pos)
    pos = np.swapaxes(pos,0,1)
    # Convert to per-particle list
    pos = [ mol for mol in pos ]

    # Analyse using Dfit
    #print(f"  Dfit.Dfit.Dcov(m=20,fz={pos},dt={dt},tmin={glsIntervalLower},tmax={glsIntervalUpper},fout={out_file},nitmax=200)")
    res = Dfit.Dfit.Dcov(m=20,fz=pos,dt=dt,tmin=int(glsIntervalLower),tmax=int(glsIntervalUpper),fout=out_file,nitmax=200)
    res.run_Dfit()
    res.analysis(tc=10*dt)

def main():
    diffusion(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

if __name__ == "__main__":
    main()
