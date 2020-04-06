#---------------------------------------------------------------------------------------
#This code takes a long MD trajectory, aligns the heavy atoms to the average structure,
#splits the trajectory into segments of different length 
#for Quasi-harmonic entropy calculation
#---------------------------------------------------------------------------------------

import numpy as np

import MDAnalysis as md

from MDAnalysis import transformations

import MDAnalysis.analysis.encore.covariance as cov

u = md.Universe('../ionized.psf','../ATATAT_production_NPT.dcd')

base = u.select_atoms('bynum 1:382 and not name H* and not (resid 1 or resid 6 or resid 7 or resid 12)')

#write the pdb file
with md.Writer("TATA.pdb",base.n_atoms) as W:
    if u.trajectory.frame == 0:
        W.write(base)


#get average structure
numframes = len(u.trajectory)
p_avg = np.zeros_like(base.positions)
for ts in u.trajectory:
    p_avg += base.positions
p_avg /= numframes

base.positions = p_avg

with md.Writer("TATA_avg.pdb",base.n_atoms) as W:
    W.write(base)
u.trajectory.rewind()

#align by removing rotation and translation
ref = md.Universe("TATA_avg.pdb")
transform = transformations.fit_rot_trans(base, ref, weights="mass")
u.trajectory.add_transformations(transform)

#write trajectories
#discard first 5 ns (each frame = 10 ps)

eqbm_time = 500 #1 ns

t = np.arange(eqbm_time,4500,100)

for time in t:
    with md.Writer("TATA_traj_%d_ps.dcd"%(int(time*10)),base.n_atoms) as W:
        for ts in u.trajectory:
            if ts.frame > eqbm_time and ts.frame < eqbm_time + time:
                W.write(base)

