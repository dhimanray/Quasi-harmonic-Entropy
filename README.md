# Quasiharmonic Entropy from Molecular Dynamics simulation

This repository contains codes to calculate conformational entropy of biomolecules using quasiharmonic approximation. The heavy atoms from each frame of MD trajectories are first aligned with respect to the average structure. Then they are split into trajectory segments of increasing length in order to check convergence of entropy. This part is accomplished by `trajectories.py` code. 
