# Quasiharmonic Entropy from Molecular Dynamics simulation

This repository contains codes to calculate conformational entropy of biomolecules using quasiharmonic approximation. The heavy atoms from each frame of MD trajectories are first aligned with respect to the average structure. Then they are split into trajectory segments of increasing length in order to check convergence of entropy. This part is accomplished by `trajectories.py` code. Then the `entropy.py` code is run to calculate quasiharmonic entropy according to either Andricioaei-Karplus scheme or Schlitter scheme. 

The theoretical underpinning of these methods are descride in the following publication:
Andricioaei, Ioan, and Martin Karplus. "On the calculation of entropy from covariance matrices of the atomic fluctuations." The Journal of Chemical Physics 115.14 (2001): 6289-6292. 
https://aip.scitation.org/doi/abs/10.1063/1.1401821
