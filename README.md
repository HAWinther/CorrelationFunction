# CorrelationFunction
Author: Hans A. Winther, 2018

Very simple code to compute the radial monopole correlation function for data from a galaxy survey / mocks.
Written for test purposes, but has the same (or faster) speed as [CUTE](https://github.com/damonge/CUTE) with OpenMP

Input: galaxyfile, randomfile, nbins, rmax and OmegaM
Output: [xi(r), err(r), DD(r), DR(r), RR(r)]

Input is assumed to have format: [RA, DEC, z, (weight)] with RA/DEC in degrees. Rmax in Mpc/h

Defines:

- WEIGHTS : Use weighted galaxies

- OMP : Parallelize using OpenMP

- BRUTEFORCE : Brute-force pair-counts instead of grid (linked list)

- OUTPUT\_PAIRS\_TO\_SCREEN : print paircount to screen as we go along


