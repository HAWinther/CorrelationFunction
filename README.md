# CorrelationFunction
Author: Hans A. Winther, 2018

Very simple code to compute the radial monopole correlation function for data from a galaxy survey / mocks. Using grids to speed up the computation. Written for test purposes and well commented so should be very easy to modify. Has the same (or faster) speed as [CUTE](https://github.com/damonge/CUTE) with OpenMP (but does not have half of the features as CUTE so use that or SWOT instead for scientific purposes).

Input: galaxyfile, randomfile, nbins, rmax and OmegaM
Output format: [r, xi(r), err(r), DD(r), DR(r), RR(r)] 
where err are Poisson errors.

Input is assumed to have format: [RA, DEC, z, (weight)] with RA/DEC in degrees. Rmax in Mpc/h

Requires the GSL library (only needed to compute the radial co-moving distance-redshift relation r(z))

All the code is for simplicity included in the file [main.c], but see [factored\_source] for a properly factored code.

Defines:

- WEIGHTS : Use weighted galaxies

- OMP : Parallelize using OpenMP

- BRUTEFORCE : Brute-force pair-counts instead of grid

- OUTPUT\_PAIRS\_TO\_SCREEN : print paircount to screen as we go along

