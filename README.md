# CorrelationFunction
Author: Hans A. Winther, 2018

A very simple library for use in C / C++ to compute the radial monopole correlation function for data from a galaxy survey / mocks using the LZ estimator or from a periodic simulations box. Using grids to speed up the computation. Written for test purposes and well commented so should be very easy to modify. Has the same (or faster) speed as [CUTE](https://github.com/damonge/CUTE) with OpenMP (but does not have half of the features as CUTE so use that or SWOT instead for scientific purposes).

See [test\_simulation] for how to use call the library for periodic simulations boxes.

See [test\_survey] for how to use call the library for survey data / galaxy mocks (NB: no masks implemented).

Input files in this survey case is assumed to either have format: [RA, DEC, z, (weight)] with RA/DEC in degrees or [X, Y, Z, (weight)] with positions in Mpc/h.
Rmax is also given in Mpc/h. One can also call the code directly by supplying an array of galaxies (see cuter\_library.h for how different call options) defined in galaxy.h This struct can be freely modified as long as it has x[3] and w (if weights is used).

The library returns a struct with [r, xi(r), poisson\_err(r), DD(r), DR(r), RR(r)] or [r, xi(r), poisson\_err(r), DD(r)] for a periodic box. 

Requires the GSL library (only needed to compute the radial co-moving distance-redshift relation r(z))

Defines:

- WEIGHTS : Use weighted galaxies

- USE\_OMP : Parallelize using OpenMP

- USE\_MPI : Parallelize using MPI (very simple: all tasks store all data)

- BRUTEFORCE : Brute-force paircounts instead of grid

Compile by running [make OPTIONS] where OPTIONS are: USE\_MPI=1 for MPI, USE\_OMP=1 for OpenMP, USE\_BRUTEFORCE=1 for brute-force and USE\_WEIGHTS=1 for including weights.

