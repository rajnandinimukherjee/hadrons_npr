# NPR production code

Computes propagators and vertex functions for non-perturbative renormalization in regularization independent momentum-subtraction schemes with MOM and SMOM kinematics based on the [NPR modules](https://github.com/aportelli/Hadrons/tree/develop/Hadrons/Modules/MNPR) of [Hadrons](https://github.com/aportelli/Hadrons).
The program can optionally compute fourquark vertex functions or `QED_L` corrections to first order in the electromagnetic coupling constant.

## Execution
Install [Grid](https://github.com/paboyle/Grid)
and [Hadrons](https://github.com/aportelli/Hadrons) on their `develop` branch.
For CPU architectures compile via:
``` bash
./bootstrap.sh
mkdir build
cd build
../configure --with-hadrons=<Hadrons install prefix> (CPPFLAGS=-DMOBIUS) (CPPFLAGS=-DBICGSTAB)
make
```
For Nvidia GPU accelerated architectures cross-compile via:
``` bash
./bootstrap.sh
mkdir build
cd build
../configure --host=aarch64-linux-gnu --with-hadrons=<Hadrons install prefix> CC=gcc CXX=nvcc CPPFLAGS="-x cu -DBICGSTAB"
make
```
by default the `WilsonExpClover` action and the `MixedPrecisionRBPrecCG` solver are used. These defaults can be altered by specifying `CPPFLAGS` as indicated above. For now the Möbius DWF fermion action and the `MixedPrecisionRBPrecBiCGSTAB` are available as alternatives. (Note that the DWF action only works with the standard CG solver)

To execute the program run
``` bash
mpirun -n N ./hadrons_npr ./input_file.xml --grid Nx.Ny.Nx.Nt --mpi Mx.My.Mz.Mt
```
where `Nx.Ny.Nx.Nt` is the geometry of the lattice, `N` is the number of MPI processes and `Mx.My.Mz.Mt` describes the MPI process decomposition such that the product of all four components equals `N`.

## Database
The application tracks the exported files in an SQLite database. To get the in- and out- going momenta for every bilinear vertex one can for example query
```sql
SELECT bi.traj, bi.filename, ex1.momentum AS momIn, ex2.momentum AS momOut
  FROM bilinear bi
       LEFT OUTER JOIN externalLeg AS ex1 ON bi.qIn=ex1.qIn
       LEFT OUTER JOIN externalLeg AS ex2 ON bi.qOut=ex2.qIn
```
