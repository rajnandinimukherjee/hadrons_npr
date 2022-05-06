# NPR production code

Computes propagators and vertex functions for non-perturbative renormalization in regularization independent momentum-subtraction schemes with MOM and SMOM kinematics.
The program can optionally compute fourquark vertex functions or QED_L corrections to first order in the electromagnetic coupling constant.

## Execution
Install [Grid](https://github.com/paboyle/Grid)
and [Hadrons](https://github.com/aportelli/Hadrons) on their `develop` branch.
Then compile via:

``` bash
./bootstrap.sh
mkdir build
cd build
../configure --with-hadrons=<Hadrons install prefix>
make
```

To execute the program run
``` bash
mpirun -n N ./hadrons_npr ./input_file.xml --grid Nx.Ny.Nx.Nt --mpi Mx.My.Mz.Mt
```
where `Nx.Ny.Nx.Nt` is the geometry of the lattice, `N` is the number of MPI processes and `Mx.My.Mz.Mt` describes the MPI process decomposition such that the product of all four components equals `N`.
