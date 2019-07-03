# polcolmc

Monte-Carlo simulations of polydisperse spheres.

For detailed build and usage instructions, see the online [documentation](https://polcolmc.readthedocs.io/en/latest).

## Minimum requirements

- A Fortran compiler implementing the 2008 standard
- CMake ≥ 2.8

## Build

In the source directory, create a build directory, change to it, and run cmake
```
mkdir build
cd build
cmake ../
make
```

Install (default location, ex. `/usr/local/bin` on linux):

```
make install
```


## Basic usage

Assuming polcolmc is in your path,
```
polcolmc
```
will run a simulation in the working directory.


## Licensing

This software is released into the public domain. See `LICENSE.txt`.
