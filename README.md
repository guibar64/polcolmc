# polcolmc

Monte-Carlo simulations of polydisperse spheres.

<!-- Detailed build and usage instructions are in the online [documentation]. -->

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

This software is released into the public domain. See `LICENSE`.
