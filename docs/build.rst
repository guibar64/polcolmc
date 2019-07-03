========
Building
========

Requirements
============

- A recent enough Fortran compiler implementing the 2008 standard (tested with GNU ≥ 7.3 , Intel ≥ 18.0.2)

- CMake ≥ 2.8 or FoBiS_
- (Optional) `fftw3 <http://fftw.org>`_

.. _FoBiS: https://github.com/szaghi/FoBiS

CMake
=====

Default
-------

Create a build directory and change to it:

.. code::

    mdkir build
    cd build


configure :

.. code::

    cmake /path/to/polcolmc/


Build:

.. code::

    make


Install (cmake's default location, e.g. ``/usr/local/bin`` on linux):

.. code:: bash
    
    make install

Options
-------

You can change the install location by changing the ``CMAKE_INSTALL_PREFIX`` variable, for example to ``/foo/bin/`` :

.. code::

    cmake /path/to/polcolmc/ -DCMAKE_INSTALL_PREFIX=/foo/ ...


Optionally, the ffw3 library can be used to compute the structure factors.
It is much faster but not critical as it is only used at the end of the simulation.
Define ``FOURIER_TRANSFORM`` to FFTW3:

.. code::

    cmake /path/to/polcolmc/ -DFOURIER_TRANSFORM=FFTW3


If the header file and library file of fftw3 are not in the default system path,
change the ``CMAKE_PREFIX_PATH`` accordingly:

.. code::

    cmake /path/to/polcolmc/  -DCMAKE_PREFIX_PATH=/path/to/fftw3/

If you want a debug build, set the ``CMAKE_BUID_TYPE`` variable to ``Debug`` (default is ``Release``):  

.. code::

    cmake /path/to/polcolmc/ -DCMAKE_BUID_TYPE=Debug


A rather experimental OpenMP support can be activated by setting the ``OPENMP`` variable to ``true``. 
It has not been found to give significant speed up when the cutoff radius of the interactions is small and the cell decomposition is used (indeed the computation may be slower).

FoBiS
=====

.. code::

   FoBiS.py build


