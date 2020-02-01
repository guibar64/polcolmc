================
Input Parameters
================


Some symbol definitions
=======================

-  `\lambda_B`:math: : Bjerrum_length
- `\lambda_D`:math: : Debye_length
- `\beta = \frac{1}{k_B T}`:math:
-  `k_B`:math: : Boltzmann constant
- `T`:math:   : temperature
- `Q_i`:math: : charge of particle type i (defined in the distribution file)
- `R_i`:math: : radius of particle type i (defined in the distribution file)
- `u_{ij}`:math: : pair potential between particle types i and j


Potentials
==========

potential_type
--------------

Can be:

- *none*: only hard-core

-   *yukawa* [default]

    .. math::
        
        \beta u_{ij}(r) = \frac{\lambda_B}{r} \frac{Q_i Q_j \exp[- (r-R_i-R_j)/\lambda_D]}{(1+R_i/\lambda_D)((1+R_j/\lambda_D)}
	    
- *yukawa2*: tabulated yukawa, faster, but may not be accurate in some ranges.
 
- *fennel*: fennel potential.

    .. math::

        \beta u_{ij}(r) = \lambda_B Q_i Q_j \left( \frac{1}{r} + \frac{r}{r_c^2} - \frac{2}{r_c} \right)
  

energy_algorithm
----------------
  
Can be:

-  *standard* [default]: loop over all particles

- *celldec* (cell decomposition, generally faster).

Debye_length
------------

Quantifies the range of electrostatic interactions.
Unit is nm, default is 4.31.

Bjerrum_length 
--------------

Defined by:

.. math::
    \lambda_B = \sqrt{\frac{\beta e^2}{4 \pi \epsilon_0 \epsilon_r }}

Unit is nm, default is 0.7105


cutoff_radius
-------------

A pair of particles which separation is above the cutoff radius is not taken into account
for calculating the potential. Unit is nm. A negative value means that
`cutoff_estimation_tolerance`_ is used. Default is -1.

.. _cutoff_estimation_tolerance:

cutoff_estimation_tolerance
---------------------------

If set to a positive value, the cutoff is adjusted to the separation where the potential (in kT)
between members of the last family (in ``distrib``) goes above this factor. This option is ignored if `cutoff_radius`_ is set to a positive value.

.. warning::

    This is convenient for e.g. Yukawa potentials where the last family interacts
    the strongest (biggest particles). If not, adjust directly `cutoff_radius`_.


delta_energy_tolerance
----------------------

The relative tolerance when the energy is recalculated and checked against the current energy of the simulation.
If the absolute relative difference between these two energies exceeds this value, the program stops. Default is 1e-5.

Structural parameters
=====================

box_type 
--------

Can be:

- *generic* [default]: triclinic box, periodic boundary conditions

- *cubic*: cubic box,  periodic boundary conditions

- *cubicXY*: cubic box, periodic conditions only along the X and Y axis
  (hard walls along the Z axis)

- *hexagonal*: hexagonal box,  periodic boundary conditions

box_length_X 
------------

Box length along the x-axis. Can be overridden by the input configuration. Default is 1.

box_length_Y
------------

Box length along the x-axis. Can be overridden by the input configuration. Default is 1.

box_length_Z 
------------

Box length along the x-axis. Can be overridden by the input configuration. Default is 1.

alpha_box_angle
---------------

(YZ) angle of the box in degrees.
Can be overridden by the input configuration. Default is 90°. 

beta_box_angle
--------------

(ZX) angle of the box in degrees.
Can be overridden by the input configuration. Default is 90°. 

gamma_box_angle
---------------

(XY) angle of the box in degrees.
Can be overridden by the input configuration. Default is 90°.

number_of_boxes 
---------------

Total number of boxes.
Can be overridden by the input configuration. Default is 1.

max_number_of_particles
-----------------------

Total number of particles expected during the simulation (ie it should not go beyond this value)
If negative or too small, this number is the total number of particles of the input distribution. Default is -1.

volume_fraction
---------------

Set volume fraction on input. (Rescales the box lengths).
Only relevant when an initial configuration has to be generated.
A negative value is ignored, that is the box sizes are not rescaled.
Default is -1.

Simulation parameters
=====================

simulation_type
---------------

Can be "NVT", "NPT", "SGC" (semi-grand canonical) or "Ostwald" (dissolution-precipitation). See :doc:`ensembles`.

temperature 
-----------

Unit is K, default is 300.

random_seed
-----------

Seed for the random number generator.
Integer. Default is 666.

production 
-----------

The main run becomes a production run is switched on and conversely.
Default is *yes*.

.. _equilibration_steps:

equilibration_steps
-------------------

If set to a positive number, perform an equilibration run of ``equilibration_steps`` steps
before the main run.
Default is 0.

number_of_steps
---------------

Number of steps for the main run. Default is 1.

initial_step
------------

Set the value of the first step (default: 0). Can be overridden by the initial configuration.


.. warning::

    By convention, *steps* are taken as MC cycles.


stat_on_cells
-------------

If set to *yes*, the number of particles per cell and its standard deviation are computed.


Output
======

.. _output_prefix:

output_prefix
-------------

The prefix prepended to output files. Default is *pcmc*.

log_period 
----------

Period at which some current information is printed in the log file. Default is 0.


snapshot_period 
---------------

Period at which the current configuration is printed to the trajectory file.
If ≤ 0, no configuration but the initial are printed.
Default is 0.


.. _energy_output_period:

energy_output_period
--------------------

Period at which the current energy and pressure is printed to the energy file.
If ≤ 0, no configuration but the initial are printed.
Default is 0.

use_binary_trajectory
---------------------

If *yes*, the trajectory is stored with a binary format
instead of an ASCII one. Default is *no*.

Translations
============

In equilibration runs, the `translation_interval`_ can be adjusted every `translation_update_period`_
so that the acceptance of translations tends to a goal acceptance,
`translation_acceptance`_.

.. note::

    Changing the amplitude of translations during a run may yield an incorrect
    Monte-Carlo sampling so this feature is deactivated during production runs.
    If an equilibration run is made before, the value is kept for the subsequent production run.

translation_update_period
-------------------------

Period at which the translation interval is updated. Default is 100.

translation_interval
--------------------

Initial amplitude of the translations. Can be overridden by the initial configuration.
Unit is nm, default is 3.

translation_acceptance 
----------------------

Goal acceptance for translation moves. Default is 0.5.

translation_interval_wall
-------------------------

Initial maximum amplitude for translations near the wall for simulations with hard boundary conditions on the z-axis.
Unit is nm.
A negative value means that the parameter is set to `translation_interval`_. Default is -1.

Histograms
==========

.. _rdf_bin_length:

rdf_bin_length
--------------

Sets the resolution of the radial distribution functions.
Unit is nm, default is 0.1.

hsp_interval
------------

Space bin used for the hard-core pressure.
Unit is nm, default is 0.05.

.. _wave_number_interval:

wave_number_interval
--------------------

Sets the resolution of the structure factor.
Unit is `\text{nm}^{-1}`:math:, default is 0.002.

.. _calculation_period:

calculation_period 
------------------

Period at which a new configuration is used to update the histograms
used for the calculation of the RDFs.
If ≤ 0, no configuration is used but the initial one is used.

.. _zdens_period:

z_density_period
----------------

Set the number of steps between each sampling of the mean density along the z-axis.
A negative value deactivate the calculation. Default is -1.

.. _energy_ergo_analysis:

energy_ergo_analysis
--------------------

Turn on an ergodicity analysis based on correlation between mean particle energy.

Swaps
=====

In a swap move two particles selected at random swaps theirs positions.

swap_probability
----------------

The probability to perform a swap.

NPT
===

reference_pressure
------------------

The pressure P in NPT.

volume_change_probability 
-------------------------

The probability to perform a volume change move.

deltavolume 
-----------

The amplitude of relative volume difference for a volume change move.

standard_volume_change 
----------------------

If set to *yes*, for a volume change the volume is multiplied by a random factor. 
If not, it is the relative logarithm.
Default is *no*.

Other moves (Gibbs...)
======================

.. note::

    By default the following moves are deactivated, that is their
    probability is set to 0.

.. _interbox_swap_proba:

interbox_swap_probability
-------------------------

Probability of a move consisting in a swap between particles of different boxes.
Default is 0.


box_to_box_probability
----------------------

Probability of a move consisting in moving a particle from one box to another.

volume_exchange_probability
---------------------------

Probability for a move consisting in decreasing the volume of a box and increasing the volume
of another box by the same amount. The total volume is preserved.

exchange_deltavolume
--------------------

The maximum of relative volume difference for a volume exchange between boxes.

volume_particle_exchange_probability 
------------------------------------

Probability of a move which is a composition of a transfer of a particle from one box
to another and a volume exchange between the same boxes.

on_the_fly_chemical_potential
-----------------------------

When particle can go from a box to another, this tries to calculate the chemical
potential by the Widom method opportunistically when such move is done.

particle_exchange_rosenbluth_probability 
----------------------------------------

Probability of a biased move aiming at transferring a particle to another box.

.. note:: 

    A series of trial positions in the new box is generated, each with a weight equal to
    `w = \exp(-\beta E(trial)-E(old))`:math:. The finally attempted position is selected
    at random according to the weights. As this procedure introduce a bias,
    weights of  the old configuration and weights of the new are computed
    and injected into the acceptance ratio to correct the bias.

number_of_rosenbluth_trials 
---------------------------

The number of trials for the special move described above. Default is 10.

maximum_family_number_rosenbluth 
--------------------------------

The maximum family number for the special move described above.
Useful for restricting this move to a subset of the particles.

impure_interbox_swap_probability 
--------------------------------

Probability of a swap between a special kind of particle in one box and another kind of particle
in another box.

.. note:: 

    This move attempts to swap a particle with an "impurity"
    (i.e. the first particle family of particle, supposedly very small)
    between two different boxes.
    This was introduced to accelerate the equilibration of dense systems.
    Indeed, a very small particle (in small quantities, this can be taken as
    as an impurity) can move more easily from one box to another.
    Once there, it can help bigger particles to pass to a different box by swapping
    their positions. In practice, did not work.

interbox_swap_impure_rosbth_probability 
---------------------------------------

Probability for the same move above with repositioning of the swapped particle.

.. note::

    Another attempt to swap a particle with an "impurity"
    (ie the first particle family of particle, supposedly very small)
    between two different boxes, with Rosenbluth trials.
    Same as the previous one, with Rosenbluth trials.
    A series of positions for the non-impurity is generated, each
    deviating from the original position of the impurity within some radius. Each position has a weight
    `w = \exp(-\beta E(trial)-E(old))`:math:. The finally attempted position is
    selected at random according to its weight. The introduced bias is corrected
    by computing the total weight of the old position and the new series of positions,
    these quantities modify the acceptance factor accordingly.


interbox_swap_impure_rosbth_radius 
----------------------------------

For the special move described above, it sets the radius around the impurity selected.

interbox_swap_impure_rosbth_2_probability 
-----------------------------------------

Probability for the same move above with repositioning of the particles around.

.. note::

    Yet another attempt to swap a particle with an "impurity"
    (i.e. the first particle family of particle, supposedly very small)
    between two different boxes, with Rosenbluth trials (second method).
    Same principle that previously,
    but in addition the neighboring particles are virtually taken out the box
    and replaced into it at random around the new position of the swapped
    particle.

Others
======

.. _configurational_temperature:

configurational_temperature
---------------------------

Turn on calculations of the configurational temperature by
setting it to *yes*. Default is *no*.

.. _full_structure_period:

full_structure_period 
---------------------

Period of calculations of the 3D structure factor.
A negative value deactivate this calculation. Default is -1.

.. _density_fluctuations:

density_fluctuations
--------------------

Turn on calculations of the density fluctuations by
setting it to *yes*. Default is *no*.

.. _fluctuations_period:

fluctuations_period
-------------------

period of calculation of fluctuations. 

.. _density_fluctuations_file:

density_fluctuations_file
-------------------------

Sets input file for density fluctuations calculations.
Default is *fluctuations_input.txt*.

