=========
Ensembles
=========

Available ensembles are the straightforward NVT, NPT and Gibbs, see :doc:`parapage`
for details of their configuration.
An experimental (ie not fully tested) semi-grand canonical ensemble is available. Also Ostwald ripening can be 
simulated, see below.

NVT
===

Activated by setting ``simulation_type`` to "NVT" (Default).

NPT
===

Activated by setting ``simulation_type`` to "NPT".
``reference_pressure`` should be set to the imposed pressure.

Gibbs
=====

Activated *de facto* by setting ``box_to_box_probability`` and ``volume_exchange_probability`` to a positive value.

Semi-grand canonical ensemble
=============================

Activated by setting ``simulation_type`` to "SGC".

The radii of the particles can vary but not the total number, this variation is regulated by a tabulated chemical potential.
In addition, if ``update_chempot_table`` is set to ``yes``, the chemical potential table is updated according to a goal distribution,
the user is expected to re-run the simulation until convergence of the mean distribution to this goal distribution.

Additional input files
----------------------

table_fit_keffs.dat
    Contains parameters to the fit of charges inverse Debye length as a function of inverse Debye length.
    It also contains data to adjust the inverse Debye length from density and distribution in accordance with the polydisperse cell model.

chempot_table.dat
    Contains the table (2-column Radius and value) of chemical potential. **Note:** Overwritten at the end of the simulation
    if``update_chempot_table`` is set to ``yes``

distrib_goal.dat
    Contains a table of the goal distribution (only necessary if ``update_chempot_table`` is set to ``yes``).

Additional output files
-----------------------

pcmc_kappa.dat
    Inverse Debye length as a function of step.

pcmc_meandist.dat
    Mean size distribution over a run.

pcmc_intermpsd.dat
    Size distribution as a function of step. (2-column (R, D) size distribution separated by a '&' line, grace
    understands it).

chempot_table.dat.old
    If ``update_chempot_table`` is set to ``yes``, it is a copy of ``chempot_table.dat``.


Additional parameters
---------------------

Additional parameters are added in a section ``[SGC]`` in the config file.

particle_resize_prob
    Probability to make a particle resize move. Default is 0.

particle_resize_cvf_prob
    Probability to make a particle resize move at constant volume fraction (involving two particles). Default is 0.


big_step_kappa_update
    Activates the periodic update of the inverse Debye length. Default is ``yes``.

kappa_update_period
    Period of the update of the inverse Debye length. Default is 100.

min_radius
    Minimum reachable particle radius. Default is 2.0.

max_radius
    Maximum reachable particle radius. Default is 30.0.

min_radius_psd
    Minimum particle radius of the calculated size distribution. Default is 2.0.

max_radius_psd
    Minimum particle radius of the calculated size distribution. Default is 30.0.

psd_dr
    Interval of the calculated size distribution. Default is 0.1.

update_kappa_radius_resize
    Activates the update of the inverse Debye length for **each** particle resize. Default is ``no``



Ostwald ripening
================


Activated by setting ``simulation_type`` to "Ostwald".

The particles can exchange “monomers” with a finite bulk. This is meant to be an *out-of-equilibrium* simulation, 
the interpretation in the context of MC simulations is left to the user.


