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

Activated by setting ``simulation_type`` to "SGVT" or the isobaric variant "SGPT".

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

Additional parameters are added in a section ``[sgc]`` in the config file.

particle_resize_prob
    Probability to make a particle resize move. Default is 0.

particle_resize_cvf_prob
    Probability to make a particle resize move at constant volume fraction (involving two particles). Default is 0.


big_step_kappa_update
    Activates the periodic update of the inverse Debye length. Default is ``yes``.

kappa_update_period
    Period of the update of the inverse Debye length. Default is 100.  Desactivate with a negative value.

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


Activated by setting ``simulation_type`` to "Ostwald-VT" or "Ostwald-PT".

The particles can exchange “monomers” with a finite bulk. This is meant to be an *out-of-equilibrium* simulation, 
the interpretation in the context of MC simulations is left to the user.


Additional input files
----------------------

table_fit_keffs.dat
    Contains parameters to the fit of charges inverse Debye length as a function of inverse Debye length.
    It also contains data to adjust the inverse Debye length from density and distribution in accordance with the polydisperse cell model.

Additional output files
-----------------------

pcmc_kappa.dat
    Inverse Debye length as a function of step.

pcmc_eq_rmean.dat
    Mean radius and others as a function of step.

    ==== ================== ================ ================== =============================
    step mean radius        critial raidus   variance of radius moves/cycle/(nb of particles)
    ==== ================== ================ ================== =============================

    The last column may be useful to recalculate an effective step
    because the nb of particles is variable.

pcmc_meandist.dat
    Mean size distribution over a run.


pcmc_intermpsd.dat
    Size distribution as a function of step. (2-column (R, D) size distribution separated by a '&' line, grace
    understands it).

pcmc_psdrel.dat
    Relative size distribution as a function of step

pcmc_dFormvsdEint.dat
    Average energies of dissolution/precipitation as a function of radius.

    = =========== ============ =========== ============
    R ΔE_int_diss ΔF_form_diss ΔE_int_prec ΔF_form_prec
    = =========== ============ =========== ============

Additional parameters
---------------------

Additional parameters are added in a section ``[sgc]`` in the config file.

particle_bulk_exchange_prob
    Probability to make move that exchange monomers between the solid particles and the bulk. Default is 0.

particle_bulk_formation_prob
    Probability to make move that precipitates/dissolves a particle from/to the bulk. Default is 0.

total_number_monos
    Sets the total number of monomers. Default is 300 000.

initial_number_dissolved_monos
    Sets the initial number of dissolved monomers if a positive value is given. It takes over ``total_number_monos``.

minimum_monos_exchange
    Minimum of monomers that can be exchanged with the bulk. Default is 10.

maximum_monos_exchange
    Maximum of monomers that can be exchanged with the bulk. Default is 100.

minimum_monos_particles
    Minimum of monomers taken/given from/to the bulk for a precipitation/dissolution. Default is 10.

maximum_monos_particles
    Maximum of monomers taken/given from/to the bulk for a precipitation/dissolution. Default is 100.

gamma0
    Surface tension (in kT) for an infinite radius and a neutral surface. Default is 17.

tolman_length
    Sets the Tolman length. Default is 1.0.

gamma_doublelayer_inf
    Sets the double layer contribution to the surface tension for an infinite radius. Default is -0.633177.

gamma_doublelayer_slope
    Sets the slope regarding 1/R of the double layer contribution to the surface tension. Default is 1.535432.


big_step_kappa_update
    Activates the periodic update of the inverse Debye length. Default is ``yes``.

kappa_update_period
    Period of the update of the inverse Debye length. Default is 100. Desactivate with a negative value.

min_radius_psd
    Minimum particle radius of the calculated size distribution. Default is 2.0.

max_radius_psd
    Minimum particle radius of the calculated size distribution. Default is 30.0.

psd_dr
    Interval of the calculated size distribution. Default is 0.1.

update_kappa_radius_resize
    Activates the update of the inverse Debye length for **each** particle resize. Default is ``no``.

min_radius_gamma
    Below that value. the surface tension is constant. Default is 0.0

pseudo_kinetics
    Toggles the use of a kinetic factor for exchanges with bulk. **Not really working**. Default is ``no``.

implicits_maximum_size
    Maximum monomer size of the smallest particles. Default is 100.

implicits_minimum_size
    Minimum monomer size of the smallest particles. Default is 2.

solid_density
    Density of the solid phase (ie, particles). Default is 22.96 nm^-3.

pKs
    Decimal logarithm of the solubility constant. Default is 2.7.

finite_monomer_reservoir
    If set to ``no``, the bulk is considered infinite. Default is ``yes``.

bulk_exchange_create_destroy
    Toggles the creation/deletion of particles in the main exchange move. Default is ``no``.
