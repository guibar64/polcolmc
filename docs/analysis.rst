========
Analysis
========

.. note::

	This page use the following conventions for output files:
	
	*  The prefix is set to `pcmc` (the default)
	*  $b is used to denote the index of a simulation box
	   (meaning there will be multiple files if more than one box)

.. _RDFcalc:

Radial distribution function
============================

Also named g(r).

Switched on by when :ref:`calculation_period` is positive.
The interval of the histograms is controlled by :ref:`rdf_bin_length`.

Output files
------------

pcmc_box$b_rdf.dat:
	Contains the global RDF. Composed of lines of format ``r g(r)``.

pcmc_box$b_rdfall.dat
	Type-Type RDFS. Composed of line of the form:

	===== ================= ================= ===== ================= ================= ===== ================= =====
	r     :math:`g_{11}(r)` :math:`g_{22}(r)` ...   :math:`g_{12}(r)` :math:`g_{13}(r)` ...   :math:`g_{23}(r)` ...
	===== ================= ================= ===== ================= ================= ===== ================= =====

.. _seffcalc:

Effective structure factor
==========================

Aka S(q)

Calculated from radial distribution functions,
assuming the particles are homogeneous spheres.
A "full" S(q) is first computed with a discrete sine transformation,
it is then interpolated (see :ref:`wave_number_interval`) to get
the final S(q).

Output files
------------

pcmc_box$b_seff.dat
	Contains the effective structure factor as a function of the wave number.
	Composed of lines of format ``q S(q)``.

pcmc_box$b_sfull.dat
	Contains the **full** effective structure factor as a function of the wave number.
	Composed of lines of format ``q S(q)``.


Swaps
=====

When a swaps between two particles is accepted, a counter is incremented
for each pair (type of particle #1, type of particle #2). This allows
to produce a "map" of relative acceptance of swaps between particle type pairs.
For example, this can allow to know if a swaps are performed between particles
of not two close sizes, otherwise swaps are not really useful for this simulation.

Output files
------------

pcmc_swapmap.dat
	Contains relative acceptance *a* of swaps between particles of two different types.
	Composed of lines of the form

	=========== =========== ==============
	:math:`R_i` :math:`R_j` :math:`a_{ij}`
	=========== =========== ==============


	where :math:`R_i` is the radius of particles of type *i*.

pcmc_swapmap.dat
	Same as ``pcmc_swapmap.dat`` but for swaps betwenn particles of *different*
	simulation boxes. See :ref:`interbox_swap_proba`.


.. _bondorderanalysis:

Bond-order parameters
=====================

Bond-order analysis is invoked from the command line like this::

    polcolmc bondorder INPUT_FILE [ --bondorder-cfg BOND_ORDER_CONFIG ]

Options
-------

-d FILE           Get distribution from file *FILE* (default: distrib)

-p FILE           Get parameters from file *FILE* (default: polcolmc.cfg)

--box BOX         (convert) Index of the simulation box to consider (default: 1). 

--bondorder-cfg   BOND_ORDER_CONFIG   Separate configuration file for bond-order analysis

The bond-order parameters (qₗ) are computed from the input trajectory file ``INPUT_FILE``
(eg ``pcmc_out.inst``) and used to determine the composition of crystal phases.

This will produce (too) many files with the ``bondorder_`` prefix, see below.

In the output, phases are referred by the following shortcuts:

========= =====================
*phase*    Description
========= =====================
fcc       FCC
bcc       BCC
hcp       HCP
big       Laves(big)
small     Laves(small)
flu       fluid
========= =====================


Configuration
-------------

In 'polcolmc.cfg'(or file set by option ``-p``), or in a separate file, 
use a [bondorder] section.

    Example:

    .. code:: ini

        [bondorder]
        # typically determined from the g(r)
        cutoff = 55.800000
        cutoff_big = 63.388800
        cutoff_small = 52.731000
        # qₗ limits to separate the phases
        # Laves
        barq6_small_sup = 0.092
        barq6_big_inf = 0.092
        barq6_big_sup = 0.177
        q6_big_sup = 0.273
        q6_small_inf = 0.0354
        # others
        barq6_bcc_inf = 0.308
        barq6_bcc_sup = 0.486
        barq4_bcc_inf = 0.
        barq4_bcc_sup = 0.0617
        barq6_hcp_inf = 0.308
        barq6_hcp_sup = 0.486
        barq4_hcp_inf = 0.0602
        barq4_hcp_sup = 0.124
        barq6_fcc_inf = 0.476
        barq6_fcc_sup = 0.638
        barq4_fcc_inf = 0.124
        barq4_fcc_sup = 3.5
        barq6_flu_inf = 0.
        barq6_flu_sup = 0.308
        barq4_flu_inf = 0.
        barq4_flu_sup = 0.0617

    Parameters of the form (bar)qX_XXX_(inf|sup) determines the bounds 
    of qₗ values for each phase.    

    ======  ===================================
     bar     average qₗ (if absent, regular qₗ)
     XXX     phase
     inf     lower bound
     sup     upper bound
    ======  ===================================

Output files
------------

bondorder_log.txt
	Summary of the main results of the computation.

bondorder_phase_fractions.dat
	first column is volume fraction, other columns give the fraction of each 
	phase:

	=== ==== === ==== ====== ========== ============
	Φ   bcc  hcp fcc  fluid  laves(big) laves(small)
	=== ==== === ==== ====== ========== ============

bondorder_phase_fractions_traj.dat
  Fraction of each phase for each step. Columns of the form:

  ======= ==== ==== === ==== ====== ========== ============
  radius  all  bcc  hcp fcc  fluid  laves(big) laves(small)
  ======= ==== ==== === ==== ====== ========== ============

bondorder_histo_qx.dat
	Contains an histogram of the *qx* parameter where *qx* can be: 
	q4, q4, bq4, bq6, the last two being the neighbor-averaged q's.

bondorder_mean_pops.dat
	Size distribution of the phases *phase*.

bondorder_sel_phase.xyz
	Coordinates of particles in the phase *phase*, in the XYZ format.


---------------------------------------

.. warning::

	The calculations described below are experimental or broken.



Energy ergo analysis
====================

Switched on by :ref:`energy_ergo_analysis`.

Each :ref:`energy_output_period`, it computes the mean particle energy for type *j* and the variance of particle energy.

The mean energy at step *s* (:math:`n_s` sampled steps) for particle *i*:

.. math::

	 \bar{e}_i =  \frac{1}{n_s} \sum_{s'=0}^{s} e_i(s'),

where :math:`e_i(s)` is the potential energy of particle *i* at step *s*.

The mean energy at steps *s* for type *j* is therefore:

.. math::

    \bar{E}_j(s) = \frac{1}{N_j} \sum_{i \in J} \bar{e}_j

The variance at step *s* is

.. math::

  \sigma_j(s) = \frac{1}{N_j} \sum_{i \in J}
  (e_i(s) - \bar{E_j}(s))^2


Output file
-----------

pcmc_box$b_nrgergo.dat
	Mean and variance of potential energy per particle (for each type) computed at several steps.

	=============== ================= ================= ===== ============= ===== 
	step (relative)	variance (type 1) variance (type 1)  ...  mean (type 1)  ...  
	=============== ================= ================= ===== ============= =====


Switched on by 

Configurational temperature
===========================

Switched on by :ref:`configurational_temperature`.

Calculates the configurational temperature :math:`T_c`.

.. math::

	T_c = \frac{\langle \vec{\nabla} V^2 \rangle}{k_B \langle \Delta V \rangle}

Output file
-----------

pcmc_conftemp.dat
	Contains the relative configurational temperature :math:`T_c/T` , squared force F² and second derivate
	potential V'' for each box and for each step. 

	===== ====================== ============= =========== ======
	step   `T_c/T`:math: (box 1)  F² (box 1)    V''(box 1)   ...
	===== ====================== ============= =========== ======

Average Z-density
=================

Calculates average particle density along the z axis for each particle type.
It is controlled by :ref:`zdens_period`.

Output files
------------

pcmc_box$b_zdens.dat
	Density as a function of the position of the interval sampled.
	Composed of lines of the form:

	===== ==== ==== ====
	z     d₁   d₂   ... 
	===== ==== ==== ====

	where dᵢ is the density for type *i*.

3D S(q)
=======

Computes directly :math:`S(\vec{q})` from input :math:`\vec{q}`.
Period of the calculation is set by :ref:`full_structure_period`.

Input file
----------

qvals.dat
	Contains lines of the form ``q_x q_y q_z``.


Output file
-----------

pcmc_bigseff.dat
	contains the structure factor as a function of :math:`\vec{q}`.
	Composed of lines of the form:

	=== === === ===== ====== 
	q_x q_y q_z  S(q)  P(q)
	=== === === ===== ======

	P(q) is the form factor.


Density fluctuations
====================

Fluctuation of the number of particles in small volumes.
Activated by :ref:`density_fluctuations`.
Period of the calculation is set by :ref:`fluctuations_period`.

Input file
----------

fluctuations_input.txt
    Can be set by the parameter :ref:`density_fluctuations_file`.
    
    Format:

    .. code::

        {random seed}

        {number of tests}

        {number of radii}
        {option (RR ⇒ use relative radius, R ⇒ no)}
        {radius 1}
        ...
        {radius {number of radii}}

        {minimum density}
        {maximum density}
        {number of densities}

Output files
------------

pcmc_denses.dat
	Grace file containing the histograms of densities (XY sets)
	for each small volume.

pcmc_rad_denses.dat
	Grace file containing the histograms of mean particle radius
	(XY sets) for each small volume.

pcmc_phi_denses.dat
	Grace file containing the histograms of mean particle volume
	fraction (XY sets) for each small volume.

pcmc_cov_denses.dat
	Contains covariances of particle numbers.

	.. code::
	 	
	 	Fluctuations of families: cov(Ni,Nj)
 		Nb radii:" ${number of small volumes}
  		Nb families: ${number of families}
  		${volume radius #1}
  		1 ${cov(1,1)}
  		2 ${cov(2,2)}
  		...
  		1 2 ${cov(1,2)}
  		....
  		2 3 ${cov(1,2)}
  		...
  		${volume radius #2}
  		...

