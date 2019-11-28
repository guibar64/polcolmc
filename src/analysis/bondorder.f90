
module bondorder
!! bond order parameter calculation and analysis
use state, only: Distrib, Box, PI
use geom, only: volume_box
use spherharm
use iso_fortran_env
implicit none


integer, private, parameter :: lOrder=6, lOrder2=4

integer, private, parameter :: max_vois = 200
type OrderData
  !! Contains data for bon-order calculations *for each* particle
  real(8) :: pos(3)
  integer :: fam
  complex(8) :: q6v(-lOrder:lOrder), q4v(-lOrder2:lOrder2), barq4v(-lOrder2:lOrder2), barq6v(-lOrder:lOrder)
  real(8) :: q6, order_moy, q4, order2_moy, barq6, barq4
  integer :: nvois,indic, qtype , npartconn
  integer, dimension(max_vois) :: voisin
  real(8), dimension(3, max_vois) :: uvois
  real(8), dimension(max_vois) :: dist
  logical :: onprend, onprend_nb
  integer :: n_sel_nb
   
  logical :: isPreBig, isPreSmall, isBig, isSmall, extracted
  integer :: n_big_nb
end type OrderData

character(*),parameter :: fmt_519="(a,1x,i6,2x,a,1x,a,1x,i5,4x,f8.3,f8.3,f8.3,1x,f5.2,f6.2)"
character(len=:), allocatable, private :: prefout
integer, private :: analog=71

real(8), private, parameter :: s4p = sqrt(4*PI)
real(8), private :: norm_fac = s4p, rc = -1.

logical, private :: alternative_calc = .false., compute_q4q6_map = .false.
real(8), private :: dq=0.005, q_min,q_max
real(8), private :: rcin=0.d0, rc_big = -1, rc_small = -1, rc_large, rc_sel = -1.0
integer, private :: mconn = 8
real(8), private :: distConn(32), CorrQThrsh=0.2
real(8), private :: vol_frac
real(8), private :: q6_big_sup=0.077*s4p,q6_small_inf=0.1*s4p
real(8), private :: barq6_small_sup=0.02597*s4p, barq6_big_inf=0.02598*s4p, barq6_big_sup=0.05*s4p
real(8), private :: barq6_bcc_inf=0.087*s4p, barq6_bcc_sup=0.137*s4p
real(8), private :: barq4_bcc_inf=0.d0*s4p, barq4_bcc_sup=0.0174*s4p
real(8), private :: barq6_hcp_inf=0.087*s4p, barq6_hcp_sup=0.137*s4p
real(8), private :: barq4_hcp_inf=0.0174*s4p, barq4_hcp_sup=0.035*s4p
real(8), private :: barq6_fcc_inf=0.135*s4p, barq6_fcc_sup=0.18*s4p, q4_fcc_inf=0.158, q4_fcc_sup=0.25
real(8), private :: barq4_fcc_inf=0.035*s4p, barq4_fcc_sup=1.0*s4p
real(8), private :: barq6_flu_inf=0.d0*s4p, barq6_flu_sup=0.087*s4p
real(8), private :: barq4_flu_inf=0.d0*s4p, barq4_flu_sup=0.0174*s4p

real(8), private :: first_sep = -1.0, bcchcp_sep_factor = 1.29d0
real(8), private :: hcp_alt_min_bq4 = 0.04d0, bcc_alt_max_bq4 = 0.064d0, fcc_alt_min_bq4 = 0.094d0
real(8), private :: flu_alt_max_bq6 = 0.276


private :: pfil, calc_vol_frac, frac, smalls_count, bigs_count
private :: is_selected, number_of_nb_sel
private :: cond_bcc, cond_hcp, cond_fcc, cond_flu, q6_cond_1, q6_cond_2
private :: barq6_cond_1, barq6_cond_2, big_q6_bq6_cond, small_q6_bq6_cond
private :: nbvois_cond_small, nbvois_cond_big
private :: makesel, calc_mean_sel_nb, calc_mean_nb, calc_all_n_sel_nb, select_crit_custom
private :: prodm3v3, pbc
private :: init_cutoffs, do_calc_sep_phases, ze_selection_big_small_small_around_big
private :: output_bigsmallbig, calc_qs, calc_averageq_around_nb, output_pop, search_neighbours
private :: output_histo, write_conf_pdb_withbarq

abstract interface
  function OrderCondition(item)
    import
    logical :: OrderCondition
    type(OrderData), intent(in) :: item
  end  function
end interface

contains


  subroutine bondorder_initialize(np, xv, dist, config)
    !! Initialize bond-order analysis
    use readparam
    implicit none
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np) !! particle data array
    type(Distrib) :: dist     !! distribution
    type(ParamList) :: config  !! parameters
    character(len=:), allocatable :: key, val
    integer :: i

    prefout = "bondorder_"
    call start_iter_keyval(config)
    do while(iter_keyval_in(config))
      call get_current_keyval(config, key, val)
      select case(key)
      case("bondorder.cutoff_in")
        read(val, *) rcin
      case("bondorder.cutoff")
        read(val, *) rc
      case("bondorder.cutoff_sel")
        read(val, *) rc_sel
      case("bondorder.cutoff_big")
        read(val, *) rc_big
      case("bondorder.cutoff_small")
        read(val, *) rc_small
      case("bondorder.dq")
        read(val, *)  dq
      case("bondorder.q6_small_inf")
        read(val, *) q6_small_inf
      case("bondorder.q6_big_sup")
        read(val, *) q6_big_sup
      case("bondorder.barq6_small_sup")
        read(val, *) barq6_small_sup
      case("bondorder.barq6_big_inf")
        read(val, *) barq6_big_inf
      case("bondorder.barq6_big_sup")
        read(val, *) barq6_big_sup
      case("bondorder.barq4_bcc_inf")
        read(val, *) barq4_bcc_inf
      case("bondorder.barq4_bcc_sup")
        read(val, *) barq4_bcc_sup
      case("bondorder.barq6_bcc_inf")
        read(val, *) barq6_bcc_inf
      case("bondorder.barq6_bcc_sup")
        read(val, *) barq6_bcc_sup
      case("bondorder.barq4_hcp_inf")
        read(val, *) barq4_hcp_inf
      case("bondorder.barq4_hcp_sup")
        read(val, *) barq4_hcp_sup
      case("bondorder.barq6_hcp_inf")
        read(val, *) barq6_hcp_inf
      case("bondorder.barq6_hcp_sup")
        read(val, *) barq6_hcp_sup
      case("bondorder.barq4_fcc_inf")
        read(val, *) barq4_fcc_inf
      case("bondorder.barq4_fcc_sup")
        read(val, *) barq4_fcc_sup
      case("bondorder.barq6_fcc_inf")
        read(val, *) barq6_fcc_inf
      case("bondorder.barq6_fcc_sup")
        read(val, *) barq6_fcc_sup
      case("bondorder.barq4_flu_inf")
        read(val, *) barq4_flu_inf
      case("bondorder.barq4_flu_sup")
        read(val, *) barq4_flu_sup
      case("bondorder.barq6_flu_inf")
        read(val, *) barq6_flu_inf
      case("bondorder.barq6_flu_sup")
        read(val, *) barq6_flu_sup
      case("bondorder.q4_fcc_inf")
        read(val, *) q4_fcc_inf
      case("bondorder.q4_fcc_sup")
        read(val, *) q4_fcc_sup
      case("bondorder.min_correlation")
        read(val, *) CorrQThrsh
      case("bondorder.min_connections")
        read(val, *) mconn
      case("bondorder.norm_factor")
        read(val, *) norm_fac
      case("bondorder.output_prefix")
        prefout = val(1:len_trim(val)) // "_"
      case("bondorder.alternative_calc")
        alternative_calc = read_logical(val)
      case("bondorder.compute_q4q6_map")
        compute_q4q6_map = read_logical(val)
      case("bondorder.separation_first_neighbors")
        read(val, *) first_sep
      case("bondorder.bcchcp_sep_factor")
      read(val, *) bcchcp_sep_factor
      case default
        if(key(1:min(len(key),10)) == "bondorder.") then
          write(error_unit, *) "Warning: unkown parameter '",&
            key(1:len_trim(key)),"'"
        end if
      end select
      call param_list_next(config)
    end do
   
   ! nbref=4*PI/3*((rc/equivLength)**3)*np
    distConn = 0.D0

    do i=1,np
      xv(i)%order_moy = 0.d0
      xv(i)%order2_moy = 0.d0
    end do

    open(newunit=analog, file=pfil("log.txt"))
    
  end subroutine bondorder_initialize

  function pfil(s)
    !! Adds a prefix to a string (typ. file name).
    character(len=:), allocatable :: pfil
    character(*) :: s
    pfil = prefout // s
  end function pfil
  
  subroutine bondorder_finalize(np, xv)
    !! End of bond-order
    integer :: np
    type(OrderData) :: xv(np)
    close(analog)
  end subroutine bondorder_finalize    
  
  subroutine bondorder_do_step(np, xv, b, dist)
    !! Calculation step
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np) !! particle data array
    type(Box), intent(in) :: b !! simulation box (to get positions etc...)
    type(Distrib) :: dist       !! particle distribution
    vol_frac = calc_vol_frac(b)

    xv(:)%onprend = .true.
    
    ! Selections from splitting q6: n_neighbours_crit
    
    call init_cutoffs(b%volume**(1.d0/3.d0), np)
    if(compute_q4q6_map) then
      call search_neighbours(np, xv, b, rc_large, 0.0D0)
      call calc_qs(np, xv, rc, rcin)
      call calc_map(np, xv, b, rc, dq)
    else
      call output_pop(pfil("pop_ref.dat"), np, xv, b, dist)
      if(alternative_calc) then
        call do_calc_sep_phases_alternative(np, xv, b, dist)
      else
        call search_neighbours(np, xv, b, rc_large, 0.0D0)
        call calc_qs(np, xv, rc, rcin)
        call calc_averageq_around_nb(np, xv, rc, rcin)
        write(*,*) "Note: bond-order parameters calculated."
        call output_histo(pfil("histo_bq6.dat"), np, xv(:)%barq6, 0.0D0, 1.D0 , 0.005D0)
        call output_histo(pfil("histo_bq4.dat"), np, xv(:)%barq4, 0.0D0, 0.2D0 , 0.001D0)
        call output_histo(pfil("histo_q6.dat"), np, xv(:)%q6, 0.0D0, 1.0D0 , 0.005D0)
        call output_histo(pfil("histo_q4.dat"), np, xv(:)%q4, 0.0D0, 1.0D0 , 0.005D0)
        
        call write_conf_pdb_withbarq(pfil("conf_withbarq.pdb"), np, xv, [b%tens(1,1),b%tens(2,2),b%tens(3,3)])

        call do_calc_sep_phases(np, xv, b, dist)
      end if
    end if      
  end subroutine
  
  subroutine init_cutoffs(L, np)
    real(8), intent(in) :: L
    integer, intent(in) :: np
    !! Initialize small, big and selection cutoff radii.

    if(L> 0.D0 .and. rc<0.D0) then
      rc=1.5D0*L/(np**(1.D0/3))
      write(analog,*) "cutoff radius = ",rc
    end if

    if(rc_small < 0) rc_small = 0.945 * rc
    if(rc_big < 0) rc_big = 1.136 * rc
    if(rc_sel < 0) rc_sel = rc

    rc_large = max(rc_small, rc_big, rc)
  end subroutine

  real(8) function calc_vol_frac(b)
    !! Computes volume fraction of a box.
    type(Box), intent(in) :: b
    real(8) :: s,v
    integer :: i
    s = 0.0
    do i=1,b%n
      s = s + b%parts(i)%rayon**3
    enddo
    v = volume_box(b%tens)
    calc_vol_frac = 4.d0/3.d0*PI*s/v
  end function
  
  subroutine do_calc_sep_phases(np, xv, b, dist)
    !! Find the phase of particles (current state) among fluid, bcc, fcc, hcp and Laves.
    !! Outputs particle distributions of each phase.
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np) !! particle data array
    type(Box) :: b            !! simulation box
    type(Distrib) :: dist     !! particle distribution
    integer :: nbig, nsma, nbcc, nhcp, nfcc, nflu,nbig_pre, nsma_pre, ufc=72, ntot, nfcc_2, nhcp_2
    real(8) :: rc_mean
    ! mean cut-off (used as well for the computation of order parameters)
    rc_mean = rc
    write(*,*) "calc sep phases" 
    write(*,*) "Laves"
    write(analog,*)  
    ! Select 'Big' particles
    write(analog,*) "Big particles (q6+bq6+nvois)"
    nbig_pre = select_crit_custom(np, xv, "critq6bq6nvois_big", big_q6_bq6_cond, nbvois_cond_big, rc_big, with_nb=.true.)
    !call output_pop(pfil("pop_big-pre.dat"), np, xv, b, dist)
    
    xv(:)%isPreBig = xv(:)%onprend
  
    write(analog,*)
    ! Select 'Small' particles
    write(analog,*) "Small particles (q6+bq6+nvois)"
    nsma_pre = select_crit_custom(np, xv, "critq6bq6nvois_small", small_q6_bq6_cond, nbvois_cond_small, rc_small, with_nb=.true.)
    !call output_pop(pfil("pop_small-pre.dat"), np, xv, b, dist)

    xv(:)%isPreSmall = xv(:)%onprend

    call ze_selection_big_small_small_around_big(np, xv, rc_mean)
    nbig = bigs_count(np,xv)
    nsma = smalls_count(np,xv)
    write(analog,*) "Small and bigs fulfilling: rCutmean & nvoispetites_autour_d_une_grosse"
    write(analog,*) "Nbigs   = ", nbig
    write(analog,*) "Nsmalls = ", nsma
    call output_bigsmallbig(np, xv, b%tens(1,1), 0)

    xv(:)%onprend = xv(:)%isBig
    call output_pop(pfil("pop_big.dat"), np, xv, b, dist)
    xv(:)%onprend = xv(:)%isSmall
    call output_pop(pfil("pop_small.dat"), np, xv, b, dist)

    ! Filter bc/fcc/etc particles...
    write(analog,*) "bcc/fcc/etc particles"
    write(*,*) "bcc/fcc/etc particles" 
    write(analog,*) "bcc particles"
    nbcc = select_crit_custom(np, xv, "bcc", cond_bcc, is_selected,  rc_sel, with_nb=.true.)
    call output_pop(pfil("pop_bcc.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "bcc", b%tens(1,1), 0)

    write(analog,*) "hcp particles"
    nhcp = select_crit_custom(np, xv, "hcp", cond_hcp, is_selected,  rc_sel, with_nb=.true.)
    call output_pop(pfil("pop_hcp.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "hcp", b%tens(1,1), 0)


    write(analog,*) "fcc particles"
    nfcc = select_crit_custom(np, xv, "fcc", cond_fcc, is_selected,  rc_sel, with_nb=.true.)
    call output_pop(pfil("pop_fcc.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "fcc", b%tens(1,1), 0)

    write(analog,*) "fcc particles (q4-bases)"
    nfcc_2 = select_crit_custom(np, xv, "fcc_2", cond_fcc_2, is_selected,  rc_sel, with_nb=.true.)
    call output_pop(pfil("pop_fcc_2.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "fcc_2", b%tens(1,1), 0)

    write(analog,*) "fcc particles (remove q4-based fcc parts)"
    nhcp_2 = select_crit_custom(np, xv, "hcp_2", cond_hcp_2, is_selected,  rc_sel, with_nb=.true.)
    call output_pop(pfil("pop_hcp_2.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "hcp_2", b%tens(1,1), 0)

    write(analog,*) "fluid-like particles"
    nflu = select_crit_custom(np, xv, "flu", cond_flu, is_selected,  rc_sel, with_nb=.true.)
    xv(:)%onprend = .not. (xv(:)%isSmall .or. xv(:)%isBig)
    call output_pop(pfil("pop_flu.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "flu", b%tens(1,1), 0)

    ntot = nbcc + nhcp + nfcc + nflu + nbig + nsma
    write(analog,*)
    write(analog,*) "tot criter = ", ntot, " / ", np

    open(newunit=ufc, file=pfil("phase_fractions.dat"))
    write(ufc, '(9ES16.7)') vol_frac, frac(nbcc,np), frac(nhcp,np), &
    frac(nfcc,np), frac(max(0,nflu-nbig-nsma),np),  frac(nbig, np), frac(nsma,np), frac(nfcc_2, np), frac(nhcp_2, np)
    close(ufc)
  end subroutine do_calc_sep_phases

  subroutine do_calc_sep_phases_alternative(np, xv, b, dist)
    !! Find the phase of particles (current state) among fluid, bcc, fcc, hcp (alternative).
    !! Outputs particle distributions of each phase.
    !!
    !! 1. r_c = 1.4*d_1ngh ⇒ extract fluid with q6
    !! 2. r_c = f2*d_1ngh ⇒ extract fcc
    !! 3. r_c = f2*d_1ngh ⇒ ̇if bq4 < XXX: bcc else if n_neigh(in S(r_c)) ⪆ 17:  
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np) !! particle data array
    type(Box) :: b            !! simulation box
    type(Distrib) :: dist     !! particle distribution
    integer :: i, nbcc, nhcp, nfcc, nflu, ufc=72, nq
    real(8) :: rc_mean
    real(8) :: qvals(np)

    if(first_sep<0) then
      stop "bondorder: please provide a 'first_step' input parameter"
    end if

    call search_neighbours(np, xv, b, 1.56*rc_large, 0.0D0) 

    call calc_qs(np, xv, 1.4d0*first_sep, 0.d0)
    call calc_averageq_around_nb(np, xv, 1.4d0*first_sep, 0.d0)

    ! mean cut-off (used as well for the computation of order parameters)
    rc_mean = rc
    write(analog,*) "fluid-like particles"

    nflu = select_crit_custom(np, xv, "flu", cond_flu_alt, is_selected, rc_sel, with_nb=.false.)
    where(xv%onprend)
      xv%extracted = .true.
    end where
    call output_pop(pfil("pop_flu.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "flu", b%tens(1,1), 0)
    print *, "flu: ", nflu, nflu, num_extracted(np,xv)
 

    call calc_qs(np, xv, bcchcp_sep_factor*first_sep, 0.d0)
    call calc_averageq_around_nb(np, xv, bcchcp_sep_factor*first_sep, 0.d0)
    
    block
      nq = 0
      do i=1,np
        if(.not. xv(i)%extracted) then
          nq = nq + 1
          qvals(nq) = xv(i)%barq4
        end if
      end do
      call output_histo(pfil("histo_bq4_notflu.dat"), nq, qvals, 0.0D0, 0.2D0 , 0.001D0)

      nq = 0
      do i=1,np
        if(xv(i)%extracted) then
          nq = nq + 1
          qvals(nq) = xv(i)%barq4
        end if
      end do
      call output_histo(pfil("histo_bq4_flu.dat"), nq, qvals, 0.0D0, 0.2D0 , 0.001D0)

      call output_histo(pfil("histo_bq4.dat"), np, xv(:)%barq4, 0.0D0, 0.2D0 , 0.001D0)
    end block
    write(analog,*) "fcc"
    nfcc = select_crit_custom(np, xv, "fcc", is_not_extracted, cond_fcc_alt, rc_sel, with_nb=.false.)
    where(xv%onprend)
      xv%extracted = .true.
    end where
    call output_pop(pfil("pop_fcc.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "fcc", b%tens(1,1), 0)
    print *, "fcc: ", nfcc, nfcc + nflu, num_extracted(np,xv), num_sel(np, xv)

    nhcp = select_crit_custom(np, xv, "hcp", is_not_extracted, cond_hcp_alt, rc_sel, with_nb=.false.)
    print *,"hcp (presel): ", nhcp 
    block
      integer :: nvsel,jv,j
      real(8) :: rcn
      rcn = 1.74*first_sep
      xv(:)%onprend_nb = .false.
      do i=1,np
        if(xv(i)%onprend) then
          if(xv(i)%barq4 > bcc_alt_max_bq4) then
            ! Should be too much to be bcc
            xv(i)%onprend_nb = .true.
          else
            nvsel = 0
            do jv=1,xv(i)%nvois
              j = xv(i)%voisin(jv)
              if(xv(j)%onprend .and. xv(i)%dist(jv) < rcn) then
                nvsel = nvsel + 1
              end if
            end do
            if(nvsel >= 17) then
              xv(i)%onprend_nb = .true.
            end if
          end if
        end if
      end do
      xv(:)%onprend = xv(:)%onprend_nb
    !   ! 🤔
!      print *, "hcp (no neighbors): ", num_sel(np, xv)
      ! do i=1,np
      !   if(xv(i)%onprend) then
      !     do jv=1,xv(i)%nvois
      !       j = xv(i)%voisin(jv)
      !       if(xv(i)%dist(jv) < rcn) then
      !         xv(j)%onprend = .true.
      !       end if
      !     end do
      !   end if
      ! end do
    end block
    nhcp = num_sel(np, xv)
    where(xv%onprend)
      xv%extracted = .true.
    end where
    call output_pop(pfil("pop_hcp.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "hcp", b%tens(1,1), 0)
    print *,"hcp: ", nhcp, nfcc + nflu + nhcp, num_extracted(np,xv)

    xv(:)%onprend = .not. xv(:)%extracted
    nbcc = num_sel(np,xv)
    where(xv%onprend)
      xv%extracted = .true.
    end where
    print *,"bcc: ", nbcc, nbcc + nfcc + nflu + nhcp, num_extracted(np,xv)
    call output_pop(pfil("pop_bcc.dat"), np, xv, b, dist)
    call output_XYZ_selection(np, xv, "bcc", b%tens(1,1), 0)
    

    open(newunit=ufc, file=pfil("phase_fractions.dat"))
    write(ufc, '(5ES16.7)') vol_frac, frac(nbcc,np), frac(nhcp,np), &
    frac(nfcc,np), frac(nflu,np)
    close(ufc)
    contains
 
  end subroutine

  integer function num_extracted(np, xv) result(res)
    integer :: np
    type(OrderData) :: xv(np)
    integer :: i
    res = 0
    do i=1,np
      if(xv(i)%extracted) res = res + 1
    end do
  end function

  integer function num_sel(np, xv) result(res)
    integer :: np
    type(OrderData) :: xv(np)
    integer :: i
    res = 0
    do i=1,np
      if(xv(i)%onprend) res = res + 1
    end do
  end function


  logical function is_not_extracted(item) result(res)
    type(OrderData), intent(in) :: item
    res = .not. item%extracted
  end function
  
  subroutine ze_selection_big_small_small_around_big(np, xv, rc_mean)
    !! Selects the set of small and big particles among Xsmall and Xbig fulfilling
    !! the criterion : r_cut_mean and n_small_neigbors_around_a_big_particle
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8) :: rc_mean
    integer :: i,j, kk, countsmall
    xv(:)%onprend = .false.
    xv(:)%isSmall = .false.
    xv(:)%isBig = .false.
    do i=1,np
      if(xv(i)%isPreBig) then
        countsmall = 0
        do kk=1, xv(i)%nvois
          j = xv(i)%voisin(kk)
          if(xv(j)%isPreSmall .and. xv(i)%dist(kk) < rc_mean) countsmall = countsmall + 1
        enddo
        if( countsmall < 13 .and. countsmall >= 10) then
          do kk=1, xv(i)%nvois
            j = xv(i)%voisin(kk)
            if(xv(j)%isPreSmall .and. xv(i)%dist(kk) < rc_mean) then
              xv(j)%onprend = .true.
              xv(j)%isSmall = .true.
            else if(xv(j)%isPreBig) then
              xv(j)%onprend = .true.
              xv(j)%isBig = .true.
            endif
          enddo
          xv(i)%onprend = .true.
          xv(i)%isBig = .true.
        endif
      endif
    enddo
  end subroutine ze_selection_big_small_small_around_big

  subroutine output_bigsmallbig(np, xv, bl, step)
    !! Outputs the positions of big and small particles of the Laves.
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: bl
    integer, intent(in) :: step
    integer :: i, nbig, nsmall, ufc=61

    nbig = bigs_count(np,xv)
    nsmall = smalls_count(np,xv)

    open(newunit=ufc, file=pfil("sel_big.xyz"))
    write(ufc, *) nbig
    write(ufc,'("L = ",ES16.7,", step = ",I10)') bl,step
    write(ufc, '(A)') "Generated by ana_order2_conf"
    do i=1,np
      if( xv(i)%isBig) then
        write(ufc,'("SIL", I0,1x, 3ES16.7)') xv(i)%fam, xv(i)%pos
      endif
    enddo
    close(ufc)

    open(newunit=ufc, file=pfil("sel_small.xyz"))
    write(ufc, *) nsmall
    write(ufc,'("L = ",ES16.7,", step = ",I10)') bl,step
    do i=1,np
      if( xv(i)%isSmall) then
        write(ufc,'("SIL", I0,1x, 3ES16.7)') xv(i)%fam, xv(i)%pos
      endif
    enddo
    close(ufc)

  end subroutine output_bigsmallbig

  real(8) function frac(n,p)
    integer, intent(in) :: n,p
    frac = real(n,8)/real(p,8)
  end function
  
  integer function smalls_count(np, xv) result(c)
    !! Returns the number of small particles
    integer :: np
    type(OrderData) :: xv(np)
    integer i
    c = 0
    do i=1, np
      if(xv(i)%isSmall) c = c + 1
    enddo
  end function

  integer function bigs_count(np, xv) result(c)
    !! Returns the number of big particles
    integer :: np
    type(OrderData) :: xv(np)
    integer i
    c = 0
    do i=1, np
      if(xv(i)%isBig) c = c + 1
    enddo
  end function
  
  subroutine calc_qs(np, xv, rc, rcin)
    !! Calculates the bond-order parameters (after a neighbor search)
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: rc, rcin
    integer :: i,j,k,m
    real(8) :: pos_u(3), d, ntvois
    do i=1,np
      xv(i)%q6v=0.d0
      xv(i)%q4v=0.d0
      ntvois = 0.0
      do k=1,xv(i)%nvois
        d = xv(i)%dist(k)
        if(d<=rc .and. d>=rcin) then
          ntvois = ntvois + 1.0
          j = xv(i)%voisin(k)
          pos_u=xv(i)%uvois(:,k)
          do m=-lOrder,lOrder
            xv(i)%q6v(m) = xv(i)%q6v(m) + SpherHarmSpecial(lOrder, m , pos_u)
          end do
          do m=-lOrder2,lOrder2
            xv(i)%q4v(m) = xv(i)%q4v(m) + SpherHarmSpecial(lOrder2, m , pos_u)
          end do
        endif
      end do
      if(ntvois > 0.0) then
        xv(i)%q6v = xv(i)%q6v / ntvois
        xv(i)%q4v = xv(i)%q4v / ntvois
      endif
      xv(i)%q6 = 0.d0
      do m=-lOrder,lOrder
        xv(i)%q6 = xv(i)%q6 + real(xv(i)%q6v(m)*conjg(xv(i)%q6v(m)), 8)
      end do
      xv(i)%q6 = sqrt(xv(i)%q6/(2*lOrder+1)) * norm_fac

      xv(i)%q4 = 0.d0
      do m=-lOrder2,lOrder2
        xv(i)%q4 = xv(i)%q4 + real(xv(i)%q4v(m)*conjg(xv(i)%q4v(m)), 8)
      end do
      xv(i)%q4 = sqrt(xv(i)%q4/(2*lOrder2+1)) * norm_fac
      
    end do
  end subroutine
  
  subroutine calc_averageq_around_nb(np, xv, rc, rcin)
    !! Calculates the neighbor-averaged bond-order parameters
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: rc, rcin
    integer :: i, m,j,k
    real(8) :: ntvois, d
    complex(8) :: barq6v(-lOrder:lOrder), barq4v(-lOrder2:lOrder2)
    do i=1,np
      barq6v = xv(i)%q6v
      barq4v = xv(i)%q4v
      ntvois = 0
      do k=1,xv(i)%nvois
        d = xv(i)%dist(k)
        if(d<=rc .and. d>=rcin) then
          ntvois = ntvois + 1
          j=xv(i)%voisin(k)
          barq6v = barq6v + xv(j)%q6v
          barq4v = barq4v + xv(j)%q4v
        endif
      end do
      barq6v = barq6v/(ntvois + 1)
      barq4v = barq4v/(ntvois + 1)
      xv(i)%barq6 = 0.d0
      do m=-lOrder,lOrder
         xv(i)%barq6 = xv(i)%barq6 + real(barq6v(m)*conjg(barq6v(m)), 8)
      end do
      xv(i)%barq6 = sqrt(xv(i)%barq6/(2*lOrder+1)) * norm_fac
      
      xv(i)%barq4 = 0.d0
      do m=-lOrder2,lOrder2
        xv(i)%barq4 = xv(i)%barq4 + real(barq4v(m)*conjg(barq4v(m)), 8)
      end do
      xv(i)%barq4 = sqrt(xv(i)%barq4/(2*lOrder+1)) * norm_fac
      
      xv(i)%barq6v = barq6v
      xv(i)%barq4v = barq4v
    end do
  end subroutine

  subroutine calc_map(np, xv, b, rc, dq)
    !! Calculates C(i) = ∑ q6(i)⋅q6(j) / |q6(i)| |q6(j)| and makes a map C(q6,q4)
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    type(Box), intent(in) :: b
    real(8), intent(in) :: rc, dq
    integer :: i, j, k, m, f
    real(8), allocatable :: C(:)
    real(8) :: d
    complex(8) :: cl
    allocate(C(np))
    do i=1,np
      C(i) = 0.0
      do k=1,xv(i)%nvois
        d = xv(i)%dist(k)
        if(d<=rc) then
          j=xv(i)%voisin(k)

          cl = cmplx(0._8, 0._8, 8)
          do m=-lOrder,lOrder
            cl = cl + xv(i)%q6v(m)*conjg(xv(j)%q6v(m))
          end do
          C(i) = C(i) + real(cl, 8)/(xv(i)%q6*xv(j)%q6)
        endif
      end do
    end do
    open(newunit=f, file=pfil("q4q6Cloud.dat"))
    write(f, '(A)') "# q4 q6 C"
    do i=1,np
      write(f, '(3ES16.7)') xv(i)%q4, xv(i)%q6, C(i)
    end do
    close(f)

    call calc_avCmap(np, xv, C, dq, pfil("q4q6Cmap.dat"))
    call calc_avCq6(np, xv, C, dq, pfil("Cvsq6.dat"))
    call calc_Chisto(np, xv, C, 0.1_8, pfil("histo_C.dat"))

    block
      integer :: neg
      neg = 0
      do i=1,np
        ! select the upper part of the map
        if(C(i)<0 .and. xv(i)%q6 > 0.5 .and. xv(i)%q4 < 0.13) then
          neg = neg + 1
          xv(i)%onprend = .true.
        else
          xv(i)%onprend = .false.
        end if
      end do
      print *, "Note : #parts(C(i)<0,q4<0.13,q6>0.5) = ", neg
      call output_XYZ_selection(np, xv, "Cinf0_topregion", b%tens(1,1), 0)
      
      neg = 0
      do i=1,np
        ! select the lower part of the map
        if(C(i)<0 .and. xv(i)%q6 < 0.3 .and. xv(i)%q4 < 0.1) then
          neg = neg + 1
          xv(i)%onprend = .true.
        else
          xv(i)%onprend = .false.
        end if
      end do
      print *, "Note : #parts(C(i)<0,q4<0.1,q6<0.3) = ", neg
      call output_XYZ_selection(np, xv, "Cinf0_bottomregion", b%tens(1,1), 0)
  
      neg = 0
      do i=1,np
        ! "crystal-like"
        if(C(i)>4) then
          neg = neg + 1
          xv(i)%onprend = .true.
        else
          xv(i)%onprend = .false.
        end if
      end do
      print *, "Note : #parts(C(i)>4) = ", neg
      call output_XYZ_selection(np, xv, "Csup4", b%tens(1,1), 0)
      
    end block
  end subroutine

  subroutine calc_avCmap(np, xv, C,  dq, filename)
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: C(np)
    real(8), intent(in) :: dq
    character(*), intent(in) :: filename

    real(8), allocatable :: avC(:,:), nsamp(:,:)
    real(8) :: q4min, q6min
    integer :: i, j, f, ih4, ih6, nh4, nh6
    
    q6min = minval(xv(:)%q6)
    nh6 = floor((maxval(xv(:)%q6)- q6min)/dq) + 1

    q4min = minval(xv(:)%q4)
    nh4 = floor((maxval(xv(:)%q4) - q4min)/dq) + 1
    allocate(avC(nh4, nh6), nsamp(nh4, nh6))
    avC = 0._8
    nsamp = 0._8
    ! print *, "nh4 = ", nh4, " nh6 = ", nh6
    do i=1, np
      ih4 = 1 + floor((xv(i)%q4-q4min)/dq)
      ih6 = 1 + floor((xv(i)%q6-q6min)/dq)

      if(ih4 >= 1 .and. ih4 <= nh4 .and. ih6 >= 1 .and. ih6 <= nh6) then
        avC(ih4, ih6) = avC(ih4, ih6) + C(i)
        nsamp(ih4, ih6) = nsamp(ih4, ih6) + 1.0_8
      end if
    end do
    where(nsamp > 0._8) avC = avC / nsamp

    open(newunit=f, file=filename)
    write(f, '(A)') "# q4 q6 C Np/N"
    do i=1,nh4
      do j=1,nh6
        write(f, '(4ES16.7)') q4min + (i-1)*dq, q6min + (j-1)*dq, avC(i,j), nsamp(i, j)/np
      end do
    end do
    close(f)
  end subroutine

  subroutine calc_avCq6(np, xv, C,  dq, filename)
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: C(np)
    real(8), intent(in) :: dq
    character(*), intent(in) :: filename

    real(8), allocatable :: avC(:), nsamp(:)
    real(8) :: q6min
    integer :: i, j, f, ih6, nh6
    
    q6min = minval(xv(:)%q6)
    nh6 = floor((maxval(xv(:)%q6)- q6min)/dq) + 1

    allocate(avC(nh6), nsamp(nh6))
    avC = 0._8
    nsamp = 0._8
    ! print *, "nh4 = ", nh4, " nh6 = ", nh6
    do i=1, np
      ih6 = 1 + floor((xv(i)%q6-q6min)/dq)

      if(ih6 >= 1 .and. ih6 <= nh6) then
        avC(ih6) = avC(ih6) + C(i)
        nsamp(ih6) = nsamp(ih6) + 1.0_8
      end if
    end do
    where(nsamp > 0._8) avC = avC / nsamp

    open(newunit=f, file=filename)
    write(f, '(A)') "# q4 q6 C Np/N"
    do j=1,nh6
      write(f, '(4ES16.7)') q6min + (j-1)*dq, avC(j), nsamp(j)/np
    end do
    close(f)
  end subroutine

  subroutine calc_Chisto(np, xv, C, dC, filename)
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: C(np)
    real(8), intent(in) :: dC
    character(*), intent(in) :: filename

    real(8), allocatable :: h(:)
    real(8) :: Cmin
    integer :: i, ih, nh, f

    Cmin = minval(C)
    nh = 1 + floor((maxval(C) - CMin)/dC)
    allocate(h(nh))
    h = 0._8
    do i=1,np
      ih = 1 + floor((C(i)-Cmin)/dC)
      if(ih >=1 .and. ih <= nh) h(ih) = h(ih) + 1
    end do
    
    open(newunit=f, file=filename)
    do i=1,nh
      write(f, '(2ES16.7)') Cmin + (i-0.5)*dC, h(i)
    end do
    close(f)
  end subroutine

  subroutine output_pop(fileName, np, pv, b, dist)
    !! Prints a distribution to a file (particle selected by field `onprend`)
    character(*), intent(in) :: fileName
    integer, intent(in) :: np
    type(Box), intent(in) :: b
    type(Distrib), intent(in) :: dist
    type(OrderData) :: pv(np)
    integer :: ufc=62,i,j
    real(8) :: popv(dist%nfam)
    popv = 0.0
    do i=1,np
      if (pv(i)%onprend) then
        j = b%parts(i)%famille
        popv(j) = popv(j) + 1
      endif
    enddo
    open(newunit=ufc, file=fileName)
    do i=1,dist%nfam
      write(ufc,'(2ES16.7)')  dist%rad(i), popv(i)
    enddo
    close(ufc)
  end subroutine
  
  subroutine search_neighbours(np, xv, b, rc, rcin)
    !! Performs a neigbor search. Result is contained in `xv`.
    integer, intent(in) :: np   !! number of particles
    type(OrderData) :: xv(np)   !! particle data array
    type(Box), intent(in) :: b  !! simulation box
    real(8), intent(in) :: rc   !! cutoff radius
    real(8), intent(in) :: rcin !! inner cutoff radius
    integer :: i,j, npris, iv, jv
    real(8) :: rc2, rcin2, tpos(3), pos(3),r2,r
    do i=1,np
      xv(i)%pos = prodm3v3(b%tens, b%parts(i)%pos)
      xv(i)%fam = b%parts(i)%famille
    end do
    rc2 = rc*rc
    rcin2 = rcin*rcin
    npris = 0
    do i=1,np
      xv(i)%nvois = 0
    enddo
    do i=1,np
      do j=i+1,np
        npris = npris + 1
        tpos = b%parts(i)%pos(1:3)-b%parts(j)%pos(1:3)
        pos=pbc(tpos, b%tens)
        r2=dot_product(pos(1:3),pos(1:3))
        if(r2<=rc2 .and. r2>rcin2) then
            xv(i)%nvois=xv(i)%nvois+1
            if(xv(i)%nvois>max_vois) then
                write(*,*) "Nombre max de voisins atteint. Désolé"
                stop 69
            end if
            xv(j)%nvois=xv(j)%nvois+1
            iv = xv(i)%nvois
            jv = xv(j)%nvois
            xv(i)%voisin(iv)=j
            xv(j)%voisin(jv)=i
            r = sqrt(r2)
            xv(i)%uvois(:,iv)= pos/r
            xv(j)%uvois(:,jv) = - xv(i)%uvois(:,iv)
            xv(i)%dist(iv) = r
            xv(j)%dist(jv) = r
        endif
      enddo
    enddo
  end subroutine
  
  logical function is_selected(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%onprend
  end function

  logical function is_not_selected(item) result(res)
    type(OrderData), intent(in) :: item
    res = .not. item%onprend
  end function
  
  
  integer function number_of_nb_sel(np, xv, i, cond, rcut) result(res)
    !! Returns number of neigbors satisfying condition `cond`.
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    integer, intent(in) :: i
    procedure(OrderCondition) :: cond
    real(8), intent(in) :: rcut
    integer :: kk,j, npris
    xv(i)%n_sel_nb = 0
    npris = 0
    do kk=1,xv(i)%nvois
      if(xv(i)%dist(kk) <= rcut) then
        j = xv(i)%voisin(kk)
        if( cond(xv(j)) ) then
          npris = npris + 1
        endif
      endif
    enddo
    res = npris
  end function
  
  ! The following procedures express the conditions
  ! for a particle to belong to a phasee

  logical function cond_bcc(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq6, bq4
    bq6 = item%barq6
    bq4 = item%barq4
    res = bq6 >= barq6_bcc_inf .and. bq6 < barq6_bcc_sup .and. bq4 >= barq4_bcc_inf .and. bq4 < barq4_bcc_sup
  end function

  logical function cond_hcp(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq6, bq4
    bq6 = item%barq6
    bq4 = item%barq4
    res = bq4 >= barq4_hcp_inf .and. bq4 < barq4_hcp_sup .and. bq6 >= barq6_hcp_inf .and. bq6 < barq6_hcp_sup
  end function

  logical function cond_fcc(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq6, bq4
    bq6 = item%barq6
    bq4 = item%barq4
    res = bq4 >= barq4_fcc_inf .and. bq4 < barq4_fcc_sup .and. bq6 >= barq6_fcc_inf .and. bq6 < barq6_fcc_sup
  end function

  logical function cond_flu(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq6, bq4
    bq6 = item%barq6
    bq4 = item%barq4
    res = bq6 >= barq6_flu_inf .and. bq6 < barq6_flu_sup .and. bq4 >= barq4_flu_inf .and. bq4 < barq4_flu_sup
  end function

  
  logical function q6_cond_1(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%q6 < q6_big_sup
  end function
  
  logical function q6_cond_2(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%q6 > q6_small_inf
  end function
  
  logical function barq6_cond_1(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%barq6 < barq6_small_sup
  end function
  
  logical function barq6_cond_2(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%barq6 > barq6_big_inf .and. item%barq6 <= barq6_big_sup
  end function

  logical function big_q6_bq6_cond(item) result(res)
    type(OrderData), intent(in) :: item
    res = q6_cond_1(item) .and. barq6_cond_2(item)
  end function
  

  logical function small_q6_bq6_cond(item) result(res)
    type(OrderData), intent(in) :: item
    res = q6_cond_2(item) .and. barq6_cond_1(item)
  end function
  

  logical function nbvois_cond_small(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%onprend .and. item%n_sel_nb < 7 .and. item%n_sel_nb >= 5
  end function

  logical function nbvois_cond_big(item) result(res)
    type(OrderData), intent(in) :: item
    res = item%onprend .and. item%n_sel_nb < 5 .and. item%n_sel_nb >= 3
  end function

  logical function cond_fcc_2(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: q4
    q4 = item%q4
    res = q4 >= q4_fcc_inf .and. q4 < q4_fcc_sup
  end function

  logical function cond_hcp_2(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: q4
    q4 = item%q4
    res = cond_hcp(item) .and. (.not. cond_fcc_2(item))
  end function

  logical function cond_fcc_alt(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq4
    bq4 = item%barq4
    res = bq4 >= fcc_alt_min_bq4
  end function

  logical function cond_hcp_alt(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq4
    bq4 = item%barq4
    res = bq4 >= hcp_alt_min_bq4
  end function

  logical function cond_flu_alt(item) result(res)
    type(OrderData), intent(in) :: item
    real(8) :: bq6
    bq6 = item%barq6
    !bq4 = item%barq4
    res = bq6 <= flu_alt_max_bq6
  end function
  
  integer function makesel(rcut, np,xv,cond, cond_nb, selected)
    !! Make a selection with a condition on particle and optionally its neighbors.
    real(8), intent(in) :: rcut !! cutoff radius
    integer, intent(in) :: np   !! number of particles
    logical, optional :: selected !! if .true., use only already selected particles
    type(OrderData) :: xv(np)   !! particle data array
    procedure(OrderCondition) :: cond !! condition procedure (logical function (type(OrderData :: x)))
    procedure(OrderCondition), optional :: cond_nb !! condition procedure for neighbors
    integer :: i,kk,j
    logical :: eff_with_nb, selected_l
    if(present(cond_nb)) then
      eff_with_nb = .true.
    else
      eff_with_nb = .false.
    endif
    if(present(selected)) then
      selected_l = selected
    else
      selected_l = .false.
    end if
    makesel = 0
    do i=1,np
      if(selected_l) then
        xv(i)%onprend = xv(i)%onprend .and. cond(xv(i))
      else
        xv(i)%onprend = cond(xv(i))
      end if
    enddo
    
    if(eff_with_nb) then
      xv(:)%onprend_nb = .false.
      do i=1,np
        if(xv(i)%onprend) then
          do kk=1,xv(i)%nvois
            if(xv(i)%dist(kk)<rcut) then
              j = xv(i)%voisin(kk)
              if(cond_nb(xv(j))) xv(j)%onprend_nb = .true.
            endif
          enddo
        endif
      enddo

    
      do i=1,np
        if(xv(i)%onprend_nb) xv(i)%onprend = .true.
      enddo
    endif
    
    do i=1,np
      if(xv(i)%onprend) makesel = makesel + 1
    enddo
  end function
  
  real(8) function calc_mean_sel_nb(np, xv)
    !! Computes the average number of **selected** neighbors of a selected particle.
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    integer :: nbcum, npris,i
    nbcum = 0
    npris = 0
    do i=1, np
      if(xv(i)%onprend) then 
        nbcum = nbcum + xv(i)%n_sel_nb
        npris = npris + 1
      endif
    enddo
    calc_mean_sel_nb = real(nbcum, 8) / real(npris, 8)
  end function

  real(8) function calc_mean_nb(np, xv)
    !! Computes the average number of neighbors of a selected particle.
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    integer :: nbcum, npris,i
    nbcum = 0
    npris = 0
    do i=1, np
      if(xv(i)%onprend) then 
        nbcum = nbcum + xv(i)%nvois
        npris = npris + 1
      endif
    enddo
    calc_mean_nb = real(nbcum, 8) / real(npris, 8)
  end function

  real(8) function calc_all_n_sel_nb(np, xv, cond, rcut)
    !! Computes the average of neighbors satisfying `cond`
    !! around a particle satisfying `cond`.
    integer, intent(in) :: np
    type(OrderData) :: xv(np)
    integer :: nv_moy, ntot
    procedure(OrderCondition) :: cond
    real(8),intent(in) :: rcut
    integer :: i
    nv_moy = 0
    ntot = 0
    do i=1,np
      if(cond(xv(i))) then
        xv(i)%n_sel_nb = number_of_nb_sel(np, xv, i, cond, rcut)
        nv_moy = nv_moy + xv(i)%n_sel_nb
        ntot = ntot + 1
      else
        xv(i)%n_sel_nb = 0
      endif
    enddo
    calc_all_n_sel_nb = real(nv_moy, 8) / real(ntot, 8)
  end function

  integer function select_crit_custom(np, xv, name, cond1, cond2, rcut, with_nb)
    !! Selects particles satisfying the two conditions `cond1` and `cond2`.
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np)  !! particle data array
    character(*), intent(in) :: name !! name of the selection
    real(8), intent(in) :: rcut  !! cutoff radius
    procedure(OrderCondition) :: cond1 !! condition 1 (logical function (type(OrderData :: item)))
    procedure(OrderCondition) :: cond2 !! condition 2 (logical function (type(OrderData :: item)))
    logical, optional :: with_nb
    logical :: eff_with_nb
    
    integer :: nnn
    real(8) :: nv_moy
    
    if(present(with_nb)) then
      eff_with_nb = with_nb
    else
      eff_with_nb = .false.
    endif
    
    write(analog,*) "Criterion :" , name
    
    nnn=0
    
 
    nnn = makesel(rcut, np, xv,cond1)
    write(analog,*) "    ", nnn, " particles pre-selected"
    nv_moy = calc_all_n_sel_nb(np, xv, is_selected, rcut) 
    write(analog,*) "    Mean number of neighbors : ", nv_moy
    
    ! Not very useful data normally.
    ! open(newunit=ufc, file=pfil("sel_" // name // "-pre.xyz"))
    ! write(ufc, *) nnn
    ! write(ufc, '(A,A,A,I0,A)') "Generated by ana_order2_conf, Criterion ", name," => ", nnn, " particles selected"
    ! do i=1,np
    !   if( xv(i)%onprend) then
    !     write(ufc,'("SPH", I0,x, 3ES16.7)') xv(i)%fam, xv(i)%pos
    !   endif
    ! enddo
    ! close(ufc)
    if(eff_with_nb) then
      nnn = makesel(rcut, np,xv,cond2, selected = .true., cond_nb=cond1)
    else
      nnn = makesel(rcut, np,xv,cond2, selected = .true.)
    end if
    write(analog,*) "    ", nnn, " particles selected"
    nv_moy = calc_all_n_sel_nb(np, xv, is_selected, rcut) 
    write(analog,*) "    Mean number of neighbors : ", nv_moy

    select_crit_custom = nnn
  end function

  subroutine output_XYZ_selection(np, xv, name, bl, step)
    !! Prints coordinates of selected particles to a file.
    integer, intent(in) :: np !! number of particles
    type(OrderData) :: xv(np)  !! particle data array
    character(*), intent(in) :: name !! name of the selection
    real(8) , intent(in) :: bl !! box length
    integer , intent(in) :: step !! step
    integer :: ufc = 199, i, nnn

    nnn = 0
    do i=1,np
      if(xv(i)%onprend) nnn = nnn + 1
    end do

    open(newunit=ufc, file=pfil("sel_" // name // ".xyz"))
    write(ufc, *) nnn
    write(ufc,'("L = ",ES16.7,", step = ",I10)') bl,step
    do i=1,np
      if( xv(i)%onprend) then
        write(ufc,'("SIL", I0, 1x, 3ES16.7)') xv(i)%fam, xv(i)%pos
      endif
    end do
    close(ufc)

  end subroutine

  function prodm3v3(m,v) result(r)
    real(8) :: r(3)
    real(8), intent(in) :: m(3,3)
    real(8), intent(in) :: v(3)
    integer :: i,j
    r = 0
    do i=1,3
      do j=1,3
        r(i) = r(i) + m(i,j)*v(j)
      end do
    end do
  end function prodm3v3
  
  function pbc(p, bt) result(rp)
    !! Vector to the first image
    !! reduced coords to real coords
    real(8) :: rp(3)
    real(8), intent(in) :: bt(3,3)
    real(8), intent(in) :: p(3)
    real(8) :: p2(3)
    integer :: i,j
    do i=1,3
      p2(i) = p(i) - int(p(i)+p(i))
    end do
    rp = 0
    do i=1,3
      do j=1,3
        rp(i) = rp(i) + bt(i,j)*p2(j)
      end do
    end do
  end function pbc
  
  subroutine output_histo(fileName, np, v, qmin, qmax, dq)
    !! Output histogram of a bond-order parameter.
    character(*) :: fileName  !! path to the output file 
    integer, intent(in) :: np !! number of particles
    real(8), intent(in) :: v(np) !! values of the bop
    real(8), intent(in) :: qmax , qmin, dq !! interval
    
    real(8), allocatable :: histo(:)
    integer :: i, nq, ig, ufc=201
    real(8) :: dqi
    
    dqi = 1.0/dq
    
    nq = floor((qmax-qmin)*dqi) + 1
    allocate(histo(nq))
    histo = 0.D0
    do i=1, np
      ig = floor((v(i)-qmin)*dqi) + 1
      if (ig>0 .and. ig <= nq) histo(ig) = histo(ig) + 1.D0
    end do
    
    histo(:) = histo/(dq*np)
    
    open(newunit=ufc, file=fileName)
    do i=1, nq
      write(ufc, '(2ES16.7)') qmin + (i-0.5)*dq, histo(i)
    enddo
    close(ufc)
    
  end subroutine
  
  subroutine write_conf_pdb_withbarq(fileName, np, xv, lengths)
    !! Prints particle positions to a file in the PDB format.
    !! Prints also mean q_6 and q_4 as a replacement of alpha,beta parameters 
    character(*) :: fileName
    integer :: np
    type(OrderData) :: xv(np)
    real(8), intent(in) :: lengths(3)
    character(3) :: ATOM_NAME
    integer :: ufc=202,i
    open(newunit=ufc,file=fileName)

    write(ufc,'(2a)') "TITLE     ","Generated by ana_order3"
    write(ufc,'(2a)') "REMARK    ", "atom charge = local bond-order (mean q6, q4)"
    write(ufc,"(a,6f8.2,3a4)") "CRYST1",10*lengths(1),10*lengths(2),10*lengths(3),90.,90.,90.,"P","1","1"
    write(ufc,"(a,8x,a)") 'MODEL', '1'
    do i=1,np
      ATOM_NAME="   "
      write(ATOM_NAME,'("S",I0)') xv(i)%fam
      write(ufc,fmt_519) "ATOM",i,"Si ",ATOM_NAME,1,10*xv(i)%pos(1),10*xv(i)%pos(2),&
            10*xv(i)%pos(3),xv(i)%barq6,xv(i)%barq4
    enddo
    close(ufc)
  end subroutine
end module bondorder

