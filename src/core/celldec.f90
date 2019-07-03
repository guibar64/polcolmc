
module celldec
  !! Implements a cell decomposition.
  !!
  !!## Implementation notes
  !! 
  !! Each cell refers to a doubly linked list of `Particule` nodes.
  !! For each simulation box a standard array of particles is defined.
  !! The saner way would have been to use *ints* to refer to particles
  !! in the cell list, but here *pointers* are used. The reason is that
  !! the latter seems faster for energy calculations (~ 7% speedup).
  !! But it means that you **cannot reallocate the array of particles** 
  !! without reconstructing the whole lists. The same is true if you
  !! want to "delete" a particle and conserve the order of the other
  !! elements (though you can swap the element with the last element
  !! of the array and then delete "easily" the last element).
  implicit none
  
  ! `Particle derivede type must be defined here or an used module because it is tied to
  ! the cell decomposition.
  type  Particule
    !! Particle structure
    real(8),dimension(3) :: pos !! position
    integer :: numero   !! number (id)
    integer :: famille  !! family
    integer :: cellule  !! cell number
    real(8) :: rayon    !! radius
    real(8) :: ech      !! (effective) charge
    type(Particule),pointer :: suivant !! points to the next particle in the cell
    type(Particule),pointer :: precedent !! points to the previous particle in the cell
   end type Particule

  type ClistCell
    !! Structure representing a cell
    type(Particule), pointer :: hoc  !! head of the list
    integer :: coord(3)            !! grid coordinates of the cells
    integer :: nvois
    integer :: voisin(27)   !! neighboring cells
  end type ClistCell

  type CellDecomp
    !! Cell decomposition
    type(ClistCell), allocatable :: cells(:) !! All cells
    integer :: icell(3)            !! number of cells along each direction
    integer ::  ncell              !! number of cells
    real(8) :: rcell(3)            !! dimensions of a cell
  end type CellDecomp

  private :: indice_cell

contains

  subroutine part_transfer_props(source, dest)
    !! transfer particle properties from `source` to `dest`
    type(Particule), intent(in) :: source
    type(Particule), intent(out) :: dest
    dest%famille = source%famille
    dest%rayon = source%rayon
    dest%ech = source%ech
    dest%numero = source%numero
  end subroutine part_transfer_props

  integer function indice_cell(cellc,icell)
    !! Returns the cell index from the grid coordinates
    implicit none
    integer :: cellc(3),icell(3)

    indice_cell = 1 + cellc(1)*icell(2)*icell(3) + icell(3)*cellc(2) + cellc(3)

  end function indice_cell

  subroutine celldec_cree(clist,n,x,r_cell)
    !! constructs `clist` from a particle array and takes the dimensions of a cell as input (`r_cell`)
    implicit none
    integer,intent(in) :: n
    type(CellDecomp) :: clist
    real(8),intent(in) :: r_cell(3)
    type(Particule),intent(inout),target :: x(*)

    integer :: i,j,k,c,idump,counter=0
    integer :: cell_co(3),icell(3),ncell
    real(8) :: rcell(3)
    counter=counter+1
      
    icell = max(floor(1.d0/r_cell),3) ! 'reduced' box lengths = 1
    ncell = icell(1)*icell(2)*icell(3)

    if(allocated(clist%cells)) then
      if(clist%ncell<ncell) then
        deallocate(clist%cells)
        allocate(clist%cells(ncell))
      end if
    else
      allocate(clist%cells(ncell))
    end if

    rcell = 1.d0/icell
    clist%icell = icell
    clist%ncell = ncell        
    clist%rcell = rcell

    ! Here we puts all neighboring cell indexes in the array `voisin`
    ! to avoid recomputing them later.
    do c=1,ncell
      k=c-1
      clist%cells(c)%coord(1)=k/(icell(2)*icell(3))
      idump = k-icell(2)*icell(3)*clist%cells(c)%coord(1)
      clist%cells(c)%coord(2)=idump/icell(3)
      idump = idump - icell(3)*clist%cells(c)%coord(2)
      clist%cells(c)%coord(3) = idump

      idump=0

      clist%cells(c)%nvois=27
      do i=-1,1
        do j=-1,1
          do k=-1,1
            idump=idump+1
            cell_co(1)=modulo(i+clist%cells(c)%coord(1),icell(1))
            cell_co(2)=modulo(j+clist%cells(c)%coord(2),icell(2))
            cell_co(3)=modulo(k+clist%cells(c)%coord(3),icell(3))
            clist%cells(c)%voisin(idump) = indice_cell(cell_co,icell)
          end do
        end do
      end do
    end do

    call celldec_update_all(clist, n, x)

  end subroutine celldec_cree

  subroutine celldec_update_all(clist, n, x)
    !! updates `clist` from a particle array.
    implicit none
    integer,intent(in) :: n   !! number of particles
    type(CellDecomp) :: clist
    type(Particule),intent(inout),target :: x(n) !! particles
   
    integer :: i, c, cell, cell_co(3)
    type(Particule), pointer :: p1, p2
    ! Lists construction
    do c=1,clist%ncell 
      clist%cells(c)%hoc => null()     
    end do

    do i=1,n
      p1 => x(i)
      cell_co=floor(p1%pos/clist%rcell)
      cell=indice_cell(cell_co,clist%icell)
      p1%cellule=cell
      p1%suivant => clist%cells(cell)%hoc
      clist%cells(cell)%hoc => p1
    enddo

    ! reverse
    do c=1,clist%ncell
      p1 => clist%cells(c)%hoc
      p2 => null()
      do while(associated(p1))
        p1%precedent => p2
        p2 => p1
        p1 => p1%suivant
      end do
    enddo
  end subroutine celldec_update_all

  subroutine celldec_update(part,clist)
    !! update `clist` after a particle translation
    implicit none
    type(Particule),target :: part
    type(CellDecomp) :: clist
    integer :: cell,cell_co(3)
    type(Particule),pointer :: p1,p2

    cell_co=floor(part%pos/clist%rcell)
    cell=indice_cell(cell_co,clist%icell)

    if(cell/=part%cellule) then

      ! Remove from previous cell
      p1 => part%precedent
      p2 => part%suivant

      if(associated(p1)) then
        p1%suivant => p2
      else
        clist%cells(part%cellule)%hoc => p2
      end if
      if(associated(p2)) then
        p2%precedent => p1
      end if

      ! put into new cell (head of list)

      part%cellule = cell
      p1 => clist%cells(cell)%hoc 
      part%suivant => p1

      if(associated(p1)) p1%precedent => part
      part%precedent => null()
      clist%cells(cell)%hoc => part

    end if

  end subroutine celldec_update

  subroutine celldec_add(part,clist)
    !! adds a particle to `clist`
    implicit none
    type(Particule),target :: part
    type(CellDecomp) :: clist
    integer :: cell,cell_co(3)
    type(Particule),pointer :: p1

    cell_co=floor(part%pos/clist%rcell)
    cell=indice_cell(cell_co,clist%icell)

    ! put into new cell (head of list)

    part%cellule = cell
    p1 => clist%cells(cell)%hoc 
    part%suivant => p1

    if(associated(p1)) p1%precedent => part
    part%precedent => null()
    clist%cells(cell)%hoc => part
    
  end subroutine celldec_add

  subroutine celldec_update_b2b(part,clist_orig,clist_dest)
    !! updates `clist_orig` and `clist_dest` after a tranfer of a particle from a box to another one.
    implicit none
    type(Particule),target :: part
    type(CellDecomp) :: clist_orig,clist_dest
    integer :: cell,cell_co(3),counter=0
    type(Particule),pointer :: p1,p2
    counter=counter+1
    cell_co=floor(part%pos/clist_dest%rcell)
    cell=indice_cell(cell_co,clist_dest%icell)


    ! Remove from previous cell
    p1 => part%precedent
    p2 => part%suivant
    if(associated(p1)) then
      p1%suivant => p2
    else
      clist_orig%cells(part%cellule)%hoc => p2
    end if
    if(associated(p2)) then
      p2%precedent => p1
    end if

    ! put into new cell (head of list)

    part%cellule = cell
    p1 => clist_dest%cells(cell)%hoc 
    part%suivant => p1
    if(associated(p1)) p1%precedent => part
    part%precedent => null()
    clist_dest%cells(cell)%hoc => part

  end subroutine celldec_update_b2b

  subroutine celldec_echange(part1,part2,clist)
    !! update `clist` after a swap between two particles
    implicit none
    type(Particule),target :: part1,part2
    type(CellDecomp) :: clist
    integer :: cell1,cell2
    type(Particule),pointer :: p1

    cell1=part1%cellule
    cell2=part2%cellule
   
    if(cell1/=cell2) then
       p1 => part1%suivant
       
       part1%suivant => part2%suivant
       if(associated(part1%suivant)) part1%suivant%precedent => part1

       part2%suivant => p1
       if(associated(p1)) p1%precedent => part2
       

       p1 => part1%precedent
       
       part1%precedent => part2%precedent
       if(associated(part2%precedent)) then
          part1%precedent%suivant => part1
       else
          clist%cells(cell2)%hoc => part1
       end if

       part2%precedent => p1
       if(associated(p1)) then
          p1%suivant => part2
       else
          clist%cells(cell1)%hoc => part2
       end if

       part1%cellule = cell2
       part2%cellule = cell1
    end if
    
  end subroutine celldec_echange

  integer function get_cell(clist,pos)
    !! position to the index of the cell it belongs
    implicit none
    type(CellDecomp),intent(in) :: clist
    real(8),intent(in) :: pos(3)  !! (reduced) position
    integer :: cell_co(3) 
    cell_co=floor(pos/clist%rcell)
    get_cell=indice_cell(cell_co,clist%icell)
  end function get_cell

  subroutine print_clist(fu,clist,n,x)
    !! prints cell decomposition
    implicit none
    type(CellDecomp) :: clist
    integer,intent(in) :: n   !! number of particles
    integer, intent(in) :: fu !! file descriptor 
    type(Particule),intent(In), target :: x(n) !! particles
    integer :: cell, nb
    type(Particule), pointer :: cpar
    write(fu,*) "Liste cellules"
    nb = 0
    do cell=1,clist%ncell
       write(fu,'(A,I4,A,I3,A1,I3,A1,I3,A)',advance='no') "cell(",cell,")(",clist%cells(cell)%coord(1),",",&
            clist%cells(cell)%coord(2),",",clist%cells(cell)%coord(3),"): "
       cpar => clist%cells(cell)%hoc
       do while(associated(cpar))
          nb = nb + 1
          write(fu,'(I5)',advance='no') cpar%numero
          cpar => cpar%suivant
       end do
       write(fu,'(A,I0)',advance='no') " , nvois= ",clist%cells(cell)%nvois
       write(fu,*)
    end do
    write(fu,*) " Nb parts : ", nb
  end subroutine

  subroutine celldec_numofpart(clist,numv,ncell)
    !! Returns the number of particles in the cells (`numv`) and the number of cells (`ncell`)
    implicit none
    type(CellDecomp) :: clist
    integer,intent(inout), allocatable :: numv(:)
    integer,intent(out) :: ncell
    type(Particule), pointer :: cpar
    integer :: cell
    if(allocated(numv)) then
      if(ubound(numv,1) < clist%ncell) then
        deallocate(numv)
        allocate(numv(ncell))
      end if
    else
      allocate(numv(ncell))
    end if

    ncell = clist%ncell
    do cell=1,clist%ncell
       cpar => clist%cells(cell)%hoc
       numv(cell)=0
       do while(associated(cpar))
          numv(cell) = numv(cell) + 1
          cpar => cpar%suivant
       end do
    end do
    
  end subroutine celldec_numofpart

end module celldec
