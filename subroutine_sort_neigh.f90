! =================================================================================
    subroutine sort_neighbours
! ====================

!----------------------------------------------------------------------------------
!
!     This subroutine sorts the order of the neighbours of each atom given in the 
!     .rmc6f file for specific atom type provided by the user.
!     Use -sort_neigh to give the atom type, and -neigh_range to give the maximum
!     range you want to calculate and -num_neigh to give the number of neighbours 
!     of each atom.
!     e.g. -sort_neigh [Fe,O] -neigh_range [1.5,3] -num_neigh 6 *.rmc6f
!	  Run this on initial config generated from data2config to get correct 
!     neighbour information. dipole subroutine will use the .neigh file from this subroutine.
!----------------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use coordinates, only: vectors

    implicit none
    
    integer,parameter :: neighfile = 88
    integer :: i,j,k,l,ii,jj,kk,ll,m,jatom,count_temp,kmax,kmin,n
    integer :: nfake,nsort,ncount,temp_test,ierror
    character(len=2) :: catom_one,catom_two
    character(len=30) :: buffer,cneigh_temp
    logical,allocatable :: lneighbour(:,:)
    integer,allocatable :: neigh_list(:,:),new_list(:,:),count_neigh(:),ntemp(:)
    double precision :: dx,dy,dz,mean_bond
    double precision,allocatable :: dxo(:,:),dyo(:,:),dzo(:,:),atemp(:),rv(:,:)
    real :: neigh_rangemin,neigh_rangemax
     
    csort_neigh = adjustl(csort_neigh)
    n=index(csort_neigh,',')
    cneigh_temp = csort_neigh(1:n-1)
    ierror = 0
    read(cneigh_temp,*,iostat=ierror) catom_one
    cneigh_temp = adjustl(csort_neigh(n+1:))
    ierror = 0
    read(cneigh_temp,*,iostat=ierror) catom_two
    
    neigh_range= adjustl(neigh_range)
    n=index(neigh_range,',')
    cneigh_temp = neigh_range(1:n-1)
    ierror = 0
    read(cneigh_temp,*,iostat=ierror) neigh_rangemin
    cneigh_temp= adjustl(neigh_range(n+1:))
    ierror = 0
    read(cneigh_temp,*,iostat=ierror) neigh_rangemax
    
    allocate(neigh_list(natoms,nnum_neigh))
    allocate(new_list(natoms,nnum_neigh))
    allocate(count_neigh(natoms))
    allocate(dxo(natoms,natoms))
    allocate(dyo(natoms,natoms))
    allocate(dzo(natoms,natoms))
    allocate(ntemp(1))
    allocate(atemp(nnum_neigh))
    allocate(lneighbour(natoms,natoms))
    allocate(rv(natoms,natoms))
    
    lneighbour = .false.
    neigh_list = 0
    new_list = 0
    nfake = 0
    nsort = 0
    ncount = 0
    count_neigh = 0
    ntemp = 0
    mean_bond = 0
    
    buffer = trim(catom_one)//'_'//trim(catom_two)//'_sort'//'.neigh'
    open(neighfile,file=trim(buffer),form='formatted',status='unknown')

    do i = 1,natoms
      ii = reference_number(i)
      if (element(atom_type(i))/=catom_one) cycle
      do jatom = 1,neighbours(i,0)
        j =  neighbours(i,jatom)
        if (element(atom_type(j))/=catom_two) cycle
        jj = reference_number(j)
        dx = xf(j) - xf(i) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(j) - yf(i) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(j) - zf(i) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo(i,j) = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo(i,j) = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo(i,j) = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        rv(i,j) = dsqrt(dxo(i,j)**2 + dyo(i,j)**2 + dzo(i,j)**2)
        kmax = int(100*rv(i,j)/neigh_rangemax) + 1
        kmin = int(100*rv(i,j)/neigh_rangemin) + 1
        if (kmax > 100 .or. kmin < 100) cycle
        count_neigh(i) = count_neigh(i) + 1
        lneighbour(i,j) = .true.
      end do
      if (count_neigh(i) /= nnum_neigh) then
        nfake = nfake+1
        ncount = ncount +1
      else 
        nsort = nsort + 1
        ncount = ncount +1
      end if
    end do  

!get neighbours of i atom saved in neigh_list    
    count_temp = 0
    do i = 1,natoms 
      if (element(atom_type(i))/=catom_one) cycle
      if(count_neigh(i) /= nnum_neigh) cycle
      k = 0
      do jatom = 1,neighbours(i,0)
        j =  neighbours(i,jatom)
        if(element(atom_type(j))/=catom_two) cycle
        if(lneighbour(i,j)) then         
          k = k + 1
          neigh_list(i,k) = j
          mean_bond = mean_bond + rv(i,j)
          count_temp = count_temp + 1
        end if
      end do
    end do
    mean_bond = mean_bond/count_temp
    
    write(neighfile,*) 'Total number of atoms: ',ncount
    write(neighfile,*) 'Number of atoms sorted: ',nsort
    write(neighfile,*) 'Number of fake atoms: ',nfake
    write(neighfile,*) 'Number of neighbours set: ',nnum_neigh
! write out average catom_one - catom_two bond length, take acount of all neighbours within 
! range and doesnt check the number of neighbours.
    write(neighfile,*) 'Average bond length of ',catom_one,'-',catom_two,':', mean_bond
    write(neighfile,*) 'Number of bonds taken account of in average: ', count_temp
    
    do m=1,natoms
      if(element(atom_type(m))/=catom_one) cycle
      if(count_neigh(m) /= nnum_neigh) cycle
      write(neighfile,*) m,element(atom_type(m)),neigh_list(m,:)
    end do
	  
    close(neighfile)
    
    return
  end subroutine sort_neighbours
