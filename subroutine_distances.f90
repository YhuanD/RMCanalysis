    subroutine distances
! =======================

!--------------------------------------------------------------------------
!
!     This subroutine calculates the distances between two atom types
!     provided by the user in [atom1 atom2] format.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations

    implicit none

    integer, parameter :: ibonds = 14,nhistogram =40
!,ineighbour = 41
    integer :: i,j,k,l,number_cells,natoms_in_cell,iatom,jatom,ii,jj
    integer, allocatable :: histogram(:,:,:),icount(:)
    logical, allocatable :: lneighbour(:,:),lhistogram(:,:)
    character(len=2000), allocatable :: chistogram(:)
    character(len=2000) :: header
    character(len=2), allocatable :: catom(:)
    character(len=50)  :: buffer
    double precision :: dx,dy,dz,dxo,dyo,dzo,rv,xave,yave,zave
    double precision, allocatable ::xfu(:),yfu(:),zfu(:),xou(:),you(:),zou(:),&
    xfuave(:),yfuave(:),zfuave(:)
   
    number_cells = 0
    natoms_in_cell = 0

    number_cells = ncell_read(1)*ncell_read(2)*ncell_read(3)
    natoms_in_cell = natoms/number_cells
 
    allocate(chistogram(nhistogram))   
    allocate(catom(2))

    do i = 1,nhistogram
      write(buffer,'(f8.3)') i*rr_max/nhistogram
      chistogram(i) = trim(adjustl(buffer))
    end do

    allocate(xfu(natoms),yfu(natoms),zfu(natoms),xou(natoms),you(natoms),zou(natoms))
    xfu = 0.0d0
    yfu = 0.0d0
    zfu = 0.0d0

    do i = 1,natoms
      xfu(i) = xf(i)*ncell_read(1)-dble(reference_cell(i,1))
      yfu(i) = yf(i)*ncell_read(2)-dble(reference_cell(i,2))
      zfu(i) = zf(i)*ncell_read(3)-dble(reference_cell(i,3))
    end do

    allocate(xfuave(natoms_in_cell),yfuave(natoms_in_cell),zfuave(natoms_in_cell))
    allocate(icount(natoms_in_cell))
    allocate(lneighbour(natoms,natoms))
    allocate(histogram(natoms_in_cell,natoms_in_cell,nhistogram))
    allocate(lhistogram(natoms_in_cell,natoms_in_cell))

    lneighbour = .false.
    lhistogram = .false.
    histogram = 0
   
    xfuave = 0.0d0
    yfuave = 0.0d0
    zfuave = 0.0d0
    icount = 0

    do j = 1, natoms
      k = reference_number(j)
        icount(k) = icount(k) +1
        xave = xfuave(k)/icount(k)
        yave = yfuave(k)/icount(k)
        zave = zfuave(k)/icount(k)
        if (icount(k) > 0) then
                xfu(j) = xfu(j) + nint(xave-xfu(j))
                yfu(j) = yfu(j) + nint(yave-yfu(j))
                zfu(j) = zfu(j) + nint(zave-zfu(j))
        end if
        xfuave(k) = xfuave(k) + xfu(j)
        yfuave(k) = yfuave(k) + yfu(j)
        zfuave(k) = zfuave(k) + zfu(j)
    end do

    xfuave = xfuave/dble(number_cells)
    yfuave = yfuave/dble(number_cells)
    zfuave = zfuave/dble(number_cells)

!Read the required start and end atoms from the argument given when running the
! program
    bond_list = adjustl(bond_list)
    catom(1) = bond_list(1:2)
    bond_list = adjustl(bond_list(3:))
    catom(2) = bond_list(1:2)
!Now we can calculate only one pair which are two atoms each time. -Juan


!    do i = 1,natoms
!      k = reference_number(i)
!      if (element(atom_type(i))/=catom(1)) cycle
!      do jatom = 1,neighbours(i,0)
!        j = neighbours(i,jatom)
!        l = reference_number(j)
!        if (element(atom_type(j))/=catom(2)) cycle
!        dx = xfuave(k) - xfuave(l) + 1.5d0
!        dx = dx - aint(dx) - 0.5d0
!        dy = yfuave(k) - yfuave(l) + 1.5d0
!        dy = dy - aint(dy) - 0.5d0
!        dz = zfuave(k) - zfuave(l) + 1.5d0
!        dz = dz - aint(dz) - 0.5d0
!        dxo = cell(1,1)*dx/ncell_read(1) + cell(1,2)*dy/ncell_read(2) + &
!        cell(1,3)*dz/ncell_read(3)
!        dyo = cell(2,1)*dx/ncell_read(1) + cell(2,2)*dy/ncell_read(2) + &
!        cell(2,3)*dz/ncell_read(3)
!        dzo = cell(3,1)*dx/ncell_read(1) + cell(3,2)*dy/ncell_read(2) + &
!        cell(3,3)*dz/ncell_read(3)
!        rv = dsqrt(dxo**2 + dyo**2 + dzo**2)
!        if (rv > rr_max) cycle
!        lneighbour(i,j) = .true.
!      end do
!    end do
!    
!    do i = 1,natoms
!      write(41,*) i,element(atom_type(i)),count(lneighbour(i,:)),neighbours(i,0)
!    end do
!    do i = 1,natoms
!      write(42,'(i0,1x,a,100(1x,i0))') i,element(atom_type(i)),(reference_number(neighbours(i,j)),j=1,neighbours(i,0))
!    end do
!stop
!    do i = 1,natoms
!      write(6,*) lneighbour(i,:)
!    end do
! Here form the histogram
    do i = 1,natoms
      ii = reference_number(i)
      if (element(atom_type(i))/=catom(1)) cycle
      do jatom = 1,neighbours(i,0)
        j =  neighbours(i,jatom)
        if (element(atom_type(j))/=catom(2)) cycle
        jj = reference_number(j)
        dx = xf(i) - xf(j) + 1.5d0   
        dx = dx - aint(dx) - 0.5d0   
        dy = yf(i) - yf(j) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(i) - zf(j) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        rv = dsqrt(dxo**2 + dyo**2 + dzo**2)
!
        k = int(nhistogram*rv/rr_max) + 1
        if (k > nhistogram) cycle
!
        histogram(ii,jj,k) = histogram(ii,jj,k) + 1
        if (.not.lhistogram(ii,jj)) lhistogram(ii,jj) = .true.
      end do
    end do

    buffer = trim(catom(1))//'_'//trim(catom(2))//'_distances'//'.csv'
    open(ibonds,file=trim(buffer),form='formatted',status='unknown')
    header = 'r'
    do k = 1,nhistogram
      iloop: do i = 1,natoms_in_cell
        jloop: do j = 1,natoms_in_cell
          if (.not.lhistogram(i,j)) cycle jloop
          if (k==1) then
            write(buffer,'(a,i0,a,i0)') ',',i,'-',j
            header = trim(header)//trim(buffer)
          end if
          write(buffer,'(a,i6)') ',',histogram(i,j,k)
          chistogram(k) = trim(chistogram(k))//trim(adjustl(buffer))
        end do jloop
      end do iloop
    end do

    write(ibonds,*) trim(header)

    do k = 1,nhistogram
      write(ibonds,'(a)') trim(chistogram(k))
    end do
    close(ibonds)
!
    return

  end subroutine distances
