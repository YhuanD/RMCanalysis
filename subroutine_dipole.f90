! =================================================================================
    subroutine dipole
! ====================

!----------------------------------------------------------------------------------
!
!     This subroutine gives all the dipole vectors and average and deviations information. 
!     This subroutine use .neigh file from -sort_neigh, and refined configuration file.
!     Output file is .dipole file
! The file dipole.neigh must exist
! e.g. ./rmcanalysis -dipole saved_configuration_15K.rmc6f
!----------------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use coordinates, only: vectors
    
    implicit none

    integer, parameter :: neighfile = 15
    integer, parameter :: dipolefile = 16
    integer :: i,j,ndipoles
    integer :: ierror
    integer,allocatable :: iatom1(:)
    integer :: nneigh
    integer, allocatable :: iatom2(:)
    character(len=50)  :: buffer
    character(len=100) :: ctemp
    character(len=2) :: catom1,catom2
    double precision :: dx,dy,dz,dxo,dyo,dzo
    double precision, allocatable :: netdipole(:,:),dimean(:),dimean_sq(:),devi(:)
    double precision, allocatable :: vec(:)
    
    buffer = 'dipole.neigh'
    open(neighfile,file=trim(buffer),form='formatted',status='unknown')
    buffer = 'dipole.mean'
    open(dipolefile,file=trim(buffer),form='formatted',status='unknown')
    
    allocate(dimean(3))
    allocate(dimean_sq(3))
    allocate(devi(3))

    dimean = 0.0d0
    dimean_sq = 0.0d0
    devi = 0.0d0
    nneigh = 0
    ctemp = ''
    ndipoles = 0
    
    read(neighfile,*) ctemp
    read(neighfile,*) ctemp,ctemp,ctemp,ctemp,ndipoles
    read(neighfile,*) ctemp
    read(neighfile,*) ctemp,ctemp,ctemp,ctemp,nneigh
    read(neighfile,*) ctemp
    read(neighfile,*) ctemp
    
    allocate(iatom2(nneigh))
    allocate(iatom1(ndipoles))
    allocate(vec(nneigh))
    allocate(netdipole(ndipoles,3))
    
    netdipole = 0.0d0
    iatom1 = 0
    iatom2 = 0
    vec = 0.0d0
    
    do i = 1,ndipoles
      read(neighfile,*) iatom1(i),ctemp,iatom2(:)
      catom1 = element(atom_type(iatom1(i)))
      catom2 = element(atom_type(iatom2(1)))     
      do j = 1,nneigh
        dx = xf(iatom2(j)) - xf(iatom1(i)) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(iatom2(j)) - yf(iatom1(i)) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(iatom2(j)) - zf(iatom1(i)) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz    
        netdipole(i,1) = netdipole(i,1) + dxo
        netdipole(i,2) = netdipole(i,2) + dyo
        netdipole(i,3) = netdipole(i,3) + dzo
      end do 
    end do

! calcualte average and deviations
    do i = 1,ndipoles
      dimean(1) = dimean(1) + netdipole(i,1)
      dimean(2) = dimean(2) + netdipole(i,2)
      dimean(3) = dimean(3) + netdipole(i,3)
      dimean_sq(1) = dimean_sq(1) + netdipole(i,1)*netdipole(i,1)
      dimean_sq(2) = dimean_sq(2) + netdipole(i,2)*netdipole(i,2)
      dimean_sq(3) = dimean_sq(3) + netdipole(i,3)*netdipole(i,3)
    end do

    dimean(1) = dimean(1)/ndipoles
    dimean(2) = dimean(2)/ndipoles
    dimean(3) = dimean(3)/ndipoles

    dimean_sq(1) = dimean_sq(1)/ndipoles
    dimean_sq(2) = dimean_sq(2)/ndipoles
    dimean_sq(3) = dimean_sq(3)/ndipoles
    
    devi(1) = dimean_sq(1) - dimean(1)*dimean(1)
    devi(2) = dimean_sq(2) - dimean(2)*dimean(2)
    devi(3) = dimean_sq(3) - dimean(3)*dimean(3)

    catom1 = trim(adjustl(catom1))
    catom2 = trim(adjustl(catom2))
    
! write out dipole file 
    write(dipolefile,'(a,a,a,a)') 'Type of dipole: ',catom1,'-',catom2
    write(dipolefile, '(a,i0)') 'Number of dipoles: ', ndipoles
    write(dipolefile,'(a,a,a,i0)') 'Number of neighbours: ',catom1,': ',nneigh
    write(dipolefile,'(a,3f10.4)') 'Mean_dipole = ', dimean(:)
    write(dipolefile,'(a)') 'Deviation = <miu(r)* miu(r)> - <miu(r)>* <miu(r)>'
    write(dipolefile,'(a,3f10.4)') 'devi = ',devi(:)
    do i = 1,ndipoles
      write(dipolefile,'(i0,2x,a,2x,3f12.6)') iatom1(i),catom1,netdipole(i,:)
    end do
    
!end write out dipole file    
    close(dipolefile)
    close(neighfile)
    return
  end subroutine dipole
