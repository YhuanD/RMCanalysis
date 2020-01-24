    subroutine structure_analysis(icif,itf)
! ==========================================

!----------------------------------------------------------------------------------
!
!     This subroutine is for analysis of the structure e.g calculate the average
!     bond length in range min to max for specific two atoms, calculate the
!     deviations from the average bond length, etc.
!
!----------------------------------------------------------------------------------
    use arguments
    use structure_data
    use cif_stuff 
    use useful_stuff
    use mod_commons_constants, only: pi
    use common_arguments, only: ldiag
    implicit none

    integer, parameter :: avevecfile = 23
    double precision, allocatable :: xfu(:),yfu(:),zfu(:),xou(:),you(:),zou(:)
    double precision, allocatable :: xfuave(:),yfuave(:),zfuave(:),asmall,bsmall,&
                                     csmall,alphasmall,betasmall,gammasmall
    double precision, allocatable :: xouave(:),youave(:),zouave(:),dxo(:),dyo(:),dzo(:),&
                                     dx(:),dy(:),dz(:)
    double precision, allocatable :: beta_DW(:,:,:),b_DWsmall(:,:,:),b_DW_buffer(:,:,:)
    double precision :: xave,yave,zave,F(3,3),b1,b2,buffera(3), beta3
    integer,allocatable :: icount(:)
    integer :: i,j,k,number_cells,natoms_in_cell,itf,icif
    character(len=2),allocatable :: el(:)

    open(itf,file='vertifluc.csv',form='formatted',status='unknown')
  if (ldiag) then
    write(main_output,*)
    write(main_output,*) 'Output from structure_analysis'
    write(main_output,*) ' Number of sites: ',natoms
    write(main_output,*)
    write(main_output,'(a)') ' === Configuration in fractional coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),&
              xf(i),yf(i),zf(i)
    end do
    write(main_output,*)
    write(main_output,'(a)') ' === Configuration in orthogonal coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),&
             xo(i),yo(i),zo(i)
    end do
    write(main_output,*)
  end if


    number_cells = 0
    natoms_in_cell = 0

    allocate(xfu(natoms),yfu(natoms),zfu(natoms),xou(natoms),you(natoms),zou(natoms))
    xfu = 0.0d0
    yfu = 0.0d0
    zfu = 0.0d0

    do i = 1,natoms
      xfu(i) = xf(i)*ncell_read(1)-dble(reference_cell(i,1))
      xfu(i) = xfu(i) - int(xfu(i))
      yfu(i) = yf(i)*ncell_read(2)-dble(reference_cell(i,2))
      yfu(i) = yfu(i) - int(yfu(i))
      zfu(i) = zf(i)*ncell_read(3)-dble(reference_cell(i,3))
      zfu(i) = zfu(i) - int(zfu(i))
    end do

    asmall = a/ncell_read(1)
    bsmall = b/ncell_read(2)
    csmall = c/ncell_read(3)

    alphasmall = alpha
    betasmall = beta
    gammasmall = gamma
    
    allocate(el(natoms))
    do i = 1,natoms
      el(i) = element(atom_type(i))
    end do

!write a CIF file for the reduced coordinates
    call write_cif(icif,'reduced.cif',asmall,bsmall,csmall,alphasmall,betasmall,&
                 gammasmall,xfu,yfu,zfu,el)
    deallocate(el)

    ! So far so good, or at least the CIF files look OK.

    number_cells = ncell_read(1)*ncell_read(2)*ncell_read(3)
    natoms_in_cell = natoms/number_cells
    
    allocate(xfuave(natoms_in_cell),yfuave(natoms_in_cell),zfuave(natoms_in_cell))
    allocate(xouave(natoms_in_cell),youave(natoms_in_cell),zouave(natoms_in_cell))
    allocate(beta_DW(natoms_in_cell,3,3),b_DWsmall(natoms_in_cell,3,3),b_DW_buffer(natoms_in_cell,3,3))
    allocate(el(natoms_in_cell))
    allocate(dxo(natoms),dyo(natoms),dzo(natoms))
    allocate(dx(natoms),dy(natoms),dz(natoms))
    allocate(icount(natoms_in_cell))

    xfuave = 0.0d0
    yfuave = 0.0d0
    zfuave = 0.0d0
    icount = 0
    beta_DW = 0.0d0
    b_DWsmall = 0.0d0
    F = 0.0d0

    do j =1, natoms
      k =reference_number(j)
        icount(k) = icount(k) +1
        xave = xfuave(k)/icount(k)
        yave = yfuave(k)/icount(k)
        zave = zfuave(k)/icount(k)
        if(icount(k) == 1) then
          el(k) = element(atom_type(j))
        end if
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

!    call write_cif(icif,'averaged.cif',asmall,bsmall,csmall,alphasmall,betasmall,&
!                 gammasmall,xfuave,yfuave,zfuave,el)

    ! So far so good: have correctly averaged the positions of each atom.
 
    natoms = natoms_in_cell
    natoms = natoms * number_cells
    do i = 1,3
        cell(:,i) = cell(:,i)/ncell_read(i)
    end do
! Convert the mean fractional coordinates into cartesian coordinates
    xouave = cell(1,1)*xfuave + cell(1,2)*yfuave + cell(1,3)*zfuave
    youave = cell(2,1)*xfuave + cell(2,2)*yfuave + cell(2,3)*zfuave
    zouave = cell(3,1)*xfuave + cell(3,2)*yfuave + cell(3,3)*zfuave
    xo = cell(1,1)*xfu + cell(1,2)*yfu + cell(1,3)*zfu 
    yo = cell(2,1)*xfu + cell(2,2)*yfu + cell(2,3)*zfu
    zo = cell(3,1)*xfu + cell(3,2)*yfu + cell(3,3)*zfu

!    ! Write an xyz file
!    open(84, file='orthogonal.xyz')
!    write(84, '(i10/)') natoms_in_cell
!    do k = 1, size(xouave)
!       write(84, '(a3, 1x, 3f10.5)') el(k), xouave(k), youave(k), zouave(k)
!    end do
!    close(84)

    ! So far so good: correctly orthogonalised.

! Calculate the fluctuations of positions

    do i = 1,natoms
      k =reference_number(i)
   
      dx(i) = xfu(i) - xfuave(k)
      dx(i) = dx(i) - nint(dx(i))
      dy(i) = yfu(i) - yfuave(k)
      dy(i) = dy(i) - nint(dy(i)) 
      dz(i) = zfu(i) - zfuave(k)
      dz(i) = dz(i) - nint(dz(i)) 

      dxo(i) = cell(1,1)*dx(i) + cell(1,2)*dy(i) + cell(1,3)*dz(i)
      dyo(i) = cell(2,1)*dx(i) + cell(2,2)*dy(i) + cell(2,3)*dz(i)
      dzo(i) = cell(3,1)*dx(i) + cell(3,2)*dy(i) + cell(3,3)*dz(i)
 
      beta_DW(k,1,1) = beta_DW(k,1,1) + dxo(i)**2
      beta_DW(k,2,2) = beta_DW(k,2,2) + dyo(i)**2
      beta_DW(k,3,3) = beta_DW(k,3,3) + dzo(i)**2
      beta_DW(k,1,2) = beta_DW(k,1,2) + dxo(i)*dyo(i)
      beta_DW(k,1,3) = beta_DW(k,1,3) + dxo(i)*dzo(i)
      beta_DW(k,2,3) = beta_DW(k,2,3) + dyo(i)*dzo(i)
      beta_DW(k,2,1) = beta_DW(k,2,1) + dxo(i)*dyo(i)
      beta_DW(k,3,1) = beta_DW(k,3,1) + dxo(i)*dzo(i)
      beta_DW(k,3,2) = beta_DW(k,3,2) + dyo(i)*dzo(i)
      
    end do

    beta_DW = beta_DW/dble(number_cells)

!    do i = 1,natoms_in_cell
!      do j = 1,3
!        do k = 1,3
!          write(itf,*) i,j,k,beta_DW(i,j,k)
!        end do
!      end do
!    end do
  
    call write_cif(icif,'mean.cif',asmall,bsmall,csmall,alphasmall,betasmall,&
                 gammasmall,xfuave,yfuave,zfuave,el,beta_DW)

    ! write in average with fluctuations in vector format. Used specifically for ScF3
    open(avevecfile,file='average.vec',form='formatted',status='unknown')
    
    write(avevecfile,'(3a)') 'Atom number in rmc6f file, ','Atom type, ',&
    'Vector from average to real' 
    write(avevecfile,*) natoms
    do i = 1,natoms
      write(avevecfile,'(a,2x,i0,3f12.6)') element(atom_type(i)),i,dxo(i),dyo(i),dzo(i)
    end do
    
    close(avevecfile)
    return
  end subroutine structure_analysis
