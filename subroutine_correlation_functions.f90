! =================================================================================



    subroutine correlation_functions
! ====================

!----------------------------------------------------------------------------------
!
!     Calculate the average correlation functions along each direction which corresponds 
!     to .neigh file
!     Results are saved in mean.corr file.
!     Use -corr flag without input files. use rmc6f file temporarily,although it 
!     doesn't use it at all!!!! --Juan
!     sort.neigh,dipole.diag as default names.
!     Please create an empty file called mean.corr first
!----------------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations

    implicit none
    
    integer,parameter :: fcorrneigh = 20
    integer,parameter :: fdipolemean = 21
    integer,parameter :: fresult = 22   
    character(len=50) :: ctemp 
    character(len=2) :: ctype_temp
    integer :: i,ndipoles,idipole
    real,allocatable :: divec_x(:,:)
    
 !   buffer = 'corr.neigh'
 !   open(sortfile,file=trim(buffer),form='formatted',status='old')
!    buffer = 'dipole.mean'
 !   open(diagfile,file=trim(buffer),form='formatted',status='old')
 !   buffer = 'all.corr'
 !   open(corrfile,file=trim(buffer),form='formatted',status='new')

    read(fdipolemean,*) 
    read(fdipolemean,*) ctemp,ctemp,ctemp, ndipoles
    do i = 2,6
      read(fdipolemean,*) 
    end do
    !allocate(divec(ndipoles,3))
    do i = 1,ndipoles
 !     read(fdipolemean,*) idipole,ctype_temp,divec(i,:)
    end do
    
    close(fcorrneigh)
    close(fdipolemean)
    close(fresult)
    return
  end subroutine correlation_functions
