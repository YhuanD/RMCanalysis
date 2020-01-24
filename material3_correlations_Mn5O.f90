program correlations
implicit none

integer :: i,ndipoles,natoms,nneighbours,ncorrs
integer,allocatable :: ineighours(:,:),idipole(:),iatom(:)
character(len=40) :: ctemp
character(len=2) :: ctype_temp,ctype_temp2
real :: divec(:,:)


open(1,'corr.neigh',form = 'formatted','status','unknown') 
open(2,'dipole.mean',form = 'formatted','status','unknown') 
open(3,'result.corr',form = 'formatted','status','unknown')

read(2,*)
read(2,*) ctemp,ctemp,ctemp,ndipoles
do i =2,6
  read(2,*)
end do

allocate(divec(ndipoles,3))
allocate(idipole(ndipoles))

do i = 1,ndipoles
  read(2,*) idipole(i),ctype_temp,divec(i,:)
end do

read(1,*) ctemp,ctemp,ctemp,ctemp,natoms
read(1,*)
read(1,*)
read(1,*) ctemp,ctemp,ctemp,ctemp,nneighbours
 
allocate(ineighbours(natoms,nneighbours))
allocate(iatom(natoms))

read(1,*)
read(1,*)

ncorrs = nneighbours*natoms

do i = 1,natoms
  read(1,*) iatom(i),ctype_temp2,ineighbours(i,:)
end do

do i = 1,natoms
  do j = 1,ndipoles
    do k = 1,nneighbours
      if(iatom(i) == idipole(j)) 
    end do
  end do

end do


close(1)
close(2)
close(3)
end program
