program sort_Y1Y2
implicit none
!./sort_dipole_Y1Y2 
! input: dipole_Y8O_Y1Y2.mean, dipole_initial.rmc6f
! output: dipole_Y1.sort dipole_Y2.sort
integer :: i,j,ndipole,natoms,natoms_Y,iatom_temp,natom_MnO,idipole
integer,allocatable :: iatom(:)
character(len=50) :: ctemp
character(len=2) :: catom
character(len=3) :: csite_temp
character(len=3),allocatable :: csite(:)
real :: vec_x,vec_y,vec_z

open(1,file='dipole_Y8O_Y1Y2.mean',form='formatted',status='unknown')
open(2,file='dipole_initial.rmc6f',form='formatted',status='unknown')
open(3,file='dipole_Y1.sort',form='formatted',status='unknown')
open(4,file='dipole_Y2.sort',form='formatted',status='unknown')

do i = 1,5
  read(2,*)
end do

read(2,*) ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,natoms_Y

allocate(iatom(natoms_Y))
allocate(csite(natoms_Y))

do i = 7,10
  read(2,*) 
end do

read(2,*) ctemp,ctemp,ctemp,natoms

do i = 12,19
  read(2,*)
end do

natom_MnO = natoms - natoms_Y
! loop over Mn and O
do i = 1,natom_MnO
  read(2,*)
end do
do i = 1,natoms_Y
  read(2,*) iatom_temp,catom,csite_temp
  catom = trim(adjustl(catom))
  if (catom == 'Y') then
     iatom(i) = iatom_temp
     csite(i) = csite_temp
  end if
end do

read(1,*)
read(1,*) ctemp,ctemp,ctemp,ndipole
do i = 3,6
  read(1,*)
end do

do i = 1,ndipole
  read(1,*) idipole,ctemp,vec_x,vec_y,vec_z
  do j = 1,natoms_Y
    if(iatom(j)==idipole) then
      if(csite(j)=='[1]') write(3,'(3f12.6)') vec_x,vec_y,vec_z
      if(csite(j)=='[2]') write(4,'(3f12.6)') vec_x,vec_y,vec_z
      exit
    end if 
  end do
end do


close(1)
close(2)
close(3)
close(4)

end program
