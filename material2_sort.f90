program sort_O
implicit none

character(len= 20) :: ctemp
integer :: i,natoms,iatom,atom_number
character(len=2) :: catomtype
real,allocatable :: frac(:)

natoms = 6600
catomtype = 'O'
allocate(frac(3))

! input .rmc6f file format for rmc
open(1,file='standard.rmc6f',form='formatted',status='unknown')
open(2,file='O3.sort',form='formatted',status='unknown')
open(10,file='O4.sort',form='formatted',status='unknown')
open(3,file='O5.sort',form='formatted',status='unknown')
open(11,file='O6.sort',form='formatted',status='unknown')

do i = 1,6619
  read(1,*) ctemp
end do
do i = 6620,13219
  read(1,*) iatom,catomtype,ctemp,frac(:),atom_number
  if(atom_number == 3) then
    write(2,'(i0,2x,a,2x,3f10.6,2x,i0)') iatom,catomtype,frac(:),atom_number
  else if(atom_number == 4) then
    write(10,'(i0,2x,a,2x,3f10.6,2x,i0)') iatom,catomtype,frac(:),atom_number
  else if (atom_number == 5) then
    write(3,'(i0,2x,a,2x,3f10.6,2x,i0)') iatom,catomtype,frac(:),atom_number
  else if (atom_number == 6) then
    write(11,'(i0,2x,a,2x,3f10.6,2x,i0)') iatom,catomtype,frac(:),atom_number
  else 
    write(6,*) 'wrong atom number'
    stop
  end if
end do

close(1)
close(2)
close(3)
end program
