program coscalculate
implicit none
! ./coscalculate BiFeO3_15K.cfg direction.vec 
integer :: npoints,npoints_test,i 
character(len=50) :: filename_cfg,filename_vec
double precision,allocatable:: spinx(:),spiny(:),spinz(:),vecx(:),vecy(:),vecz(:)
double precision,allocatable :: cosresult(:)

call get_command_argument(1,filename_cfg)
call get_command_argument(2,filename_vec)

if(command_argument_count()==0) then
  write(6,*) 'You forgot to give the input file.'
  stop
end if

open(1,file=trim(adjustl(filename_cfg)),form='formatted',status='unknown')
open(2,file=trim(adjustl(filename_vec)),form='formatted',status='unknown')

! read in cfg file
read(1,*)
read(1,*) npoints
allocate(spinx(npoints))
allocate(spiny(npoints))
allocate(spinz(npoints))
allocate(vecx(npoints))
allocate(vecy(npoints))
allocate(vecz(npoints))
allocate(cosresult(npoints))

do i = 1,npoints
  read(1,*) spinx(i),spiny(i),spinz(i)
end do

! read in vec file
read(2,*) npoints_test
if(npoints_test /= npoints) then
  write(6,*) 'Number of points in the input files are not equal.'
  stop
end if
do i = 1,npoints
  read(2,*) vecx(i),vecy(i),vecz(i)
end do

! calculation
do i = 1,npoints
  cosresult(i) = spinx(i)*vecx(i) + spiny(i)*vecy(i) + spinz(i)*vecz(i)
end do

! write out result
open(3,file='cosresult.csv',form='formatted',status='unknown')
write(3,'(i8)') npoints
do i = 1,npoints
  write(3,'(i8,a,F12.6)') i,',',cosresult(i)
end do

end program coscalculate
