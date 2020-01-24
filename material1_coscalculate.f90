program coscalculate
implicit none

integer :: npoints,npoints_test,i,ctemp,n1,n2,n,ierror
character(len=50) :: filename_vecbond,filename_vec,buffer
character(len=10) :: cbond,cdirection,cdirectionx,cdirectiony,cdirectionz
real,allocatable:: vecbondx(:),vecbondy(:),vecbondz(:),vecx(:),vecy(:),vecz(:)
real,allocatable :: cosresult(:),lengthbond(:)
real :: directionx,directiony,directionz
real :: length

ierror = 0

call get_command_argument(1,filename_vecbond)
call get_command_argument(2,cdirection)

! e.g. ./coscalculate *.vec [0,0,1]
cdirection = adjustl(cdirection)
n1 = index(cdirection,'[') ; n2 = index(cdirection(n1+1:),']')
if ((n1>0).and.(n2>0)) then
  cdirection = trim(adjustl(cdirection(n1+1:n2)))
end if
n = index(cdirection,',')
cdirectionx = cdirection(1:n-1)
cdirection = adjustl(cdirection(n+1:))
n = index(cdirection,',')
cdirectiony = cdirection(1:n-1)
cdirectionz = adjustl(cdirection(n+1:))
read(cdirectionx,*,iostat=ierror) directionx
read(cdirectiony,*,iostat=ierror) directiony
read(cdirectionz,*,iostat=ierror) directionz

if(command_argument_count()==0) then
  write(6,*) 'You forgot to give the input file.'
  stop
end if

open(1,file=trim(adjustl(filename_vecbond)),form='formatted',status='unknown')

! read in vec bond file and write out direction.vec file
read(1,*) cbond,npoints
allocate(vecbondx(npoints))
allocate(vecbondy(npoints))
allocate(vecbondz(npoints))
allocate(vecx(npoints))
allocate(vecy(npoints))
allocate(vecz(npoints))
allocate(cosresult(npoints))
allocate(lengthbond(npoints))

write(6,'(i6,4x,a,f2.0,a,f2.0,a,f2.0,a)') npoints,'[',directionx,',',directiony,',',directionz,']'
length = sqrt(directionx*directionx + directiony*directiony + directionz*directionz)
directionx = directionx/length
directiony = directiony/length
directionz = directionz/length
write(6,'(a,3f2.0)')'Normalised direction: ',directionx,directiony,directionz

do i = 1,npoints
  read(1,*) ctemp,ctemp,vecbondx(i),vecbondy(i),vecbondz(i)
  lengthbond(i) = sqrt(vecbondx(i)*vecbondx(i)&
 +vecbondy(i)*vecbondy(i)+vecbondz(i)*vecbondz(i))
  vecbondx(i) = vecbondx(i)/lengthbond(i)
  vecbondy(i) = vecbondy(i)/lengthbond(i)
  vecbondz(i) = vecbondz(i)/lengthbond(i)
end do

cosresult = 0

! calculation
do i = 1,npoints
  cosresult(i) = vecbondx(i)*directionx&
+ vecbondy(i)*directiony + vecbondz(i)*directionz
end do

! write out result
open(3,file='cosresult.csv',form='formatted',status='new')
write(3,'(i8,8x,3f2.0)') npoints,directionx,directiony,directionz
do i = 1,npoints
  write(3,'(i8,a,F12.6)') i,',',cosresult(i)
end do

close(1)
close(2)
close(3)
end program coscalculate
