program perpendicular
! The calculation doesnt take account of atoms at the edges so, for ScF3, F are
! shared by two bonds, vertical value for two bonds are both considered and
! accumulated.

implicit none

character(len=20) :: ctemp
character(len=40) :: bondvecfile,flucvecfile
integer :: nbonds,i,j,nfluc
integer,allocatable :: ibond1(:),ibond2(:),ifluc(:)
real,allocatable :: bondx(:),bondy(:),bondz(:),flucx(:),flucy(:),flucz(:)&
,crossx(:),crossy(:),crossz(:),modcross(:),modbond(:),vertical(:)

call get_command_argument(1,bondvecfile)
call get_command_argument(2,flucvecfile)

if(command_argument_count()==0)then
  write(6,*) 'You forgot to give the filename'
  stop
end if

open(1,file=trim(adjustl(bondvecfile)),form='formatted',status='unknown')
open(2,file=trim(adjustl(flucvecfile)),form='formatted',status='unknown')
open(3,file='vertifluc.csv',form='formatted',status='unknown')

! read in bond file
read(1,*) ctemp,nbonds

allocate(ibond1(nbonds))
allocate(ibond2(nbonds))
allocate(crossx(nbonds))
allocate(crossy(nbonds))
allocate(crossz(nbonds))
allocate(modcross(nbonds))
allocate(modbond(nbonds))
allocate(vertical(nbonds))
allocate(bondx(nbonds))
allocate(bondy(nbonds))
allocate(bondz(nbonds))

do i = 1,nbonds
  read(1,*) ibond1(i),ibond2(i),bondx(i),bondy(i),bondz(i)
end do 

! read in fluc vec file
read(2,*) ctemp
read(2,*) nfluc

allocate(ifluc(nfluc))
allocate(flucx(nfluc))
allocate(flucy(nfluc))
allocate(flucz(nfluc))

do i =1,nfluc
  read(2,*) ctemp,ifluc(i),flucx(i),flucy(i),flucz(i) 
end do

! cross product to get the result
do i = 1,nbonds
  do j = 1,nfluc
    if(ibond2(i) == ifluc(j)) then
      crossx(i) = bondy(i)*flucz(j)-bondz(i)*flucy(j)
      crossy(i) = bondz(i)*flucx(j)-bondx(i)*flucz(j)
      crossz(i) = bondx(i)*flucy(j)-bondy(i)*flucx(j)
      modcross(i) = sqrt(crossx(i)*crossx(i)+crossy(i)*crossy(i)&
+crossz(i)*crossz(i))
      modbond(i) = sqrt(bondx(i)*bondx(i)+bondy(i)*bondy(i)&
+bondz(i)*bondz(i))
      vertical(i) = modcross(i)/modbond(i)
    end if
  end do
end do

! write out result
write(3,'(2a)') 'atom number in bondfile, ','vertical fluctuations'
write(3,'(a,i0)') 'Number of atoms calculated: ',nbonds
do i = 1,nbonds
  write(3,'(i0,a,f12.6)') i,',',vertical(i)
end do

close(1)
close(2)
close(3)
end program perpendicular
