program dipole_O4Cu
implicit none
!input file *.mean file
integer :: i,j,iatom_34,iatom_56,count_34,count_56
integer,allocatable :: jatom(:)
character(len=50) :: filename
character(len=30) :: ctemp
real,allocatable :: dipolevec(:,:),mean_34(:),mean_56(:),sqmean_34(:),sqmean_56(:)
real,allocatable :: dev_34(:),dev_56(:)

call get_command_argument(1,filename)
if(command_argument_count()==0) then
  write(6,*) 'You forgot to give the filename'
  stop
end if

allocate(dipolevec(6600,3))
allocate(mean_34(3))
allocate(mean_56(3))
allocate(sqmean_34(3))
allocate(sqmean_56(3))
allocate(dev_34(3))
allocate(dev_56(3))
allocate(jatom(6600))

open(1,file=trim(adjustl(filename)),form='formatted',status='unknown')
open(2,file='O34.sort',form='formatted',status='unknown')
open(3,file='O56.sort',form='formatted',status='unknown')
open(4,file='dipole_O4Cu.meandev',form='formatted',status='unknown')
! 6600 O
do i = 1,6
read(1,*) ctemp
end do

mean_34 = 0
mean_56 = 0
sqmean_34 = 0
sqmean_56 = 0
count_34 = 0
count_56 = 0
jatom = 0
dipolevec = 0
dev_34 = 0
dev_56 = 0

do i = 1,6600
  read(1,*) jatom(i),ctemp,dipolevec(i,:)
end do
do i = 1,3300
  read(2,*) iatom_34
  read(3,*) iatom_56
  do j = 1,6600
    if(jatom(j) == iatom_34) then
      mean_34(1) = mean_34(1) + dipolevec(j,1)
      mean_34(2) = mean_34(2) + dipolevec(j,2)
      mean_34(3) = mean_34(3) + dipolevec(j,3)
      sqmean_34(1) = sqmean_34(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_34(2) = sqmean_34(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_34(3) = sqmean_34(3) + dipolevec(j,3)*dipolevec(j,3)
      count_34 = count_34 + 1
    else if(jatom(j) == iatom_56) then
      mean_56(1) = mean_56(1) + dipolevec(j,1)
      mean_56(2) = mean_56(2) + dipolevec(j,2)
      mean_56(3) = mean_56(3) + dipolevec(j,3)
      sqmean_56(1) = sqmean_56(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_56(2) = sqmean_56(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_56(3) = sqmean_56(3) + dipolevec(j,3)*dipolevec(j,3)
      count_56 = count_56 + 1
    end if
  end do
end do
mean_34(1) = mean_34(1)/count_34
mean_34(2) = mean_34(2)/count_34
mean_34(3) = mean_34(3)/count_34

mean_56(1) = mean_56(1)/count_56
mean_56(2) = mean_56(2)/count_56
mean_56(3) = mean_56(3)/count_56

sqmean_34(1) = sqmean_34(1)/count_34
sqmean_34(2) = sqmean_34(2)/count_34
sqmean_34(3) = sqmean_34(3)/count_34

sqmean_56(1) = sqmean_56(1)/count_56
sqmean_56(2) = sqmean_56(2)/count_56
sqmean_56(3) = sqmean_56(3)/count_56

dev_34(1) = sqmean_34(1) - mean_34(1)*mean_34(1)
dev_34(2) = sqmean_34(2) - mean_34(2)*mean_34(2)
dev_34(3) = sqmean_34(3) - mean_34(3)*mean_34(3)
dev_56(1) = sqmean_56(1) - mean_56(1)*mean_56(1)
dev_56(2) = sqmean_56(2) - mean_56(2)*mean_56(2)
dev_56(3) = sqmean_56(3) - mean_56(3)*mean_56(3)


write(4,*) 'O3,4-Cu1,2,7,8'
write(4,'(a,3f12.6)') 'mean_34: ',mean_34
write(4,'(a,3f12.6)') 'dev_34 = sqmean-meansq: ',dev_34
write(4,*) 'O5,6-Cu1,2,7,8'
write(4,'(a,3f12.6)') 'mean_56: ',mean_56
write(4,'(a,3f12.6)') 'dev_56 = sqmean-meansq: ',dev_56

close(1)
close(2)
close(3)
close(4)
end program
