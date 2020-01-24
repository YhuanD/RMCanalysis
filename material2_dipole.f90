program dipole_O4Cu
implicit none
!input file *.mean file
integer :: i,j,iatom_3,iatom_4,iatom_5,iatom_6,count_3,count_4,count_5,count_6
integer,allocatable :: jatom(:)
character(len=50) :: filename
character(len=30) :: ctemp
real,allocatable :: dipolevec(:,:),mean_3(:),mean_4(:),mean_5(:),mean_6(:),sqmean_3(:),sqmean_4(:),sqmean_5(:),sqmean_6(:)
real,allocatable :: dev_3(:),dev_4(:),dev_5(:),dev_6(:)

call get_command_argument(1,filename)
if(command_argument_count()==0) then
  write(6,*) 'You forgot to give the filename'
  stop
end if

allocate(dipolevec(6600,3))
allocate(mean_3(3))
allocate(mean_4(3))
allocate(mean_5(3))
allocate(mean_6(3))
allocate(sqmean_3(3))
allocate(sqmean_4(3))
allocate(sqmean_5(3))
allocate(sqmean_6(3))
allocate(dev_3(3))
allocate(dev_4(3))
allocate(dev_5(3))
allocate(dev_6(3))
allocate(jatom(6600))

open(1,file=trim(adjustl(filename)),form='formatted',status='unknown')
open(2,file='O3.sort',form='formatted',status='unknown')
open(10,file='O4.sort',form='formatted',status='unknown')
open(3,file='O5.sort',form='formatted',status='unknown')
open(11,file='O6.sort',form='formatted',status='unknown')
open(4,file='dipole_O4Cu.meandev',form='formatted',status='unknown')
! 6600 O
do i = 1,6
read(1,*) ctemp
end do

mean_3 = 0
mean_4 = 0
mean_5 = 0
mean_6 = 0
sqmean_3 = 0
sqmean_4 = 0
sqmean_5 = 0
sqmean_6 = 0
count_3 = 0
count_4 = 0
count_5 = 0
count_6 = 0
jatom = 0
dipolevec = 0
dev_3 = 0
dev_4 = 0 
dev_5 = 0
dev_6 = 0

do i = 1,6600
  read(1,*) jatom(i),ctemp,dipolevec(i,:)
end do
do i = 1,1650
  read(2,*) iatom_3
  read(10,*) iatom_4
  read(3,*) iatom_5
  read(11,*) iatom_6
  do j = 1,6600
    if(jatom(j) == iatom_3) then
      mean_3(1) = mean_3(1) + dipolevec(j,1)
      mean_3(2) = mean_3(2) + dipolevec(j,2)
      mean_3(3) = mean_3(3) + dipolevec(j,3)
      sqmean_3(1) = sqmean_3(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_3(2) = sqmean_3(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_3(3) = sqmean_3(3) + dipolevec(j,3)*dipolevec(j,3)
      count_3 = count_3 + 1
    else if(jatom(j) == iatom_4) then
      mean_4(1) = mean_4(1) + dipolevec(j,1)
      mean_4(2) = mean_4(2) + dipolevec(j,2)
      mean_4(3) = mean_4(3) + dipolevec(j,3)
      sqmean_4(1) = sqmean_4(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_4(2) = sqmean_4(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_4(3) = sqmean_4(3) + dipolevec(j,3)*dipolevec(j,3)
      count_4 = count_4 + 1
    else if(jatom(j) == iatom_5) then
      mean_5(1) = mean_5(1) + dipolevec(j,1)
      mean_5(2) = mean_5(2) + dipolevec(j,2)
      mean_5(3) = mean_5(3) + dipolevec(j,3)
      sqmean_5(1) = sqmean_5(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_5(2) = sqmean_5(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_5(3) = sqmean_5(3) + dipolevec(j,3)*dipolevec(j,3)
      count_5 = count_5 + 1
    else if(jatom(j) == iatom_6) then
      mean_6(1) = mean_6(1) + dipolevec(j,1)
      mean_6(2) = mean_6(2) + dipolevec(j,2)
      mean_6(3) = mean_6(3) + dipolevec(j,3)
      sqmean_6(1) = sqmean_6(1) + dipolevec(j,1)*dipolevec(j,1)
      sqmean_6(2) = sqmean_6(2) + dipolevec(j,2)*dipolevec(j,2)
      sqmean_6(3) = sqmean_6(3) + dipolevec(j,3)*dipolevec(j,3)
      count_6 = count_6 + 1
    end if
  end do
end do
mean_3(1) = mean_3(1)/count_3
mean_3(2) = mean_3(2)/count_3
mean_3(3) = mean_3(3)/count_3

mean_4(1) = mean_4(1)/count_4
mean_4(2) = mean_4(2)/count_4
mean_4(3) = mean_4(3)/count_4

mean_5(1) = mean_5(1)/count_5
mean_5(2) = mean_5(2)/count_5
mean_5(3) = mean_5(3)/count_5

mean_6(1) = mean_6(1)/count_6
mean_6(2) = mean_6(2)/count_6
mean_6(3) = mean_6(3)/count_6

sqmean_3(1) = sqmean_3(1)/count_3
sqmean_3(2) = sqmean_3(2)/count_3
sqmean_3(3) = sqmean_3(3)/count_3

sqmean_4(1) = sqmean_4(1)/count_4
sqmean_4(2) = sqmean_4(2)/count_4
sqmean_4(3) = sqmean_4(3)/count_4

sqmean_5(1) = sqmean_5(1)/count_5
sqmean_5(2) = sqmean_5(2)/count_5
sqmean_5(3) = sqmean_5(3)/count_5

sqmean_6(1) = sqmean_6(1)/count_6
sqmean_6(2) = sqmean_6(2)/count_6
sqmean_6(3) = sqmean_6(3)/count_6

dev_3(1) = sqmean_3(1) - mean_3(1)*mean_3(1)
dev_3(2) = sqmean_3(2) - mean_3(2)*mean_3(2)
dev_3(3) = sqmean_3(3) - mean_3(3)*mean_3(3)

dev_4(1) = sqmean_4(1) - mean_4(1)*mean_4(1)
dev_4(2) = sqmean_4(2) - mean_4(2)*mean_4(2)
dev_4(3) = sqmean_4(3) - mean_4(3)*mean_4(3)

dev_5(1) = sqmean_5(1) - mean_5(1)*mean_5(1)
dev_5(2) = sqmean_5(2) - mean_5(2)*mean_5(2)
dev_5(3) = sqmean_5(3) - mean_5(3)*mean_5(3)

dev_6(1) = sqmean_6(1) - mean_6(1)*mean_6(1)
dev_6(2) = sqmean_6(2) - mean_6(2)*mean_6(2)
dev_6(3) = sqmean_6(3) - mean_6(3)*mean_6(3)

write(4,*) 'O3,4-Cu1,2,7,8'
write(4,'(a,3f12.6)') 'mean_3: ',mean_3
write(4,'(a,3f12.6)') 'dev_3 = sqmean-meansq: ',dev_3
write(4,'(a,3f12.6)') 'mean_4: ',mean_4
write(4,'(a,3f12.6)') 'dev_4 = sqmean-meansq: ',dev_4
write(4,*) 'O5,6-Cu1,2,7,8'
write(4,'(a,3f12.6)') 'mean_5: ',mean_5
write(4,'(a,3f12.6)') 'dev_5 = sqmean-meansq: ',dev_5
write(4,'(a,3f12.6)') 'mean_6: ',mean_6
write(4,'(a,3f12.6)') 'dev_6 = sqmean-meansq: ',dev_6
close(1)
close(2)
close(3)
close(4)
close(11)
close(10)
end program
