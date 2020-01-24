program bonds_CuO
implicit none
! ./bonds_CuO *.rmc6f
integer :: i,j,filename_len,file_err,natoms,n1,bond_count,k,m,n
integer,allocatable :: iatom(:),cell_nx(:),cell_ny(:),cell_nz(:),atom_number(:)
integer :: ncellx,ncelly,ncellz 
character(len=100) :: filename,filename_root,filename_ex,ctemp
character(len=2),allocatable :: catom(:)
real,allocatable :: xf(:),yf(:),zf(:)
real :: cell(3,3),dotvec,modvec1,modvec2,cosbond
real :: dis,dx,dy,dz,dxo,dyo,dzo
logical,allocatable :: lbond(:,:)
real,allocatable :: distance(:,:)
real,allocatable :: vec(:,:,:)

call get_command_argument(1,filename)
if(command_argument_count()==0) then
  write(6,*) 'You forgot to give the filename'
  stop
end if

filename = trim(adjustl(filename))
filename_len = len_trim(filename)

if (filename_len == 0) stop
i = index(filename,'.')

if(i /= 0) then
  filename_root = filename(1:i-1)
  filename_ex = filename(i+1:filename_len)
end if

open(1,file=trim(adjustl(filename)),form='formatted',status='unknown')

do i = 1,6 
  read(1,*) 
end do
read(1,'(a)') ctemp
n1 = index(ctemp,':')
ctemp = trim(adjustl(ctemp(n1+1:)))
read(ctemp,*) natoms
read(1,*)
read(1,*) ctemp,ctemp,ncellx,ncelly,ncellz 
read(1,*)
read(1,*)
do i = 1,3
  read(1,*) cell(i,:)
end do
read(1,*)

allocate(catom(natoms))
allocate(xf(natoms))
allocate(yf(natoms))
allocate(zf(natoms))
allocate(cell_nx(natoms))
allocate(cell_ny(natoms))
allocate(cell_nz(natoms))
allocate(iatom(natoms))
allocate(atom_number(natoms))
allocate(lbond(natoms,natoms))
allocate(distance(natoms,natoms))
allocate(vec(natoms,natoms,3))

do i = 1,natoms
  read(1,*) iatom(i),catom(i),xf(i),yf(i),zf(i),atom_number(i),cell_nx(i),cell_ny(i),cell_nz(i)
end do
! '[101]'
! 'Cu-O,1-4 1-5,2-3 2-6,7-4 7-5,8-3 8-6'
! '[10-1]'
! 'Cu-O,1-3 1-6,2-4 2-5,7-3 7-6,8-4 8-5' 

lbond = .false.
! Calculations:
do i = 1,6600
  do j = 6601,13200
      dx = xf(j) - xf(i) + 1.5d0
      dx = dx - aint(dx) - 0.5d0
      dy = yf(j) - yf(i) + 1.5d0
      dy = dy - aint(dy) - 0.5d0
      dz = zf(j) - zf(i) + 1.5d0
      dz = dz - aint(dz) - 0.5d0
      dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
      dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
      dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
      dis = sqrt(dxo**2 + dyo**2 + dzo**2)
      if (dis >2.16) cycle
      if(dis>0) then
        lbond(i,j) = .true.
        distance(i,j) = dis
!vec j - i
      vec(i,j,1) = dxo
      vec(i,j,2) = dyo
      vec(i,j,3) = dzo
      end if
  end do
end do
! write out
open(2,file='CuO145.dis',form='formatted',status='unknown')
open(3,file='CuO136.dis',form='formatted',status='unknown')
open(4,file='CuO236.dis',form='formatted',status='unknown')
open(18,file='CuO245.dis',form='formatted',status='unknown')
open(19,file='CuO745.dis',form='formatted',status='unknown')
open(20,file='CuO736.dis',form='formatted',status='unknown')
open(8,file='CuO836.dis',form='formatted',status='unknown')
open(9,file='CuO845.dis',form='formatted',status='unknown')
open(10,file='CuOCu147.ang',form='formatted',status='unknown')
open(11,file='CuOCu157.ang',form='formatted',status='unknown')
open(12,file='CuOCu137.ang',form='formatted',status='unknown')
open(13,file='CuOCu167.ang',form='formatted',status='unknown')
open(14,file='CuOCu238.ang',form='formatted',status='unknown')
open(15,file='CuOCu268.ang',form='formatted',status='unknown')
open(16,file='CuOCu248.ang',form='formatted',status='unknown')
open(17,file='CuOCu258.ang',form='formatted',status='unknown')

do i = 1,6600
  do j = 6601,13200
      if(lbond(i,j) .eqv. .false.) cycle
      if(atom_number(i) == 1) then
        if((atom_number(j)==4).or.(atom_number(j)==5)) then
          write(2,'(f12.6)') distance(i,j)
        end if
        if((atom_number(j)==3).or.(atom_number(j)==6)) then
          write(3,'(f12.6)') distance(i,j)
        end if
      end if

      if(atom_number(i) == 2) then
        if((atom_number(j)==3).or.(atom_number(j)==6)) then
          write(4,'(f12.6)') distance(i,j)
        end if
        if((atom_number(j)==4).or.(atom_number(j)==5)) then
          write(18,'(f12.6)') distance(i,j)
        end if
      end if
     
      if(atom_number(i) == 7) then
        if((atom_number(j)==4).or.(atom_number(j)==5)) then
          write(19,'(f12.6)') distance(i,j)
        end if
        if((atom_number(j)==3).or.(atom_number(j)==6)) then
          write(20,'(f12.6)') distance(i,j)
        end if
      end if

      if(atom_number(i) == 8) then
        if((atom_number(j)==3).or.(atom_number(j)==6)) then
          write(8,'(f12.6)') distance(i,j)
        end if
        if((atom_number(j)==4).or.(atom_number(j)==5)) then
          write(9,'(f12.6)') distance(i,j)
        end if
      end if
  end do
end do

do i = 1,6600
  do j = 6601,13200
    if(lbond(i,j) .eqv. .false.) cycle
    do k = 1,6600
      if((lbond(k,j) .eqv. .false.) .or.(lbond(i,j) .eqv. .false.)) cycle
      dotvec = vec(i,j,1)*vec(k,j,1)+vec(i,j,2)*vec(k,j,2)+vec(i,j,3)*vec(k,j,3)
      modvec1 = sqrt(vec(i,j,1)**2+vec(i,j,2)**2+vec(i,j,3)**2)
      modvec2 = sqrt(vec(k,j,1)**2+vec(k,j,2)**2+vec(k,j,3)**2)
      cosbond = dotvec/(modvec1*modvec2)
      if((atom_number(i)==1).and.(atom_number(k)==7)) then
        if(atom_number(j)==4) write(10,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==5) write(11,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==3) write(12,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==6) write(13,'(f12.6)') acos(cosbond)*180/3.1416
      end if  
      if((atom_number(i)==2).and.(atom_number(k)==8)) then
        if(atom_number(j)==3) write(14,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==6) write(15,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==4) write(16,'(f12.6)') acos(cosbond)*180/3.1416
        if(atom_number(j)==5) write(17,'(f12.6)') acos(cosbond)*180/3.1416
      end if
    end do
  end do
!if(mod(i,100)==0) write(6,'(a)',advance='no')'.'
end do

close(1)
close(2)
close(3)
close(4)
close(18)
close(19)
close(20)
close(8)
close(9)
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
end program
