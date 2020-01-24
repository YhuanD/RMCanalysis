! This program calculate the distribution of displacement of F from its
! neighbouring Sc-Sc bond
! e.g. ./anharmonicity Sc_F*.vec Sc_Sc*.vec
program anharmonicity
implicit none

character(len=40) :: ScFvecfile,ScScvecfile
character(len=20) :: ctemp
integer :: nbonds1,nbonds2,i,j,k
integer :: icount,ntemp
integer,allocatable :: nF(:),nSc1(:),nSc2(:),ScSc1(:),ScSc2(:),Sc(:),f(:)
double precision,allocatable :: ScFvec(:,:),ScScvec(:,:)
double precision :: distemp,dotpro,frac,shadow,modScSc,fracvecx,fracvecy,fracvecz
double precision :: fracvecx_new,fracvecy_new,fracvecz_new,ScFvec_newx &
,ScFvec_newy,ScFvec_newz,modtemp
double precision :: m11,m12,m13,m21,m22,m23,m31,m32,m33,u,v,w,sintheta,costheta
double precision,allocatable :: disx(:),disy(:),disz(:),dismod(:)
logical :: ltemp

call get_command_argument(1,ScFvecfile)
call get_command_argument(2,ScScvecfile)

if(command_argument_count()==0)then
  write(6,*) 'You forgot to give the filename'
  stop
end if

open(1,file=trim(adjustl(ScFvecfile)),form='formatted',status='unknown')
open(2,file=trim(adjustl(ScScvecfile)),form='formatted',status='unknown')
open(3,file='displacement.csv',form='formatted',status='unknown')

! read in ScFvecfile
read(1,*) ctemp,nbonds1
read(2,*) ctemp,nbonds2
write(3,'(4a)') 'F,','Sc,','Sc,','disx,disy,disz,dismod'

allocate(nF(nbonds1))
allocate(nSc1(nbonds1))
allocate(nSc2(nbonds1))
allocate(disx(nbonds1))
allocate(disy(nbonds1))
allocate(disz(nbonds1))
allocate(dismod(nbonds1))
allocate(Sc(nbonds1))
allocate(f(nbonds1))
allocate(ScFvec(nbonds1,3))
allocate(ScSc1(nbonds2))
allocate(ScSc2(nbonds2))
allocate(ScScvec(nbonds2,3))

nF = 0
nSc1 = 0
nSc2 = 0
Sc = 0
f = 0
disx = 0
disy = 0
disz = 0
dismod = 0
icount = 0

do i = 1,nbonds1
  read(1,*) Sc(i),F(i),ScFvec(i,:)
end do

do i = 1,nbonds2
  read(2,*) ScSc1(i),ScSc2(i),ScScvec(i,:)
end do

do i = 1,nbonds1
  do j = i,nbonds1
    ltemp = (F(i) == F(j)).and.(Sc(i)/=Sc(j))
    if (ltemp) then
      do k = 1,nbonds2
        ltemp = (Sc(i) == ScSc1(k)).and.(Sc(j) == ScSc2(k))
        if(ltemp) then
          modScSc = sqrt(ScScvec(k,1)**2 + ScScvec(k,2)**2 + ScScvec(k,3)**2)
          dotpro = ScFvec(i,1)*ScScvec(k,1)&
                   +ScFvec(i,2)*ScScvec(k,2)&
                   +ScFvec(i,3)*ScScvec(k,3)
          shadow = dotpro/modScSc
          fracvecx = shadow*ScScvec(k,1)/modScSc
          fracvecy = shadow*ScScvec(k,2)/modScSc
          fracvecz = shadow*ScScvec(k,3)/modScSc
!write(6,*) 'frac,shadow,fracScScx,y,z: ',shadow/modScSc,shadow,fracvecx,fracvecy,fracvecz
! change coordinates, Sc-Sc bond to be z axis
          costheta = ScScvec(k,3)/modScSc
          sintheta = sqrt(ScScvec(k,1)**2 + ScScvec(k,2)**2)/modScSc
!write(6,*) 'costheta,sintheta',costheta,sintheta
          u = ScScvec(k,2)
          v = -ScScvec(k,1)
          w = 0.0
          modtemp = sqrt(u**2+v**2)
          u = u/modtemp
          v = v/modtemp
!write(6,*) 'u,v,w: ',u,v,w
          m11 = costheta + (u * u) * (1 - costheta)
          m12 = u * v * (1 - costheta) - w * sintheta
          m13 = u * w * (1 - costheta) + v * sintheta
          m21 = u * v * (1 - costheta) + w * sintheta
          m22 = costheta + v * v * (1 - costheta)
          m23 = w * v * (1 - costheta) - u * sintheta
          m31 = u * w * (1 - costheta) - v * sintheta
          m32 = v * w * (1 - costheta) + u * sintheta
          m33 = costheta + w * w * (1 - costheta)
!write(6,*) 'm11,m12,m13; m21,m22,m23; m31,m32,m33:&
!',m11,m12,m13,';',m21,m22,m23,';',m31,m32,m33
          fracvecx_new = m11*fracvecx + m12*fracvecy + m13*fracvecz
          fracvecy_new = m21*fracvecx + m22*fracvecy + m23*fracvecz
          fracvecz_new = m31*fracvecx + m32*fracvecy + m33*fracvecz
!write(6,*) 'fracnewxyz: ',fracvecx_new,fracvecy_new,fracvecz_new
          ScFvec_newx = m11*ScFvec(i,1) + m12*ScFvec(i,2) + m13*ScFvec(i,3)
          ScFvec_newy = m21*ScFvec(i,1) + m22*ScFvec(i,2) + m23*ScFvec(i,3)
          ScFvec_newz = m31*ScFvec(i,1) + m32*ScFvec(i,2) + m33*ScFvec(i,3)
          disx(i) = ScFvec_newx - fracvecx_new
          disy(i) = ScFvec_newy - fracvecy_new
          disz(i) = ScFvec_newz - fracvecz_new 
          dismod(i) = sqrt(disx(i)**2+disy(i)**2+disz(i)**2)
!write(6,*) 'modtemp: ',modtemp,'dismod(i): ',dismod(i)
          icount = icount + 1
          nF(i) = F(i)
          nSc1(i) = Sc(i)
          nSc2(i) = Sc(j)
        end if
      end do
    end if
  end do
end do

write(3,'(i0)') icount

do i = 1,nbonds1
ltemp = (disx(i)/=0.0) .or. (disy(i)/=0.0)
  if(ltemp) then
    do j = i+1,nbonds1
      if(nF(i) > nF(j)) then
        ntemp = nF(i)
        nF(i) = nF(j)
        nF(j) = ntemp
        
        ntemp = nSc1(i)
        nSc1(i) = nSc1(j)
        nSc1(j) = ntemp
        
        ntemp = nSc2(i)
        nSc2(i) = nSc2(j)
        nSc2(j) = ntemp
        
        distemp = disx(i)
        disx(i) = disx(j)
        disx(j) = distemp
        distemp = disy(i)
        disy(i) = disy(j)
        disy(j) = distemp
        distemp = disz(i)
        disz(i) = disz(j)
        disz(j) = distemp 
        distemp = dismod(i)
        dismod(i) = dismod(j)
        dismod(j) = distemp
      end if
    end do
  end if
end do

do i = 1,nbonds1 
ltemp = (disx(i)/=0.0) .or. (disy(i)/=0.0)
  if (ltemp) then
    write(3,'(i0,a,i0,a,i0,a,f12.6,a,f12.6,a,f12.6,a,f12.6)') nF(i),',',nSc1(i),',',nSc2(i)&
          ,',',disx(i),',',disy(i),',',disz(i),',',dismod(i)
  end if
end do

close(1)
close(2)
close(3)
end program
