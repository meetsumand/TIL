!**** Code: to compute total number of peaks with n mutations as a function of x.
!**** Input: Total mutation number N (in main program)
!**** Output: To files with names of the form "pkflN20n12.d", where N= total n of mutations, n=number of of mutations in peak
!**** Output: Column 1: x, Column 2: mean number of peaks

module strings
contains
  
function bitstr(i,N)
  integer,intent(in)::N
  integer::i,j,bitstr(1:N)
  
  do j=1,N
     k=N-j+1
     bitstr(j)=(mod(i,2**(k))-mod(i,2**(k-1)))
     bitstr(j)=bitstr(j)/(2**(k-1))
  enddo


end function bitstr



function insertbit(str,i,ii,N)
  integer,intent(in)::N
  integer::str(1:N),i,ii,insertbit(1:N)

  insertbit=str
  
  do j=2,i
     k=str(j)
     insertbit(j-1)=k
  enddo

  insertbit(i)=ii

end function insertbit


function dec(str,N)
  integer,intent(in)::N
  integer::str(1:N),i,dec
!common N

  dec=0
  do i=1,n
     dec=dec+str(n-i+1)*(2**(i-1))
  enddo


end function dec

function neighbor(str,i,N)
    integer,intent(in)::N,i
    integer::str(1:N),neighbor(1:N)

    neighbor=str
    neighbor(i)=abs(1-str(i))


  end function neighbor

  function size(str,N)
    integer,intent(in)::N
    integer::str(1:N),i
size=0.d0
    do i=1,N
       size=size+1.d0*str(i)
    enddo

  end function size

  
    
       
FUNCTION ran2(idum)
!    gerador de números aleatórios baseado na combinacao
!    de dois geradores lineares congruenciais tendo um periodo
!    maior do que 2x10^18. A saida e' "baralhada" 
!     -> valor inicial de idum=-1  (Numerical recipes 2a. edicao)
IMPLICIT NONE
INTEGER, PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,&
     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
     ir2=3791,ntab=32,ndiv=1+imm1/ntab
DOUBLE PRECISION , PARAMETER ::   am=1.d0/im1,eps=1.d-14,rnmx=1.d0-eps
DOUBLE PRECISION :: ran2
INTEGER, DIMENSION(ntab) :: iv
INTEGER :: idum,idum2,j,k,iy

save iv,iy,idum2
data idum2/123456789/,iv /ntab*0/,iy /0/
      
if(idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=ntab+8,1,-1
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-ir1*k
    if(idum.lt.0) idum=idum+im1
    if(j.le.ntab) iv(j)=idum
   end do
   iy=iv(1)
endif

k=idum/iq1
idum=ia1*(idum-k*iq1)-ir1*k
if(idum.lt.0) idum=idum+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-ir2*k

if(idum2.lt.0) idum2=idum2+im2

j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1)iy=iy+imm1
ran2=min(am*iy,rnmx)

END FUNCTION ran2


end module










!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM











program main
  
  use strings
  
  implicit none

  integer,parameter::N=12,NN=10,kexp=1
  integer::i,j,ii,m,a,kk,jj,irea,rea=1000,fg,kp=100,numpks,pklist(1:2**N)=0,szint=0
  double precision::setot=0.d0,benf=0.d0,benftot=0.d0,semax=0.d0,signep(1:N)=0.d0,exct=0.d0,benft=0.d0,kp2,reaf(1:NN)=0.d0
  integer::str(1:N),seed=-78974946,npk=0,nstr(1:N)=0
  double precision::rw=1.d0,rsingle(1:N)=0.d0,micw=1.d0,micsingle(1:N)=0.d0,r0(0:2**N-1)=1.d0,r(0:2**N-1)=0.d0
  double precision::xu=0.d0,xl=0.d0,dx=.001d0,rnum2=0.d0,cs2=0.d0,kn(1:NN,0:N)=0.d0,logmu=0.d0,xx(1:NN)=0.d0
  double precision::gam1=0.d0,gam2=0.d0
double precision::mic0(0:2**N-1)=1.d0,x=0.d0,rr(0:1)=0.d0,eps=.0000001d0,rnum,pk(1:NN)=0.d0,cdav=0.d0,cs(1:NN)=0.d0,pk2=0.d0
double precision::gam=0.d0,alph=1.d0,pktl(1:NN)=0.d0
integer::alpha=2,beta=2,nup=2

integer::k=11
double precision::a1=1.5d0,a2=.1d0,r0a=.5d0,b=.95d0,ae=.8d0,aa=1.d0,x0=1.d0


30 format(f48.9,3x,f16.7)
character(len=12)::flnm
character(len=9)::flnm2
xl=0.d0
xu=20.d0

dx=(2.5d0*N)/(1.d0*NN)
!kp2=((xu-xl)/dx)!**(1.d0/3.d0)
!kp=kp2
pk=0.d0
cs=0.d0
reaf=0.d0

do irea=1,rea
   numpks=1
   pklist(1)=0
  !************* Assigning growth rate to mutations



!************ Assigning MICs and r's to mutations

   do i=1,N
gam=0.d0

     ! assigning r
     gam1=0.d0
     do ii=1,alpha
        gam1=gam1-log(1.d0-ran2(seed))
     enddo

     gam2=0.d0
     do ii=1,beta
        gam2=gam2-log(1.d0-ran2(seed))
     enddo
rsingle(i)= gam1/(gam1+gam2)


     ! assigning MIC
     do ii=1,nup
        gam=gam-alph*log(1.d0-ran2(seed))
     enddo

     micsingle(i)=1.d0/sqrt(rsingle(i))+gam


     
end do

  !**************** Assigning fitness and mic to genotypes
  mic0=1.d0
  r0=1.d0
 do i=0,2**n-1
  
     str=bitstr(i,N)
     do m=1,n
        if(str(m)==1) r0(i)=r0(i)*rsingle(m)
     enddo


do m=1,n
   if (str(m)==1) mic0(i)=mic0(i)*micsingle(m)
enddo


enddo
!***************** end

do kk=1,NN
   x=2.d0**(xl+kk*dx)-1.d0
   xx(kk)=x
   npk=0
   !***** Assigning growth rates to genotypes at current x
do i=0,2**n-1
     r(i)=r0(i)*1.d0/(1.d0+(x/mic0(i))**2)
!r(i)=r0(i)*exp(-x/mic0(i))
end do
!**** End

  


  !*************** finding peaks

fg=0
pk2=0.d0
cs2=0.d0

do j=0,2**(n)-1
npk=0
str=bitstr(j,N)
   do ii=1,N
      nstr=neighbor(str,ii,N)
      a=dec(nstr,N)
      if (r(j).le.r(a)) then
         exit
      else
         npk=npk+1
      endif

   enddo
        
     if (npk==N) then
        pk(kk)=pk(kk)+1.d0
        pk2=pk2+1.d0
        cs2=cs2+size(str,N)
        szint=size(str,N)
        kn(kk,szint)=kn(kk,szint)+1.d0
        fg=0
        if (irea==rea) then
        do ii=1,numpks
           if (pklist(ii)==j) then
              fg=1
              exit
           endif

        enddo

        

        if (fg==0) then
           numpks=numpks+1
           pklist(numpks)=j
        endif

     endif

  endif

  

  enddo
  if (pk2>0.d0) then
     cs(kk)=cs(kk)+cs2/pk2
     reaf(kk)=reaf(kk)+1.d0
  endif

  
  enddo

!************* end
  enddo



!do kk=1,NN
!   print*,2.d0**(xl+kk*dx)-1.d0,pk(kk)/(1.d0*rea),cs(kk)/(1.d0*reaf(kk))
!enddo

!print*,numpks

  pktl=0.d0
do kk=1,NN

   do j=0,N
      pktl(kk)=pktl(kk)+kn(kk,j)/(1.d0*rea)
   enddo

   
enddo

WRITE(flnm2,fmt='(A5,I2.2,A2)') 'pkfTN',N,'.d'
OPEN(UNIT=50,file=flnm2,status='unknown')

do kk=1,NN
 write(unit=50,fmt=30) xx(kk),pktl(kk)
enddo




do j=0,N
WRITE(flnm,fmt='(A5,I2.2,A1,I2.2,A2)') 'pkflN',N,'n',j,'.d'
OPEN(UNIT=10,file=flnm,status='unknown')  

do kk=1,NN
      write(unit=10,fmt=30) xx(kk),kn(kk,j)/(1.d0*rea)
   enddo

   close (unit=10)
   close (unit=50)

   enddo
end program main




  
  

  
