!***** Code: number of backgrounds in which mutation is beneficial
!***** Input: Number of mutations N, in main program
!***** Output: Column 1: x, Column 2: n_b (total), Column 2+i: n_b for backgrounds with i mutations. 

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
  integer,intent(in)::N,str(1:N)
  integer::i,dec
!common N

  dec=0
  do i=1,n
     dec=dec+str(n-i+1)*(2**(i-1))
  enddo


end function dec

function neighbor(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::neighbor(1:N)

    neighbor=str
    neighbor(i)=abs(1-str(i))


  end function neighbor

function bitflip(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::bitflip(1:N)

    bitflip=str
    bitflip(i)=abs(1-str(i))

  end function bitflip


  ! function ifpeak(str,r,N)
  !   integer,intent(in)::N,str(1:N)
  !   double precision,intent(in)::r(0:2**N-1)
  !   integer::i,j,npk,ifpeak,ii,nstr(1:N)
    
  !  do ii=1,N
  !     nstr=neighbor(str,ii,N)
  !     i=dec(nstr,N)
  !     j=dec(str,n)
  !     if (r(j)>r(i)) then
  !        exit
  !     else
  !        npk=npk+1
  !     endif
  !  enddo
   
  !     if (npk==N) ifpeak=1

  !    end function ifpeak
  



     function sizestr(str,N)
    integer,intent(in)::N,str(1:N)
    integer::i,sizestr
 sizestr=0
    do i=1,N
       sizestr=sizestr+str(i)
    enddo

  end function sizestr

  
     
       
       
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


 double precision function choose(n,m)
  integer,intent(in)::m,n
  integer::i,j

  choose=1.d0*fac(n)/(1.d0*fac(m)*fac(n-m))

end function choose

double precision function fac(n)
  integer,intent(in)::n
fac=1.d0
  if (n>0) then

  do i=1,n
     fac=fac*1.d0*i
  enddo

  endif

end function fac
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

  integer,parameter::N=5,NN=40,kexp=1,npop=1000
  integer::i,j,ii,m,a,kk,jj,irea,rea=1600,reacount(1:NN,1:2**n)=0,fg=0,kp=100
integer::numpks,pklist(1:2**N)=0,szint=0,ij,prev,nxt,ik,ifpeak=0,bt=0
  double precision::setot=0.d0,benf=0.d0,benftot=0.d0,semax=0.d0,signep(1:N)=0.d0,exct=0.d0,benft=0.d0,kp2,reaf(1:NN)=0.d0,s1=0.d0
  integer::str(1:N)=0,seed=-18974946,npk=0,nstr(1:N)=0,pkindicator(0:2**n-1)=0,str1(1:N)=0,str2(1:N)=0
  double precision::rw=1.d0,rsingle(1:N)=0.d0,micw=1.d0,micsingle(1:N)=0.d0,r0(0:2**N-1)=1.d0,r(0:2**N-1)=1.d0
double precision::rch(1:NN,1:2**n)=0.d0,ptot(1:NN)=0.d0,rch2(1:NN,1:2**n)=0.d0,recptot(1:NN)=0.d0
  double precision::xu=0.d0,xl=0.d0,dx(1:NN)=.001d0,rnum2=0.d0,cs2=0.d0,kn(1:NN,0:N)=0.d0,logmu=0.d0,xx(1:NN)=0.d0,u1,u2,meanft
double precision::mic0(0:2**N-1)=1.d0,x=0.d0,rr(0:1)=0.d0,eps=.0000001d0,rnum,pk(1:NN)=0.d0,cdav=0.d0,cs(1:NN)=0.d0,pk2=0.d0
double precision::gam=0.d0,gam2=0.d0,gam1=0.d0,bt2=0.d0,acrea(1:NN)=0.d0
integer::k=11,kj
integer::pkx,peakrank(1:2**n),jk=0,pmax=0,grank(1:2**n)=0

double precision::nbtot(1:NN)=0.d0,nbn(1:NN,0:N-1)=0.d0,septot(1:NN)=0.d0,sepn(1:NN,0:N-1)=0.d0,nbn1(0:N-1)=0.d0,nbtot1=0.d0

double precision::a1=1.5d0,a2=.1d0,r0a=.5d0,b=.95d0,ae=.8d0,aa=1.d0,x0=1.d0,pi=acos(-1.d0),kin=0.d0,sizefittest(1:NN)=0.d0
double precision::alph=1.d0
integer::alpha=2,beta=2,nup=2

30 format(f48.9,3x,f16.7)
character(len=12)::flnm
xl=0.d0
xu=20.d0

!dx=(1.5d0*N)/(1.d0*NN)
do i=1,NN
dx(i)=N*.25445d0*sqrt(.1d0)
enddo
ptot=0.d0
rch=0.d0
recptot=0.d0
nbn=0.d0
!************************ loop over realizations
do irea=1,rea
!************ Assigning MICs and r's to mutations
  
  do i=1,N
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
   !  rsingle(i)=ran2(seed)
! assigning MIC
 gam=0.d0
     do ii=1,2
        gam=gam-alph*log(1.d0-ran2(seed))
     enddo

     micsingle(i)=1.d0/(sqrt(rsingle(i)))+gam
     
end do



!**************** Assigning fitness and mic to genotypes
r0=1.d0
mic0=1.d0

  do i=1,2**n-1
  
     str=bitstr(i,N)
     do m=1,n
       if (str(m)==1)  r0(i)=r0(i)*rsingle(m)
     enddo


     do m=1,n
         if (str(m)==1) mic0(i)=mic0(i)*micsingle(m)
     enddo


  enddo
!***************** end


!***************** LOOP OVER CONCENTRATIONS 
  do kj=1,NN
     !   print*," "
     x=2.d0**(kj*dx(kj))-1.d0
     xx(kj)=x

     !***** Assigning growth rates to genotypes at current x
     do i=0,2**n-1
        r(i)=r0(i)*1.d0/(1.d0+(xx(kj)/mic0(i))**2)
     end do
     !**** End    
     
nbtot1=0.d0
nbn1=0.d0
     do i=1,2**N-1,2
!print*,xx(kj),i,sizestr(bitstr(i-1,N),N)
        if (r(i)>r(i-1)) then
           nbtot(kj)=nbtot(kj)+1.d0
           nbn(kj,sizestr(bitstr(i-1,N),N))=nbn(kj,sizestr(bitstr(i-1,N),N))+1.d0
!print*,xx(kj),i,size(bitstr(i-1,N),N)
            nbtot1=nbtot1+1.d0
           nbn1(sizestr(bitstr(i-1,N),N))=nbn1(sizestr(bitstr(i-1,N),N))+1.d0
           endif
     enddo

     septot(kj)=septot(kj)+(1.d0*nbtot1)/(1.d0*2**(N-1))*(1.d0-nbtot1/(1.d0*2**(N-1)))
    

     do i=0,n-1
  sepn(kj,i)=sepn(kj,i)+nbn1(i)/choose(n-1,i)*(1.d0-nbn1(i)/choose(n-1,i))   
enddo


enddo
!*************** conc. loop ends
  
  
enddo
!************** realization loop ends
do kj=1,NN
     do i=0,n-1
  sepn(kj,i)=sepn(kj,i)
enddo
enddo
do jj=1,NN


print*,xx(jj), nbtot(jj)/(1.d0*rea),nbn(jj,:)/(1.d0*rea),
!do i=1,N
!print*,i,nbn(jj,i)
!enddo


enddo

end program main




  
  

  
