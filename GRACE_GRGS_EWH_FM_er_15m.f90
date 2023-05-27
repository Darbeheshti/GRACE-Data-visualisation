 !updated on 13 aug 2012

!change in surface mass density from changes in the spherical harmonic geoid coefficients [Wahr et al., 1998] eq. 30
! we use grgs data, noneed for smoothin, so we do not have 2*pi and Wl in eq.30

! EWT versus Forwadrd model 

!a time avarage of 293 solutions is removed to each 10-days solution

!EWT in cm

PROGRAM GRACE_GRGS_EWH_FM_er_15m

use Somearrays

implicit none


!declare variables

real, parameter ::pi = 3.141592653589793,ro_ave=5517,ro_w=1000,R_smoothing=500000,deg2rad=pi/180,jazr=sqrt(2.0)

!         AE                  1/F                 GM                 OMEGA
!0.63781364600000E+070.29825765000000E+030.39860044150000E+150.72921150000000E-04
!EARTH 0.3986004415E+15 0.6378136460E+07

!gives EWT in cm

double precision, parameter ::a_earth=0.6378136460E+07*100

integer,parameter ::lmax=50

double precision :: sum_s,sum_e,b_wahr,cons_wahr,vcos
real*8 ::	         sum_i,minu

INTEGER :: i,j,l_s,m_s,il,k,nj

integer, dimension(:,:), allocatable :: degree,order

double precision, dimension(:,:), allocatable :: cnm,snm,ecnm,esnm,delta_sigma,pleg,uncer

double precision, dimension(:), allocatable :: landa,phi

integer, dimension(:), allocatable :: n_lov,m,n,ma,na

double precision, dimension(:), allocatable :: wahr,c_lov,s_lov,k_lov,c,s,ca,sa,ec,es


 CHARACTER(len=80) :: dummy
 character(len=6), dimension(:), allocatable :: chara
 character(len=45), dimension(:), allocatable :: names,namesw
 
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables for subroutine date_and_time

      integer time_array_0(8), time_array_1(8)
      real start_time, finish_time

! Mark the beginning of the program

      call date_and_time(values=time_array_0)
      start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate (names(293))
allocate (namesw(293))

OPEN(9,file='GRACE_names.DAT',form='formatted')

OPEN(10,file='GRACE_EWT_names_15m.DAT',form='formatted')

DO j=1,293

read (9, 109) names(j)
read (10, 109) namesw(j)

109 format(a45)


END DO

 close(9)

 close(10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!we do not need to use b_wahr, because we are using GRGS data
b_wahr=log(2.0)/(1-cos(R_smoothing/a_earth))

cons_wahr=a_earth*ro_ave/3
cons_wahr=cons_wahr/ro_w
   
allocate (c_lov(1:1025),s_lov(1:1025),n_lov(1:1025),k_lov(1:1025))
 
 k_lov=0
 n_lov=0
 c_lov=0
 s_lov=0



OPEN(19,file='Load_Love2_CM.dat',form='formatted')

DO i=1,1025

read (19, 196) n_lov(i),c_lov(i),s_lov(i),k_lov(i)

196 format(i6,1x,d17.12,1x,d17.12,1x,d17.12)


END DO

 close(19)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!number of lines of 10 days GRGS Grace coefficients file 1326

!   allocate arrays
    allocate (cnm(51,51),snm(51,51),ecnm(51,51),esnm(51,51),degree(51,51),order(51,51))
    allocate (c(1326),s(1326),n(1326),m(1326),ca(1326),sa(1326),na(1326),ma(1326),ec(1326),es(1326))
    allocate (chara(1326))


 cnm=0
 snm=0
 degree=0
 order=0

 n=0
 m=0
 c=0
 s=0

 na=0
 ma=0
 ca=0
 sa=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!making a grid of 15m by 15m for whole world

nj=718
allocate (phi(nj),landa(nj*2))
minu=.25

sum_i = 1

DO i=1,nj*2
landa(i)=sum_i
sum_i=landa(i)+minu
END DO

sum_i = -90

DO i=1,nj
phi(nj+1-i)=sum_i
sum_i=phi(nj+1-i)+minu
end DO


landa=deg2rad*landa

phi=deg2rad*phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!opening and reading avarage of 293 files

OPEN(111,file='avarage.DAT',form='formatted')

DO i=1,1326

read (111, 1) na(i),ma(i),ca(i),sa(i)

1 format(i5,i5,d19.12,d19.12)

 END DO

 close(111)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (delta_sigma(2*nj,nj))
allocate (uncer(2*nj,nj))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!opening and reading 293 files

DO k=1,293

print*,'file no',k
OPEN(k+300,file=names(k),form='formatted',status='old')

do i=1,5

read(k+300,*) dummy
    
END DO



DO i=1,1326

read (k+300, 101) chara(i),n(i),m(i),c(i),s(i),ec(i),es(i)
!to test reading
!if (k==1) print*, chara(i),n(i),m(i),c(i),s(i),ec(i),es(i)

101 format(a6,i5,i5,d19.12,d19.12,d11.4,d11.4)

 c(i)=c(i)-ca(i)
 s(i)=s(i)-sa(i)

        cnm(n(i)+1,m(i)+1) = c(i)
        snm(n(i)+1,m(i)+1) = s(i)
        ecnm(n(i)+1,m(i)+1) = ec(i)
        esnm(n(i)+1,m(i)+1) = es(i)

!it seems i never used degree and order in this program
        degree(n(i)+1,m(i)+1)=n(i)
        order(n(i)+1,m(i)+1)=m(i)

END DO

 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



OPEN(k+600,file=namesw(k),form='formatted',status='unknown')

!print*,names(k)

DO i=1,2*nj

DO j=1,nj

vcos=sin(phi(j))


delta_sigma(i,j)=0
uncer(i,j)=0

DO l_s=0,50


Plm=0

sum_s=0
sum_e=0

 call  legendre(vcos,l_s,1501)

il=l_s*(l_s+1)/2+1

!!!!!!!!!!!!!!!!!

DO m_s=0,l_s

if (m_s /= 0) Plm(il)=Plm(il)*jazr

sum_s=sum_s+Plm(il)*(cnm(l_s+1,m_s+1)*cos(m_s*landa(i))+snm(l_s+1,m_s+1)*sin(m_s*landa(i)))
sum_e=sum_e+(Plm(il)**2)*((ecnm(l_s+1,m_s+1)**2)*(cos(m_s*landa(i))**2)+(esnm(l_s+1,m_s+1)**2)*(sin(m_s*landa(i))**2))
il=il+1

end DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (l_s == 0) k_lov(l_s+1)=0.000
if (l_s == 1) k_lov(l_s+1)=0.027

delta_sigma(i,j)=delta_sigma(i,j)+(((2*l_s)+1)/(1+k_lov(l_s+1)))*sum_s
uncer(i,j)=uncer(i,j)+((((2*l_s)+1)/(1+k_lov(l_s+1))))**2*sum_e

end DO 


delta_sigma(i,j)=delta_sigma(i,j)*cons_wahr
uncer(i,j)=uncer(i,j)*(cons_wahr**2)

write (k+600, 119) landa(i)/deg2rad,phi(j)/deg2rad,delta_sigma(i,j),sqrt(uncer(i,j))
119 format(f5.1,1x,f5.1,1x,d19.12,1x,d19.12)

end DO

end DO

 close(k+300)
 close(k+600)
 
 
end DO

 

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 call date_and_time(values=time_array_1)
      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)







write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)/60 



end Program GRACE_GRGS_EWH_FM_er_15m







!subroutine legendre(X,ndeg,maxnm,Plm)
subroutine legendre(X,ndeg,maxnm)
!
!     calculate the normalised associated Legendre function values at x up to degree ndeg
!
!     The normalisation is such that the integral of Pnm^2 from -1 to 1 is 2.
!     This is slightly different from the normalisation used in the older
!     sea-level programs.  In those programs, if m is not equal to 0, the
!     normalisation had an extra factor sqrt(2), so that integral over the
!     sphere of Pnm cos or sin m phi is 4 pi.  With this normalisation
!     the integral over the sphere of Pnm exp(i m phi) is equal to 4 pi.
!
! PT070326: compute only to the maximum degree of the input file, ndeg, rather
!           than to the maximum dimensioned degree.
!

use Somearrays
     
      implicit none


! input cos(co-latitude) in radians
 ! real*8 X
double precision X

! local variables
  !real*8 aa,as,s,sum,c,d
double precision aa,as,s,sum,c,d

      integer i,j,l,m,l0,i0,i1,i2,ndeg,maxnm
 !     real*8 Plm(maxnm)

      Plm(1)=1
      Plm(2)=X
      i1=1
      i2=2
      do i=2,ndeg
        i0=i1
        i1=i2
        i2=i*(i+1)/2+1
        Plm(i2)=((2*i-1)*X*Plm(i1)-(i-1)*Plm(i0))/i 
         if(i2.eq.2)print*,'Plm(2)',Plm(2)
!        print*,'legendre: ',i,i2,Plm(i2)
      enddo

      do i=1,ndeg
        i2=i*(i+1)/2+1
        Plm(i2)=dsqrt(2.d0*i+1)*Plm(i2)
!         if(i2.eq.2)print*,'2 Plm(2)',Plm(2)
!        print*,'legendre 2: ',i,i2,Plm(i2),earthrad
      enddo

      AA=1.d0-X*X
      IF(AA.eq.0) then
        DO 200 l=0,ndeg
         l0=l*(l+1)/2+1
         do 200 m=1,l
!         if(m.eq.2)print*,'Plm(2)',Plm(2)
200       Plm(m+l0)=0
        RETURN    
      endif
!      print*,'still in legendre'
      AS=dSQRT(AA)
      Plm(3)=AS*dSQRT(1.5D0)
      Plm(5)=3*AS*X*dSQRT(2.5D0/3.d0)
      Plm(6)=3*AA*dSQRT(5.D0/24.d0)

      DO l=3,ndeg
        SUM=1.d0
        DO  J=1,l
          S=float(2*J-1)/(2.d0*J)
          SUM=SUM*dSQRT(AA*S) 
        enddo
        l0=l*(l+1)/2+1
        Plm(l+l0)=SUM*dSQRT(2.D0*l+1.d0)
!         if(l+l0.eq.2)print*,'3 Plm(2)',Plm(2)
!        print*,'legendre 3: ',l,l+l0,Plm(l+l0),earthrad
        SUM=1.d0
        DO  J=2,l
          S=float(2*J-1)/(2*J-2.d0)
          SUM=SUM*dSQRT(AA*S)
        enddo
        plm(l-1+l0)=X*dSQRT(2.D0*l+1.d0)*SUM
!         if(l-1+l0.eq.2)print*,'Plm(2)',Plm(2)
!        print*,'legendre 4: ',l,Plm(l-1+l0),earthrad
        DO  M=l-2,1,-1
          C=2.d0*(M+1)/dSQRT((l+M+1.D0)*(l-M))
          C=C*X/AS
          D=(l+M+2)*(l-M-1)/((l+M+1.D0)*(l-M))
          D=dSQRT(D)
          Plm(M+l0)=C*Plm(M+1+l0)-D*Plm(M+2+l0)  
         if(M+l0.eq.2)print*,'Plm(2)',Plm(2)
!        print*,'legendre 5: ',l,m,m+l0,Plm(m+l0),earthrad
!        if(m.eq.19.and.l.eq.34)then
!          print*,C,Plm(M+1+l0),D,Plm(M+2+l0),l0
!          stop
!        endif
        enddo
      enddo

!          print*,'end of legendre Plm(2)',Plm(2)
      end
























