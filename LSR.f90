

PROGRAM LSR
!reviewed on 22 jun 
!updated on 16 May 2011

!opening 293 EWT files, computed from GRACE_GRGS_EWH_FM.f90

implicit none


double precision :: sums,sumc


INTEGER :: i,j

 CHARACTER(len=80) :: dummy


 character(len=6), dimension(:), allocatable :: chara
 character(len=45), dimension(:), allocatable :: names

integer, dimension(:), allocatable :: m,n


double precision, dimension(:), allocatable :: landa,phi,delta_sigma,days

 double precision, dimension(:,:), allocatable :: cnm,snm,a,at,ata,iata,iataat,b,lc
  
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


OPEN(10,file='GRACE_EWT_names.DAT',form='formatted')

DO j=1,293


read (10, 100) names(j)


100 format(a45)



END DO

 close(10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!number of lines of 10 day Grace GRGS coefficients file 1326


    allocate (landa(65160),phi(65160),delta_sigma(65160))
    allocate (cnm(293,65160))
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!opening and reading 293 files

DO j=1,293


OPEN(j+300,file=names(j),form='formatted',status='old')

DO i=1,65160


read (j+300, 101) landa(i),phi(i),delta_sigma(i)


101 format(f5.1,1x,f5.1,1x,d19.12)

 cnm(j,i)=delta_sigma(i)
 

END DO



 close(j+300)



END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate (a(293,2))
allocate (at(2,293))
allocate (ata(2,2))
allocate (iata(2,2))
allocate (iataat(2,293))

allocate (b(293,1))

allocate (lc(2,1))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (days(293))


OPEN(20,file='epoch.DAT',form='formatted')

DO j=1,293


read (20, 200) days(j)



200 format(d7.0)

END DO

 close(10)

DO j=1,12
days(j)=(days(j)-2002000)/365.25
END DO

DO j=13,40
days(j)=(days(j)-2003000)/365.25+1
END DO

DO j=41,75
days(j)=(days(j)-2004000)/365.25+2
END DO

DO j=76,112
days(j)=(days(j)-2005000)/365.25+3
END DO

DO j=113,148
days(j)=(days(j)-2006000)/365.25+4
END DO

DO j=149,185
days(j)=(days(j)-2007000)/365.25+5
END DO

DO j=186,221
days(j)=(days(j)-2008000)/365.25+6
END DO

DO j=222,258
days(j)=(days(j)-2009000)/365.25+7
END DO

DO j=259,293
days(j)=(days(j)-2010000)/365.25+8
END DO


DO j=1,293
a(j,1)=j
!a(j,1)=1
a(j,2)=days(j)
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call  transp(a,at,293,2)
 call  matmult(at,a,ata,2,293,2)
 call  invert(ata,iata,2)
 call  matmult(iata,at,iataat,2,2,293)



open(unit=19, file="EWT_rate0.DAT")


DO i=1,65160


DO j=1,293

b(j,1)=cnm(j,i)

END DO


 call  matmult(iataat,b,lc,2,293,1)

write (19, 119) landa(i),phi(i),lc(1,1)


119 format(f5.1,1x,f5.1,1x,d19.12)



END DO
 close(19)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)/60



 
end PROGRAM LSR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine invert(a,b,n)

!  subroutine to invert matrix A into matrix B. Both are n x n matrices

      implicit none

      integer n,i,j,k
      real*8 a(n,n),b(n,n),z

      do  i=1,n
        do  j=1,n
          b(i,j) = 0.0
        enddo
        b(i,i) = 1.0
      enddo

      do k = 1,n
        do i = 1,n
	       if(i.ne.k)then
	         z = a(i,k)/a(k,k)
            do j=1,n
              a(i,j) = a(i,j)-a(k,j)*z
              b(i,j) = b(I,j) - b(k,j)*z
            enddo
	       endif
	     enddo
	     z = a(k,k)
	     do j=1,n
	       a(k,j) = a(k,j)/z
	       b(k,j) = b(k,j)/z
	     enddo
	   enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matmult(M1,M2,M3,l,m,n)

!  multiplies a   l x m   matrix M1 by a   m x n  matrix M2 to give a   l x n   matrix M3
!  written by Paul Tregoning
!  9th June 1992

        implicit none
        integer i,j,k,l,m,n
        real*8 M1(l,m),M2(m,n),M3(l,n),temp

        temp = 0.d0

        do i=1,l
          do k=1,n
            do  j=1,m
              temp = M2(j,k) * M1(i,j) + temp  
            enddo

            M3(i,k) = temp
            temp = 0.d0
          enddo    
        enddo   

        return
        end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transp(R,Rt,m,n)

!  transposes a mxn matrix
!  Written by Paul Tregoning
!  10/6/92

        implicit none
        integer i,j,m,n
        real*8 R(m,n),Rt(n,m)

!  set Rt to zero initially
        do i=1,n
          do j=1,m
             Rt(i,j)=0.0
          enddo  
        enddo  

!  perform the transpose
        do i=1,m
          do j=1,n
            Rt(j,i) = R(i,j)
          enddo  
        enddo   

        return
        end






