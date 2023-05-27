

PROGRAM time_point

!updated on 31 Aug 2011

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (days(293))


OPEN(20,file='epoch.DAT',form='formatted')

DO j=1,293


read (20, 200) days(j)


200 format(d7.0)

END DO

 close(10)

! Region IV  (147,-32) 26485
! Region III  (147,-43) 26474

! Region VI  (117,-32) 21055
! Region VII  (117,-25) 21062
! Region XII  (125,-25) 22510
! Region VIII  (130,-17) 23423
! Region X  (140,-25) 25225
! Region IX  (141,-19) 25412
! Region I  (147,-23) 26494
! Region II  (150,-35) 27025

! Region XI  (143,-29) 25764
! Region V  (137,-31) 24676
!open(unit=19, file="time21055.txt")
!open(unit=191, file="time21062.txt")
!open(unit=192, file="time22510.txt")
!open(unit=193, file="time23423.txt")
!open(unit=194, file="time25225.txt")
!open(unit=195, file="time25412.txt")
!open(unit=196, file="time26494.txt")
!open(unit=197, file="time27025.txt")

open(unit=196, file="time25764.txt")
open(unit=197, file="time24676.txt")
open(unit=198, file="time24676i.txt")

DO j=1,293


!write (19, 119) days(j)/1000,cnm(j,21055)
!write (191, 119) days(j)/1000,cnm(j,21062)
!write (192, 119) days(j)/1000,cnm(j,22510)
!write (193, 119) days(j)/1000,cnm(j,23423)
!write (194, 119) days(j)/1000,cnm(j,25225)
!write (195, 119) days(j)/1000,cnm(j,25412)
!write (196, 119) days(j)/1000,cnm(j,26494)
!write (197, 119) days(j)/1000,cnm(j,27025)

write (196, 119) days(j)/1000,cnm(j,25764)
write (197, 119) days(j)/1000,cnm(j,24676)

write (198, 120) j,cnm(j,24676)

119 format(f12.6,1x,f12.6)
120 format(i5,1x,f12.6)


END DO
 close(19)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)/60



 
end PROGRAM time_point

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






