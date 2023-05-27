program rate4onepoint

!       updated on 28 may 2012, we estimted erros from error coeeficients
!       for point in polygon
!       epoch was corrected
!       on mac:calling 293 EWT files when spherical represenation of 120 regions has been don by mass_loss_oz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This program will (1) create a set of unevenly spaced
!       data points from a spherical harmonic file, and then
!       (2) expand these into spherical harmonics using a least-squares
!       inversion if there are more data points than spherical harmonic
!       coefficients (overdetermined case), or using a minimum norm
!       solution if there are more coefficients than data points
!       (underdetermined case).
!
!       Written by Neda Darbeheshti (June 2004)
!
!       Copyright (c) 2005, Mark A. Wieczorek
!       All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       

        implicit none
        integer, parameter ::   dmax = 65160
        real*8 ::               d(dmax), lat(dmax), lon(dmax),obs(dmax)
                                      
        integer ::              lmax, l, i,j,jj,kk,ii
        character*120 ::                infile,output_file

      
       character(len=45), dimension(:), allocatable :: names, namesg,namesn
     
       CHARACTER(len=32) :: arg

! Variables to extract the epoch from the filename
          integer year(293), idoy(293)
      real*8 epoch(293)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables for subroutine date_and_time

      integer time_array_0(8), time_array_1(8)
      real start_time, finish_time



! Mark the beginning of the program

      call date_and_time(values=time_array_0)
      start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate (namesg(293))


OPEN(10,file='GRACE_EWT_names.DAT',form='formatted')

DO j=1,293

read (10, 100) namesg(j)

100 format(a45)

! extract the epoch from the filename
           read(namesg(j)(7:10),'(i4)')year(j)
           read(namesg(j)(11:13),'(i3)')idoy(j)
           epoch(j) = float(year(j)) + float(idoy(j)+5)/365.d0
!           print*,'epoch = ',epoch(j)
       
END DO

 close(10)

! stop

       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!opening and reading 293 files

OPEN(20,file='EWT315_064',form='formatted')


DO jj=1,293


OPEN(jj+500,file=namesg(jj),form='formatted',status='old')

DO i=1,65160


read (jj+500, 101) lon(i),lat(i),obs(i),d(i)


101 format(f5.1,1x,f5.1,1x,d19.12,1x,d19.12)

 END DO

 close(jj+500)


write (20, *) epoch(jj),obs(56989),d(56989)

 END DO

 close(20)


 call date_and_time(values=time_array_1)
      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
           + time_array_1 (7) + 0.001 * time_array_1 (8)

write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)/60




end program rate4onepoint

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


                                                              
                                                                  


                                                                                           
