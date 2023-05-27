

PROGRAM LSR

!updated on 4 May 2011

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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!time series plot according to 5 points of Chen 2009 paper

!A 240 -74 43276
!B 300 -66 54144
!C 300 -74 54136
!D 60E -66 10704
!E 350 -74 63186



open(unit=19, file="time_E.DAT")



DO j=1,293

write (19, 119) days(j)/1000,cnm(j,63186)

119 format(f12.6,1x,f12.6)


END DO
 close(19)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)/60


 end PROGRAM LSR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!















