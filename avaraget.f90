PROGRAM avarage


!updated on 4 May 2011
!opening 293 GRGS files, computing avarage of Cnm and Snm

!use Somearrays

implicit none





double precision :: sums,sumc


INTEGER :: i,j

 CHARACTER(len=80) :: dummy




 character(len=6), dimension(:), allocatable :: chara
 character(len=45), dimension(:), allocatable :: names,namesw

integer, dimension(:), allocatable :: m,n





double precision, dimension(:), allocatable :: c,s

 double precision, dimension(:,:), allocatable :: cnm,snm
  





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

OPEN(10,file='GRACE_EWT_names.DAT',form='formatted')

DO j=1,293

read (9, 109) names(j)
read (10, 109) namesw(j)

109 format(a45)
print*,names(j)
print*,namesw(j)
END DO

 close(9)

 close(10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!number of lines of 10 day Grace GRGS coefficients file 1326




    allocate (c(1326),s(1326),n(1326),m(1326),chara(1326))
    allocate (cnm(293,1326),snm(293,1326))


 !cnm=0
 !snm=0
 !degree=0
 !order=0
 n=0
 m=0
 c=0
 s=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!opening and reading 293 files

DO j=1,293


OPEN(j+10,file=names(j),form='formatted',status='old')

do i=1,5

read(j+10,*) dummy
    
END DO



DO i=1,1326




read (j+10, 101) chara(i),n(i),m(i),c(i),s(i)

101 format(a6,i5,i5,d19.12,d19.12)

 cnm(j,i)=c(i)
 snm(j,i)=s(i)

END DO



 close(j+10)



END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!test if the reading is all right

open(unit=220, file="test.DAT")



DO i=1,1326



write (220, 230) cnm(1,i),cnm(2,i),cnm(3,i),cnm(4,i)

230 format(d19.12,d19.12,d19.12,d19.12)

END DO

 close(220)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!avarage 0f 293 Cnm Snm

 c=0
 s=0

DO i=1,1326

DO j=1,293

 c(i)=cnm(j,i)+c(i)
 s(i)=snm(j,i)+s(i)



  END DO  
END DO

 c=c/293
 s=s/293
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Writing the Cnm Snm avarages 



open(unit=200, file="avarage.DAT")



DO i=1,1326



write (200, 201) n(i),m(i),c(i),s(i)

201 format(i5,i5,d19.12,d19.12)

END DO

 close(200)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



write (6, '(8x, 1a, 1f16.3)') 'elapsed wall clock time:', &
           (finish_time - start_time)



 
end PROGRAM avarage




