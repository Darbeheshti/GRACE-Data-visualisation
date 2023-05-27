
program point_in_polygreen

!	updated on 15 Sep 2011

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	implicit none
        logical ::     point_in_polygon
        logical ::     logi
	integer, parameter ::	dmax = 65341, nb=245
        integer ::		nmax(nb),no1(nb)
	real ::		d(dmax), lat(dmax), lon(dmax) 
	real, dimension(18300) :: latr, lonr,dd
				
	integer ::		lmax, l, i,j,jj,kk
	character*80 ::		infile

       character(len=80), dimension(:), allocatable :: names,namesn
      

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables for subroutine date_and_time

      integer time_array_0(8), time_array_1(8)
      real start_time, finish_time

! Variables to extract the epoch from the filename
      integer year, idoy
      real*8 epoch

! Mark the beginning of the program

      call date_and_time(values=time_array_0)
      start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
           + time_array_0 (7) + 0.001 * time_array_0 (8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!opening the file of 1*1 of the globe

OPEN(502,file='/Users/Neda/Desktop/SHTOOLS/examples/poin_in polygon/regi1h01.txt',form='formatted',status='old')

!Users/Neda/Desktop/SHTOOLS/examples/oz_case/regi1h01_oz4indot5oz.txt
!dmax =74906
!/Users/Neda/Desktop/SHTOOLS/examples/oz_case/regi1h01_oz4indot4oz.txt
!dmax = 81963

DO i=1,dmax

read (502, *) lon(i),lat(i),d(i)

 END DO

 close(502)

d=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!devided Greenland into 16 regions (Wouters et al 2008)

allocate (names(nb))

!Wouters 16 boundry regions gridded with GMT and points outside the region set to zero
OPEN(11,file='basin_names245.txt',form='formatted')


DO j=1,nb

read (11, 100) names(j)
100 format(a45)

END DO

 close(11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!16 files to write 0&1 grid over the globe

allocate (namesn(nb))

OPEN(12,file='basinh01_names245.txt',form='formatted')


DO j=1,nb

read (12, 100) namesn(j)

END DO

 close(12)

!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(400,file='number_of.txt',form='formatted')

DO i=1,nb
read (400, *) l,nmax(i)
END DO
 close(400)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!opening and reading 16 files (regions)


do j=1,nb
      d=0
OPEN(j+300,file=names(j),form='formatted',status='old')
OPEN(j+400,file=namesn(j),form='formatted')


DO i=1,nmax(j)

read (j+300, *) lonr(i),latr(i),dd(i)

END DO


jj=0
DO i=1,dmax



  logi= point_in_polygon(lonr, latr, lon(i), lat(i), nmax(j))
if (logi ) d(i)=1

if (logi ) jj=jj+1

 END DO


print*,jj
no1(j)=jj


!writing  the file of 1*1 of the globe with regions assigned to 1
DO i=1,dmax

write (j+400, *) lon(i),lat(i),d(i)

 END DO


 close(j+400)

 close(j+300)

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(500,file='number_ofone.txt',form='formatted')

DO i=1,245
write (500, *) i,no1(i)
END DO
 close(500)

end program point_in_polygreen


! Written by Evan Gowan, last updated August 22, 2011


logical function point_in_polygon(x_boundary, y_boundary, x, y, number_points)

! this function determines whether a give point (x,y) is within a polygon defined by x_boundary and y_boundary
! number_points is the number of points in the polygon

	
	implicit none

	integer, intent(in) :: number_points
	real, intent(in) :: x, y
	real, dimension(number_points), intent(in) :: x_boundary, y_boundary

	integer :: current_point, next_point, last_point, crossover_counter
	logical :: found_first, found_last, inside

	found_first = .false.
	found_last = .false.
	inside = .false.

	current_point = 1
	search_boundary: do

		next_point = current_point + 1
	
		if (next_point > number_points) THEN
			next_point = 1
			found_last = .true.
		endif

! even-odd rule algorithm to determine if the point is inside or outside

		if (min(y_boundary(current_point), y_boundary(next_point)) < y .and.&
		    max(y_boundary(current_point), y_boundary(next_point)) >= y) THEN

			if (x_boundary(current_point) + (y - y_boundary(current_point)) /&
			    (y_boundary(next_point) - y_boundary(current_point)) * &
			    (x_boundary(next_point) - x_boundary(current_point)) < x) THEN

				inside = .not.(inside)

			endif

		endif

		current_point = current_point + 1

		if (found_last) THEN
			exit search_boundary
		endif
		
	
	end do search_boundary

	point_in_polygon = inside

	return
end function point_in_polygon








