subroutine init_marker

USE marker_data

use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

double precision :: a(3,2), b(3,2), points(9,2)
double precision, parameter :: half = 0.5d0
double precision, parameter :: onesixth = 0.1666666666666666666666d0
double precision, parameter :: fivesixth = 0.8333333333333333333333d0

nphase_counter = 0
ntopmarker(:) = 0
itopmarker(:,:) = 0

! define euler coordinate of the markers
! Distribute evenly first then randomize the distribution
! to start 9 markers per elements
nmarkers = 0

! zones with 9 markers per elements
! calculate the id (element number) of the zones of high res

!call random_seed
!write(333,*) 'Call to random_seed(), result may be stochastic'

do i = 1 , nx-1
    do j = 1 , nz-1
        ! Alog the edge of an element; a and b are the nodes
        !   a - o -- x -- v - b
        !
        ! x is located at the midpint
        !   x = a / 2 + b / 2
        !
        ! o is located at a distance of 1/6 length
        !   o = 5/6 * a + 1/6 * b
        !
        ! v is located at a distance of 5/6 length
        !   o = 1/6 * a + 5/6 * b
        !
        ! Considering two elements
        !   a - o1 -- x1 -- v1 - b - o2 -- x2 -- v2 - c
        !
        ! o1, x1, v1, o2, x2, v2 will be equi-distant
        !

        a(1,:) = cord(j,i,:)*fivesixth + cord(j+1,i,:)*onesixth
        a(2,:) = cord(j,i,:)*half + cord(j+1,i,:)*half
        a(3,:) = cord(j,i,:)*onesixth + cord(j+1,i,:)*fivesixth

        b(1,:) = cord(j,i+1,:)*fivesixth + cord(j+1,i+1,:)*onesixth
        b(2,:) = cord(j,i+1,:)*half + cord(j+1,i+1,:)*half
        b(3,:) = cord(j,i+1,:)*onesixth + cord(j+1,i+1,:)*fivesixth

        points(1,:) = a(1,:)*fivesixth + b(1,:)*onesixth
        points(2,:) = a(1,:)*half + b(1,:)*half
        points(3,:) = a(1,:)*onesixth + b(1,:)*fivesixth

        points(4,:) = a(2,:)*fivesixth + b(2,:)*onesixth
        points(5,:) = a(2,:)*half + b(2,:)*half
        points(6,:) = a(2,:)*onesixth + b(2,:)*fivesixth

        points(7,:) = a(3,:)*fivesixth + b(3,:)*onesixth
        points(8,:) = a(3,:)*half + b(3,:)*half
        points(9,:) = a(3,:)*onesixth + b(3,:)*fivesixth

        dx = cord(j,i+1,1) - cord(j,i,1)
        dy = cord(j+1,i,2) - cord(j,i,2)

! randomize the new coordinates inside the element
        l = 1
        do while (l .le. 9)
            ! position of the marker
            call random_number(rx)
            call random_number(ry)
            rx = 0.5 - rx
            ry = 0.5 - ry
            ddx = dx*rx/3
            ddy = dy*ry/3
            xx = points(l,1) + ddx
            yy = points(l,2) + ddy

            ! phase of the marker
            ! smooth transition of marker phase
            do n = 1, nzone_age
                if (i<ixtb1(n) .or. i>ixtb2(n)) cycle
                if (n /= 1) then
                    if (iph_col_trans(n-1) == 1) cycle
                endif
                ycol1 = hc1(n)
                ycol2 = hc2(n)
                ycol3 = hc3(n)
                ycol4 = hc4(n)
                if (iph_col_trans(n) == 1) then
                    i1 = ixtb1(n)
                    i2 = ixtb2(n)
                    r = (cord(1,i,1) - cord(1,i1,1)) / (cord(1,i2,1) - cord(1,i1,1))
                    ycol1 = hc1(n) + (hc1(n+1) - hc1(n)) * r
                    ycol2 = hc2(n) + (hc2(n+1) - hc2(n)) * r
                    ycol3 = hc3(n) + (hc3(n+1) - hc3(n)) * r
                    ycol4 = hc4(n) + (hc4(n+1) - hc4(n)) * r
                endif

                ! layer
                yyy = yy * (-1e-3)

                if (yyy.lt.ycol1) then
                    kph = iph_col1(n)
                else if (yyy.lt.ycol2) then
                    kph = iph_col2(n)
                else if (yyy.lt.ycol3) then
                    kph = iph_col3(n)
                else if (yyy.lt.ycol4) then
                    kph = iph_col4(n)
                else
                    kph = iph_col5(n)
                end if
                exit
            end do

            call add_marker(xx, yy, kph, 0.d0, nmarkers, j, i, inc)
            !call add_marker(xx, yy, iphase(j,i), 0.d0, nmarkers, j, i, inc)
            if(inc.eq.0) cycle

            l = l + 1
            !print *, xx, yy, mark(kk)%a1, mark(kk)%a2, mark(kk)%ntriag
        enddo
    enddo
enddo

!   Put initial heterogeneities
do i = 1,inhom
    if (inphase(i) < 0) cycle  ! skip

    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        call newphase2marker(iy1(i),iy2(i),ix1(i),ix2(i),inphase(i))
    endif

    ! Gauss shape:
    if (igeom(i).eq.1.or.igeom(i).eq.2) then
        ! symmetric case:
        if (igeom(i).eq.1) then
            ixc  = (ix1(i)+ix2(i))/2
            iwidth = (ix2(i)-ix1(i))
        else
            ixc    = ix1(i)
            iwidth = ix2(i)-ix1(i)
        endif

        iamp = iy2(i)-iy1(i)

        do j = ix1(i),ix2(i)
            itop = itop_geom(j,ixc,iwidth,iamp)
            do k = iy1(i),iy2(i)
                if (k .ge. (iy2(i)-itop)) then
                    call newphase2marker(k,k,j,j,inphase(i))
                endif
            end do
        end do
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = nint(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            call newphase2marker(k,k,j,j,inphase(i))
        end do
    endif

    ! Weak zone in accumulated plastic strain at 45 degree
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k1 = floor(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            k2 = ceiling(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            if( inphase(i) .ge. 0) then
                call newphase2marker(k1,k2,j,j,inphase(i))
            endif
        end do
    endif
end do

write(333,*) '# of markers', nmarkers

call marker2elem

return
end subroutine init_marker


subroutine add_marker(x, y, iph, age, kk, j, i, inc)
  ! Add a marker at physical coordinate (x, y), with phase iph and age, to
  ! element (j, i). The current (before adding thsi marker) marker size
  ! is kk. If (x, y) is not within the element, inc is set to 0 and
  ! marker not added. Otherwise, marker is added to "mark" array and kk
  ! incremented by 1.
  !
  ! *** This subroutine is not thread-safe. DON'T CALL IT WITHIN
  ! *** OPENMP/OMP SECTION.

  USE marker_data
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'

  character*200 msg

  call check_inside(x , y, bar1, bar2, ntr, i, j, inc)
  if(inc.eq.0) return

  if(j == 1) then
      if(ntopmarker(i) == max_markers_per_elem) then
          write(msg,*) 'Too many markers at surface elements:', i, ntopmarker(i)
          call SysMsg(msg)
          call SysMsg('Marker skipped, not added!')
          return
      endif
      ! recording the id of markers belonging to surface elements
      ntopmarker(i) = ntopmarker(i) + 1
      itopmarker(ntopmarker(i), i) = kk + 1
  end if

  kk = kk + 1

  mark(kk)%x = x
  mark(kk)%y = y
  mark(kk)%dead = 1
  mark(kk)%ID = kk
  mark(kk)%a1 = bar1
  mark(kk)%a2 = bar2
  mark(kk)%age = age
  mark(kk)%ntriag = ntr
  mark(kk)%phase = iph

  nphase_counter(iph,j,i) = nphase_counter(iph,j,i) + 1

  if(kk > max_markers) then
      call SysMsg('ADD_MARKER: # of markers exceeds max. value. Please increase mark array size.')
      stop 15
  endif

end subroutine add_marker
