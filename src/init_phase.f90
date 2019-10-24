!--------------------------------------------------------------
! Initialization of phases in each zone 
!--------------------------------------------------------------

subroutine init_phase
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


! Main phase
iphase = mphase

if (irphase .gt. 0) then
    !  Read distribution of the phases from the dat file
    open(12, file=phasefile)
    read(12,*) nphasl
    do k=1,nphasl
        read(12,*) lphase(k)
    enddo
    do i=1,nx-1
        do j=1,nz-1
            read(12,*) ii,jj,iphase(j,i)
        enddo
    enddo
    close(12)
else
    ! phases in horizontal layers
    do k = 1,nphasl
        do i = 1,nx-1
            do j = ltop(k),lbottom(k)
                iphase(j,i) = lphase(k)
            end do
        end do
    end do

    ! Case with iynts = 2 or 10 for continental and collision
    if (iynts.eq.2 .or. iynts.eq.10) then
        do n = 1, nzone_age
            do i = ixtb1(n), ixtb2(n)
                do j = 1, nz-1
                    y = (cord(1,i,2)-cord(j,i,2))*1.e-3
                    if (y.lt.hc1(n)) then
                        iphase(j,i) = iph_col1(n)
                    else if (y.lt.hc2(n)) then
                        iphase(j,i) = iph_col2(n)
                    else if (y.lt.hc3(n)) then
                        iphase(j,i) = iph_col3(n)
                    else if (y.lt.hc4(n)) then
                        iphase(j,i) = iph_col4(n)
                    else
                        iphase(j,i) = iph_col5(n)
                    end if
                enddo
            enddo
        enddo
    endif
endif

write(333,*) 'Phases of horizontal layers:'
do i = 1,nphasl
   write(333,*) i,lphase(i)
enddo

!   Put different rheologies for inclusions 
do i = 1,inhom
    ! Rectangular shape:
    if (igeom(i) .eq.0) then
        do j = ix1(i),ix2(i)
            do k = iy1(i),iy2(i)
                if( inphase(i) .ge. 0) then
                    iphase(k,j) = inphase(i)
                endif
                aps(k,j)=xinitaps(i)
            end do
        end do
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
                    iphase(k,j) = inphase(i)
                endif
            end do
        end do
    endif

    ! weak zone at 45 degree
    if (igeom (i) .eq.3) then
        do j = ix1(i),ix2(i)
            k = nint(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            iphase(k,j) = inphase(i)
        end do
    endif
    
    ! Weak zone in accumulated plastic strain at 45 degree        
    if (igeom (i).eq.4) then
        do j =ix1(i),ix2(i)
            k1 = floor(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            k2 = ceiling(float(iy2(i)-iy1(i))/float(ix2(i)-ix1(i))*(j-ix1(i))) + iy1(i)
            aps(k1:k2,j)=xinitaps(i)
            if( inphase(i) .ge. 0) then
                iphase(k1:k2,j) = inphase(i)
            endif
        end do
    endif
end do

! Check if viscous rheology present
ivis_present = 0
do i = 1,nx-1
    do j = 1, nz-1
        iph = iphase(j,i)
        if( irheol(iph).eq.3 .or. irheol(iph).ge.11 ) ivis_present = 1
    end do
end do
    

return
end


!==========================================================
! Gauss perturbation at the top of heterogenity 
function itop_geom(j,ixc,iwidth,iamp) 
    
    itop_geom = int(iamp*exp(-(float(j-ixc)/(0.5*float(iwidth)))**2.))

return
end 


