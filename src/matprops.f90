!==============================================
! Density
function Eff_dens( j, i)
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  include 'phases.inc'

  zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

  ! 660-km phase change with clapeyron slope 600C/50km
  if (zcord < -660e3 + (tmpr - t_bot) * 500./6) then
      Eff_dens = 3450 * (1 - 2e-5 * tmpr)  ! ~6% density jump
      return
  endif


  press = 0
  do ii = 1, 4
      press = press - (stress0(j,i,1,ii)+stress0(j,i,2,ii)+stress0(j,i,4,ii))
  enddo
  press = press / 12

  if (iint_marker.ne.1) then
      iph = iphase(j,i)
      dens = den(iph) * ( 1 - alfa(iph)*tmpr + beta(iph)*press )

!!$      ! Effect of melt
!!$      fmelt = Eff_melt(iph, tmpr)
!!$      dens = dens * ( 1.-0.1*fmelt )
!!$      Eff_dens = dens
  else
      Eff_dens = 0.
      do k = 1, nphase
          ratio = phase_ratio(k,j,i)
          ! when ratio is small, it won't affect the density
          if(ratio .lt. 0.01) cycle

          dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

          !press = 0
          !press = stressI(j,i)
          press = dens*g*zcord

          dens = den(k) * ( 1 - alfa(k)*tmpr + beta(k)*press )

          if(k==ksed1 .or. k==ksed2) then
              sed_min_density = 2400.
              delta_den = den(k) - sed_min_density
              zefold = 6000.
              if (j==1 .and. cord(j,i,2)>0.) then
                  ! sediment above ground is already compacted
                  dens = den(k) * (1.-alfa(k)*tmpr)
              else
                  ! sediment below sea level is uncompacted and less dense
                  dens = (den(k) - delta_den*exp((zcord-0.5*(cord(1,i,2)+cord(1,i+1,2)))/zefold)) &
                       * ( 1 - alfa(k)*tmpr + beta(k)*press )
              endif
              if (dens < sed_min_density) dens = sed_min_density
          endif

          ! Effect of melt
          !fmelt = Eff_melt(k, tmpr)
          !dens = dens * ( 1.-0.1*fmelt )

          Eff_dens = Eff_dens + ratio*dens

      enddo
  endif
  return
end function Eff_dens


!!$!==============================================
!!$! Melt fraction
!!$!==============================================
!!$function Eff_melt(iph, tmpr)
!!$include 'precision.inc'
!!$include 'params.inc'
!!$
!!$if( tmpr .lt. ts(iph) ) then
!!$    ! below solidus
!!$    fm = 0.
!!$elseif( tmpr .lt. tk(iph) ) then
!!$    fm = fk(iph)/(tk(iph)-ts(iph)) * (tmpr-ts(iph))
!!$    fm = min(max(fm, 0.), 1.)
!!$elseif( tmpr .lt. tl(iph) ) then
!!$    fm = (1.-fk(iph))/(tl(iph)-tk(iph))*(tmpr-tk(iph)) + fk(iph)
!!$    fm = min(max(fm, 0.), 1.)
!!$else
!!$    fm = 1.
!!$endif
!!$
!!$Eff_melt = fm
!!$
!!$return
!!$end function Eff_melt
!!$

!=================================================
! Effective Heat Capacity incorporating latent heat
!=================================================
function Eff_cp( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


HeatLatent = 420000.

iph = iphase(j,i)
Eff_cp = cp(iph)


!!$tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!!$if( tmpr .lt. ts(iph) ) then
!!$    Eff_cp = cp(iph)
!!$elseif( tmpr .lt. tk(iph) ) then
!!$    Eff_cp = cp(iph) + HeatLatent * fk(iph)/(tk(iph)-ts(iph))
!!$elseif( tmpr .lt. tl(iph) ) then
!!$    Eff_cp = cp(iph) + HeatLatent * (1.-fk(iph))/(tl(iph)-tk(iph))
!!$else
!!$    Eff_cp = cp(iph)
!!$endif
!!$
!!$
!!$! HOOK
!!$! Intrusions - melting effect - see user_ab.f90
!!$if( if_intrus .eq. 1 ) then
!!$    HeatLatent = 420000.
!!$
!!$    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
!!$    if( tmpr .lt. ts(iph) ) then
!!$        Eff_cp = cp(iph)
!!$    elseif( tmpr .lt. tl(iph)+1 ) then
!!$        Eff_cp = cp(iph) + HeatLatent/(tl(iph)-ts(iph))
!!$    else
!!$        Eff_cp = cp(iph)
!!$    endif
!!$endif

return
end function Eff_cp


!=================================================
! Effective Thermal Conductivity
!=================================================
function Eff_conduct( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


if (iint_marker.ne.1) then
    iph = iphase(j,i)
    cond = conduct(iph)

    !if( den(iph) .lt. 3000. ) then  ! for crustal material
    !    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
    !    if( tmpr.lt.25 ) tmpr = 25.
    !    Eff_conduct = -0.38*dlog(tmpr) + 4.06
    !endif

    ! HOOK
    ! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
    if( if_hydro .eq. 1 ) then
        cond = HydroDiff(j,i)*den(iph)*cp(iph)
    endif
    Eff_conduct = cond

else
    Eff_conduct = 0.
    do k = 1 , nphase

        ! when ratio is small, it won't affect the density
        if(phase_ratio(k,j,i) .lt. 0.01) cycle

        cond = conduct(k)

        !if( den(k) .lt. 3000. ) then  ! for crustal material
        !    tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
        !    if( tmpr.lt.25 ) tmpr = 25.
        !    Eff_conduct = -0.38*dlog(tmpr) + 4.06
        !endif

        ! HOOK
        ! Hydrothermal alteration of thermal diffusivity  - see user_luc.f90
        if( if_hydro .eq. 1 ) then
            cond = HydroDiff(j,i)*den(k)*cp(k)
        endif
        Eff_conduct = Eff_conduct + phase_ratio(k,j,i)*cond
    enddo
endif

!write(*,*) Eff_conduct, cond
return
end function Eff_conduct



!=================================================
! Non-Newtonian viscosity
!=================================================

! from Chen and Morgan (1990)
! Uses A value in MPa and but gives viscosity in (Pa * s) 
! Therefore there is a coefficient 1.e+6 

function Eff_visc( j, i )
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
include 'phases.inc'

Eff_visc = 0.
r=8.31448e0
tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))

zcord = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
if (zcord < -660e3) then
    ! T-dep. viscosity
    ! eta=1e22 when T=1330 Celsius, inc. 10x when T is 200 deg lower
    Eff_visc = 1e22 * exp( -0.0110 * (tmpr - 1330))
    vis = min(v_max, max(v_min, vis))
    return
elseif (zcord < -410e3) then
    Eff_visc = 1e20 * exp( -0.0110 * (tmpr - 1330))
    vis = min(v_max, max(v_min, vis))
    return
elseif (i .ge.nx-4) then
    ! left = MOR (chii)
    Eff_visc = 1e20
    vis = min(v_max, max(v_min, vis))
    return
endif

if (zcord < -200e3) then
    chamber(j,i)=0
endif

srat = e2sr(j,i)
if( srat .eq. 0 ) srat = vbc/rxbo

if (iint_marker.ne.1) then
    iph = iphase(j,i)

    pow  =  1./pln(iph) - 1.
    pow1 = -1./pln(iph)

    vis = 0.25 * srat**pow*(0.75*acoef(iph))**pow1* &
         exp(eactiv(iph)/(pln(iph)*r*(tmpr+273.)))*1.e+6

!!$    ! Effect of melt
!!$    fmelt_crit = 0.05
!!$    fmelt = Eff_melt(iph, tmpr)
!!$    if( fmelt .gt. 0. ) then
!!$        if( fmelt .lt. fmelt_crit ) then
!!$            vislog = fmelt/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
!!$            vis = 10.**vislog
!!$        else
!!$            vis = v_min
!!$        endif
!!$    endif


    ! limiting from above (quasi-Peierls)
    !sIImax = 5.e+8
    !vis_peierls = sIImax / srat / 2
    !if( vis .gt. vis_peierls ) vis = vis_peierls

    ! Final cut-off (max)
    if (vis .gt. v_max) vis = v_max

    !Effect of melt (chii)

    weak = weaken_max
    if (chamber(j,i) > weaken_launch) then
           if (weaken_saturated > chamber(j,i)) then
                vis_down1 = log(weak)*(chamber(j,i)-weaken_launch)/(weaken_saturated-weaken_launch)
                vis_down2 = exp(vis_down1)
           else
                vis_down2 = weak
           endif
           vis =vis * vis_down2

    endif

    ! Final cut-off (min)
    if (vis .lt. v_min) vis = v_min


    Eff_visc = vis

else

    do k = 1, nphase
        if(phase_ratio(k,j,i) .lt. 0.01) cycle

        pow  =  1./pln(k) - 1.
        pow1 = -1./pln(k)

        vis = 0.25 * srat**pow*(0.75*acoef(k))**pow1* &
             exp(eactiv(k)/(pln(k)*r*(tmpr+273.)))*1.e+6

        ! Effect of melt
        !fmelt_crit = 0.05
        !fmelt = Eff_melt(k, tmpr)
        !if( fmelt .gt. 0. ) then
        !    if( fmelt .lt. fmelt_crit ) then
        !        vislog = fmelt/fmelt_crit*dlog10(v_min/vis) + dlog10(vis)
        !        vis = 10.**vislog
        !    else
        !        vis = v_min
        !    endif
        !endif


        ! limiting from above (quasi-Peierls)
        !sIImax = 5.e+8
        !vis_peierls = sIImax / srat / 2
        !if( vis .gt. vis_peierls ) vis = vis_peierls

        !write(*,*) vis,srat,pln(k),eactiv(k)

        ! Final cut-off
        if (vis .lt. v_min) vis = v_min
        if (vis .gt. v_max) vis = v_max

        ! harmonic mean
        Eff_visc = Eff_visc + phase_ratio(k,j,i) / vis
        !write(*,*) i,j, Eff_visc, vis, tmpr,phase_ratio(k,j,i)
    enddo

    Eff_visc = 1 / Eff_visc

    !Effect of melt (chii)
    weak = weaken_max
    if (chamber(j,i) > weaken_launch) then
           if (weaken_saturated > chamber(j,i)) then
                vis_down1 = log(weak)*(chamber(j,i)-weaken_launch)/(weaken_saturated-weaken_launch)
                vis_down2 = exp(vis_down1)
           else
                vis_down2 = weak
           endif
           Eff_visc = Eff_visc * vis_down2
    endif

    ! Final cut-off (min)
    if (Eff_visc .lt. v_min) Eff_visc = v_min

endif

return
end function Eff_visc
