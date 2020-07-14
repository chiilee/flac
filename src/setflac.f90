
! Setup some parameters (rmass,amass,initial stress,vel,viscosity)

subroutine setflac
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
! thickness of new crust
integer:: j
real*8, parameter :: new_crust_thickness = 7.e3


nloop = 0
time = 0.
! Mesh generator
call init_cord

! Initial accumulated plastic strain
aps = 0

! Initial velocity
vel = 0

dvol = 0
strain = 0

! Phases in the mesh
call init_phase

! Setup markers
if (iint_marker.eq.1) then
call init_marker
endif
! Setup tracers
if (iint_tracer.eq.1) call init_tracer

! Inverse Areas of triangles
call init_areas

! Initiate temperature field
call init_temp

! Calculation of the initial STRESSES (as hydrostatic)
call init_stress

! Setup boundary conditions
call init_bc

! Distribution of REAL masses to nodes
call rmasses

! Initialization of viscosity
if( ivis_present.eq.1 ) call init_visc

! Inertial masses and time steps (elastic and maxwell)
call dt_mass
dt = min( dt_elastic, dt_maxwell )
time_t = 0

! Search the element for melting
do j = 1, nz-1
   ! search for crustal depth
   dep = 0.25*(cord(j,1,2)+cord(j+1,1,2)+cord(j,2,2)+cord(j+1,2,2))
   if (cord(1,1,2) - dep >= new_crust_thickness) exit
end do
new_crust_thickness_index = j
open(12,file='newOC.0')
write(12,*) new_crust_thickness_index
close(12)


return
end
