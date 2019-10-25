! -*- F90 -*-

!!! maximum number of ELEMENTS !!!
parameter(mnz=300, mnx=700)
parameter(max_markers_per_elem=32)

common /phase_arrays/ iphase(mnz,mnx), nphase_counter(20,mnz,mnx), phase_ratio(20,mnz,mnx)

common /topo_arrays/ dtopo(mnx+1), dhacc(mnx+1), extrusion(mnx), &
     andesitic_melt_vol(mnx), extr_acc(mnx), &
     ntopmarker(mnx), itopmarker(max_markers_per_elem,mnx)

common /all_arrays/ strainr(3,4,mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx),irheol_fl(mnz,mnx)

common /temp_arrays/ temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bcstress((mnz+mnx)*2,3), &
     nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

