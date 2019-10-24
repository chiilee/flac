#!/usr/bin/env python

'''Convert the binary output of flac to VTK (vts) files.
'''

import sys, os
import zlib, base64
import numpy as np
import flac

ndims = 2
nstress = 4

# whether write VTK files in compressed binary (base64 encoded) or uncompressed plain ascii
output_in_binary = True

def main(path, start=1, end=-1):

    # changing directory
    os.chdir(path)

    fl = flac.Flac()

    nex = fl.nx - 1
    nez = fl.nz - 1
    if end == -1:
        end = fl.nrec

    for i in range(start, end+1):
        print 'Writing record #%d, model time=%.3e' % (i, fl.time[i-1])
        fvts = open('flac.%06d.vts' % i, 'w')
        vts_header(fvts, nex, nez)

        # node-based field
        fvts.write('  <PointData>\n')

        vx, vz = fl.read_vel(i)
        # VTK requires vector field (velocity, coordinate) has 3 components.
        # Allocating a 3-vector tmp array for VTK data output.
        tmp = np.zeros((fl.nx, fl.nz, 3), dtype=vx.dtype)
        tmp[:,:,0] = vx
        tmp[:,:,1] = vz
        # vts requires x-axis increment fastest, swap axes order
        vts_dataarray(fvts, tmp.swapaxes(0,1), 'Velocity', 3)

        a = fl.read_temperature(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Temperature')

        fvts.write('  </PointData>\n')

        # element-based field
        fvts.write('  <CellData>\n')

        # logrithm of strain rate 2nd invariant
        a = fl.read_srII(i)
        srat = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'Strain rate')

        a = fl.read_eII(i)
        eii = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'eII')

        a = fl.read_density(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Density')

        a = fl.read_aps(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Plastic strain')

        a = fl.read_sII(i)
        sii = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'Stress')

        sxx = fl.read_sxx(i)
        vts_dataarray(fvts, sxx.swapaxes(0,1), 'Sxx')

        szz = fl.read_szz(i)
        vts_dataarray(fvts, szz.swapaxes(0,1), 'Szz')

        sxz = fl.read_sxz(i)
        vts_dataarray(fvts, sxz.swapaxes(0,1), 'Sxz')

        pressure = fl.read_pres(i)
        vts_dataarray(fvts, pressure.swapaxes(0,1), 'Pressure')

        # compression axis of stress
        a = compute_p_axis(sxx, szz, sxz)
        vts_dataarray(fvts, a.swapaxes(0,1), 'P-axis', 3)

        a = fl.read_diss(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Dissipation')

        # logrithm of effective viscosity
        eff_visc = np.log10(a + 1e-45) + 8 - srat
        vts_dataarray(fvts, eff_visc.swapaxes(0,1), 'Eff. Visc')

        a = fl.read_visc(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Viscosity')

        a = fl.read_phase(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Phase')

        # Work done by stress
        a = sii * 1e8 * eii
        vts_dataarray(fvts, a.swapaxes(0,1), 'Work')

        fvts.write('  </CellData>\n')

        # coordinate
        x, z = fl.read_mesh(i)
        tmp[:,:,0] = x
        tmp[:,:,1] = z
        fvts.write('  <Points>\n')
        vts_dataarray(fvts, tmp.swapaxes(0,1), '', 3)
        fvts.write('  </Points>\n')


        vts_footer(fvts)
        fvts.close()
    return


def compute_p_axis(sxx, szz, sxz):
    mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
    theta = 0.5 * np.arctan2(2*sxz, sxx-szz)
    cost = np.cos(theta)
    sint = np.sin(theta)

    p_axis_x = mag * sint
    p_axis_z = mag * cost

    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    nx, nz = sxx.shape
    tmp = np.zeros((nx, nz, 3), dtype=sxx.dtype)
    tmp[:,:,0] = p_axis_x
    tmp[:,:,1] = p_axis_z
    return tmp


def vts_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (int, np.int32, np.int_):
        dtype = 'Int32'
    elif data.dtype in (float, np.single, np.double, np.float,
                        np.float32, np.float64, np.float128):
        dtype = 'Float32'
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header)

    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tostring()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header))
        f.write(base64.standard_b64encode(b))
    else:
        data.tofile(f, sep=' ')
    f.write('\n</DataArray>\n')
    return


def vts_header(f, nex, nez):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<StructuredGrid WholeExtent="0 {0} 0 {1} 0 0">
<Piece Extent="0 {0} 0 {1} 0 0">
'''.format(nex, nez))
    return


def vts_footer(f):
    f.write(
'''</Piece>
</StructuredGrid>
</VTKFile>
''')
    return


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print '''usage: flac2vtk.py path [step_min [step_max]]

Processing flac data output to VTK format.
If step_max is not given, processing to latest steps
If both step_min and step_max are not given, processing all steps'''
        sys.exit(1)

    path = sys.argv[1]

    start = 1
    end = -1
    if len(sys.argv) >= 3:
        start = int(sys.argv[2])
        if len(sys.argv) >= 4:
            end = int(sys.argv[3])

    main(path, start, end)
