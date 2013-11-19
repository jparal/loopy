
from __future__ import division, print_function, absolute_import
import numpy as np
from numpy import ix_, sin, cos
from scipy.interpolate import griddata
from pyhdf.SD import SD, SDC
from loopy.types import struct
from loopy.io import loadhdf4

__all__ = ['MHDFile', 'read_mhd_var', 'read_mhd_efield', 'read_mhd_grid']

class MHDFile(object):
    def __init__(self, fname):
        self._fname = fname
        self._bfield = None
        self._efield = None
        self._velocity = None
        self._plasma = None
        self._grid = None


    def cart2sph(self, vec):
        """Converts a vector from Cartesian to spherical coordinates useing internal
           grid object."""
        tht = self.grid.tht
        phi = self.grid.phi

        sph = struct()
        sph.rad = vec.x * sin(tht) * cos(phi) + \
          vec.y * sin(tht) * sin(phi) + vec.z * cos(tht)
        sph.tht = vec.x * cos(tht) * cos(phi) + \
          vec.y * cos(tht) * sin(phi) - vec.z * sin(tht)
        sph.phi = -vec.x * sin(phi) + vec.y * cos(tht)

        return sph

    def interp(self, var, xyzmesh, method='linear'):
        x = self.grid.cx.flatten()
        y = self.grid.cy.flatten()
        z = self.grid.cz.flatten()

        points = np.vstack( [x, y, z] ).transpose()
        values = var.flatten()
        xi, yi, zi = xyzmesh
        loc = np.array([xi.flatten(),yi.flatten(),zi.flatten()]).transpose()
        dat = griddata(points, values, loc, method=method)

        return dat.reshape(xi.shape).astype(var.dtype)

    @property
    def B(self):
        """Get magnetic field (nT)"""
        return self._get_bfield()

    @property
    def E(self):
        """Get electric field (mV/m)"""
        return self._get_efield()

    @property
    def V(self):
        """Get bulk velocity of plasma (km/s)"""
        return self._get_velocity()

    @property
    def dn(self):
        """Get numer density (#/cm^2)"""
        return self._get_plasma().dn

    @property
    def cs(self):
        """Get sound speed (km/s)"""
        return self._get_plasma().cs

    @property
    def pt(self):
        """Get total presure (keV/cm^3)"""
        return self._get_plasma().pt

    @property
    def grid(self):
        """Get total presure (keV/cm^3)"""
        return self._get_grid()

    @property
    def fname(self):
        """Get total presure (keV/cm^3)"""
        return self._fname

    @staticmethod
    def get_fname_format():
        return '{prefix}_mhd_{date}.hdf'


    def _get_bfield(self):
        if not self._bfield:
            b = struct()
            b.x = read_mhd_var(self._fname, 'bx_', True) * 1.e5 # (nT)
            b.y = read_mhd_var(self._fname, 'by_', True) * 1.e5 # (nT)
            b.z = read_mhd_var(self._fname, 'bz_', True) * 1.e5 # (nT)
            self._bfield = b

        return self._bfield

    def _get_efield(self):
        if not self._efield:
            self._efield = read_mhd_efield(self._fname, self._get_grid())

        return self._efield

    def _get_velocity(self):
        if not self._velocity:
            v = struct()
            v.x = read_mhd_var(self._fname, 'vx_', True) * 1e-5 # (km/s)
            v.y = read_mhd_var(self._fname, 'vy_', True) * 1e-5 # (km/s)
            v.z = read_mhd_var(self._fname, 'vz_', True) * 1e-5 # (km/s)
            self._velocity = v

        return self._velocity

    def _get_plasma(self):
        if not self._plasma:
            dn = read_mhd_var(self._fname, 'rho_', True)
            cs = read_mhd_var(self._fname, 'c_', True)
            self._plasma = struct()
            self._plasma.dn = dn*4.7619e23 # (#/cm^3)
            self._plasma.cs = cs*1e-5      # (km/s)
            self._plasma.pt = dn*cs**2*3.75e8 # Total pressure (keV/cm^3)

        return self._plasma

    def _get_grid(self):
        if not self._grid:
            self._grid = read_mhd_grid(self._fname)

        return self._grid



############################################################
#   Help functions
############################################################

def read_mhd_var(fname, varname, close_grid = False):
    """ Read a single variable from LFM HDF files and permute indecies
    to follow LFM/Fortran codes"""
    data = loadhdf4(fname,varname)
    data = data.transpose((2,1,0))

    if close_grid:
        data = _close_grid(data, True)

    return data


def read_mhd_efield(fname, grid=None):
    if not grid:
        grid = read_mhd_grid(fname)

    reinv = 6.38e8
    x = grid.x * reinv # (convert back to cm)
    y = grid.y * reinv
    z = grid.z * reinv
    ei = read_mhd_var(fname, 'ei_', False)
    ej = read_mhd_var(fname, 'ej_', False)
    ek = read_mhd_var(fname, 'ek_', False)
    dtype = ei.dtype

    [nip1,njp1,nkp1] = ei.shape
    ni = nip1-1
    nj = njp1-1
    nk = nkp1-1;
    njp2 = njp1+1
    ss = [ni,njp2,nkp1]

    i = np.arange(ni)
    j = np.arange(nj)
    k = np.arange(nk)

    outdim = [ni,nj,nk,3]
    onedim = [ni,nj,nk,1]

    et = np.zeros(outdim, dtype)
    xt = np.zeros(outdim, dtype)
    yt = np.zeros(outdim, dtype)
    zt = np.zeros(outdim, dtype)

    i0j0k0 = ix_(i  ,j  ,k  )
    i0j0k1 = ix_(i  ,j  ,k+1)
    i0j1k0 = ix_(i  ,j+1,k  )
    i0j1k1 = ix_(i  ,j+1,k+1)
    i1j0k0 = ix_(i+1,j  ,k  )
    i1j0k1 = ix_(i+1,j  ,k+1)
    i1j1k0 = ix_(i+1,j+1,k  )
    i1j1k1 = ix_(i+1,j+1,k+1)

    et[ix_(i,j,k,[0])] = 0.25 * \
      (ei[i0j0k0]+ei[i0j0k1]+ei[i0j1k0]+ei[i0j1k1]).reshape(onedim)
    et[ix_(i,j,k,[1])] = 0.25 * \
      (ej[i0j0k0]+ej[i0j0k1]+ej[i1j0k0]+ej[i1j0k1]).reshape(onedim)
    et[ix_(i,j,k,[2])] = 0.25 * \
      (ek[i0j0k0]+ek[i1j0k0]+ek[i0j1k0]+ek[i1j1k0]).reshape(onedim)

    xt[ix_(i,j,k,[0])] = 0.25 * \
      ( (x[i1j0k0]+x[i1j1k0]+x[i1j0k1]+x[i1j1k1])-\
        (x[i0j0k0]+x[i0j1k0]+x[i0j0k1]+x[i0j1k1]) ).reshape(onedim)
    yt[ix_(i,j,k,[0])] = 0.25 * \
      ( (y[i1j0k0]+y[i1j1k0]+y[i1j0k1]+y[i1j1k1])-\
        (y[i0j0k0]+y[i0j1k0]+y[i0j0k1]+y[i0j1k1]) ).reshape(onedim)
    zt[ix_(i,j,k,[0])] = 0.25 * \
      ( (z[i1j0k0]+z[i1j1k0]+z[i1j0k1]+z[i1j1k1])-\
        (z[i0j0k0]+z[i0j1k0]+z[i0j0k1]+z[i0j1k1]) ).reshape(onedim)

    xt[ix_(i,j,k,[1])] = 0.25 * \
      ( (x[i0j1k0]+x[i1j1k0]+x[i0j1k1]+x[i1j1k1])-\
        (x[i0j0k0]+x[i1j0k0]+x[i0j0k1]+x[i1j0k1]) ).reshape(onedim)
    yt[ix_(i,j,k,[1])] = 0.25 * \
      ( (y[i0j1k0]+y[i1j1k0]+y[i0j1k1]+y[i1j1k1])-\
        (y[i0j0k0]+y[i1j0k0]+y[i0j0k1]+y[i1j0k1]) ).reshape(onedim)
    zt[ix_(i,j,k,[1])] = 0.25 * \
      ( (z[i0j1k0]+z[i1j1k0]+z[i0j1k1]+z[i1j1k1])-\
        (z[i0j0k0]+z[i1j0k0]+z[i0j0k1]+z[i1j0k1]) ).reshape(onedim)

    xt[ix_(i,j,k,[2])] = 0.25 * \
      ( (x[i0j0k1]+x[i1j0k1]+x[i0j1k1]+x[i1j1k1])-\
        (x[i0j0k0]+x[i0j1k0]+x[i1j0k0]+x[i1j1k0]) ).reshape(onedim)
    yt[ix_(i,j,k,[2])] = 0.25 * \
      ( (y[i0j0k1]+y[i1j0k1]+y[i0j1k1]+y[i1j1k1])-\
        (y[i0j0k0]+y[i0j1k0]+y[i1j0k0]+y[i1j1k0]) ).reshape(onedim)
    zt[ix_(i,j,k,[2])] = 0.25 * \
      ( (z[i0j0k1]+z[i1j0k1]+z[i0j1k1]+z[i1j1k1])-\
        (z[i0j0k0]+z[i0j1k0]+z[i1j0k0]+z[i1j1k0]) ).reshape(onedim)

    dat = struct()
    dat.x = np.zeros(ss, dtype)
    dat.y = np.zeros(ss, dtype)
    dat.z = np.zeros(ss, dtype)

    def dot(a, b):
        return np.sum(a*b, axis=a.ndim-1)

    # 1e-3 is conversion from cm/s x Gauss to mV/m
    det = 1.e-3 / dot(xt,np.cross(yt,zt))
    dat.x[ix_(i,j+1,k)] = dot(et, np.cross(yt,zt))*det
    dat.y[ix_(i,j+1,k)] = dot(xt, np.cross(et,zt))*det
    dat.z[ix_(i,j+1,k)] = dot(xt, np.cross(yt,et))*det

    dat.x = - _close_grid(dat.x, False);
    dat.y = - _close_grid(dat.y, False)
    dat.z = - _close_grid(dat.z, False)

    return dat



def read_mhd_grid(fname):

    reinv=1./6.38e8

    x = read_mhd_var(fname, 'X_grid', False)*reinv # (Re)
    y = read_mhd_var(fname, 'Y_grid', False)*reinv # (Re)
    z = read_mhd_var(fname, 'Z_grid', False)*reinv # (Re)

    [nip1,njp1,nkp1]=x.shape;
    ni=nip1-1
    nj=njp1-1
    nk=nkp1-1
    njp2=njp1+1
    ss = [ni,njp2,nkp1]

    #    CELL CENTRES
    ########################################################

    i=np.arange(ni)
    j=np.arange(nj)
    k=np.arange(nk)

    # Allocate memmory in advance
    dat = struct()
    dat.x = x
    dat.y = y
    dat.z = z
    dat.cx = np.zeros(ss, x.dtype)
    dat.cy = np.zeros(ss, x.dtype)
    dat.cz = np.zeros(ss, x.dtype)

    # X,Y,Z are already in Re units
    dat.cx[ix_(i,j+1,k)] = 0.125*( \
    x[ix_(i  ,j  ,k  )] + x[ix_(i+1,j  ,k  )] + x[ix_(i  ,j+1,k  )] + \
    x[ix_(i+1,j+1,k  )] + x[ix_(i  ,j  ,k+1)] + x[ix_(i  ,j+1,k+1)] + \
    x[ix_(i+1,j  ,k+1)] + x[ix_(i+1,j+1,k+1)] )

    dat.cy[ix_(i,j+1,k)] = 0.125*( \
    y[ix_(i  ,j  ,k  )] + y[ix_(i+1,j  ,k  )] + y[ix_(i  ,j+1,k  )] + \
    y[ix_(i+1,j+1,k  )] + y[ix_(i  ,j  ,k+1)] + y[ix_(i  ,j+1,k+1)] + \
    y[ix_(i+1,j  ,k+1)] + y[ix_(i+1,j+1,k+1)] )

    dat.cz[ix_(i,j+1,k)] = 0.125*( \
    z[ix_(i  ,j  ,k  )] + z[ix_(i+1,j  ,k  )] + z[ix_(i  ,j+1,k  )] + \
    z[ix_(i+1,j+1,k  )] + z[ix_(i  ,j  ,k+1)] + z[ix_(i  ,j+1,k+1)] + \
    z[ix_(i+1,j  ,k+1)] + z[ix_(i+1,j+1,k+1)] )

    dat.cx = _close_grid(dat.cx, False)
    dat.cy = _close_grid(dat.cy, False)
    dat.cz = _close_grid(dat.cz, False)

    # Spherical coordinates of grid centers
    dat.rad = np.sqrt(dat.cx**2 + dat.cy**2 + dat.cz**2) # Radius
    dat.tht = np.arccos(dat.cz / dat.rad)
    dat.phi = np.arctan2(dat.cx, dat.cy)

    return dat


def _close_grid(var, isorigvar = True):
    """
This function assumes that variable is of size: [nip1,njp1,nkp1]
If "isorigvar" is true (default) [ni,njp2,nkp1]
If "isorigvar" is false, the function creates grid of size [ni,njp2,nkp1] which
is same as the grid definition in hdftake-mpi.F
If "isorigvar" is false (needed in lfm_rdhdf_grid) we expect that grid
is already of correct size [ni,njp2,nkp1] and we only close the values.
"""

    if isorigvar:
        [nip1,njp1,nkp1] = var.shape
        ni=nip1-1
        nj=njp1-1
        nk=nkp1-1
        njp2=njp1+1
        ss = [ni,njp2,nkp1]

        dat = np.zeros(ss, dtype=var.dtype)
        i=np.arange(ni)
        j=np.arange(nj)
        k=np.arange(nk)
        dat[ix_(i,j+1,k)] = var[ix_(i,j,k)]
    else:
        [ni,njp2,nkp1]=var.shape;
        nk=nkp1-1
        nj=njp2-2
        njp1=njp2-1
        dat = var

    i=np.arange(ni)
    j=np.arange(0,njp2,njp1)
    k=np.arange(nk)

    jj0 = np.maximum(1, np.minimum(nj, j))
    var_av = np.sum(dat[ix_(i,jj0,k)], 2)/nk
    dat[ix_(i,j,k)] = np.tile(var_av[:,:,np.newaxis], [1,1,nk])

    i=np.arange(ni)
    j=np.arange(njp2)
    dat[ix_(i,j,[nk])] = dat[ix_(i,j,[0])]

    return dat
