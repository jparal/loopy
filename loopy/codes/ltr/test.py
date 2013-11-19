
from __future__ import division, print_function, absolute_import

from IPython.parallel import Client
from scipy.interpolate import griddata
from os.path import exists
import tables as pytab

import loopy as lpy
from loopy.codes import ltr
from loopy.io import writehdf5

######################################################################

dbeg = [2012, 10, 07]
dend = [2012, 10, 07, 00, 0]
step = 30 # (sec)
prefix = '~/data/121007_Baker13/baker'
oname = 'sample_eqplane_lr.h5'

# if exists(oname):
#     print('Loading previous state from: %s' % oname)
#     with tables.file.openFile(oname, 'r') as f:
#         dat = it

dat=lpy.types.struct()
dat.radrng = arange(4,12,0.1, dtype=float32)
dat.phirng = linspace(0, 2*pi, 101)

rad, phi, z = meshgrid(dat.radrng, dat.phirng, [0])
x = rad * cos(phi)
y = rad * sin(phi)
xyzmesh = [x,y,z]

sim = ltr.FileIter(prefix, dbeg, dend, step)

def interpolate_mhd(mhd):
    import loopy as lpy
    dat = lpy.struct()
    esph = mhd.cart2sph(mhd.E)
    dat.Bz = mhd.interp(mhd.B.z, xyzmesh).squeeze()
    dat.Ephi = mhd.interp(esph.phi, xyzmesh).squeeze()
    dat.Erad = mhd.interp(esph.rad, xyzmesh).squeeze()
    return dat

client = Client()
view = client[:]
# Push the variables
view['xyzmesh'] = xyzmesh
retval = view.map_sync(interpolate_mhd, sim)
client.close()

# dat = map(interpolate_mhd, sim)

dat.Bz = array( [d.Bz for d in retval] )
dat.Ephi = array( [d.Ephi for d in retval] )
dat.Erad = array( [d.Erad for d in retval] )

writehdf5(oname, '/', dat, new=True)

plt.close('all')
ax = plt.subplot(221, polar=True)
im = ax.contourf(phirng, radrng, transpose((Ephi[0,:,:])), 100)
ax.set_rmin(0.)
colorbar(im)
