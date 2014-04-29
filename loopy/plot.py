
import os
from numpy import array, ix_
from matplotlib.pyplot import imshow, colorbar, gca
from matplotlib.pyplot import clim as mclim
from matplotlib.pyplot import savefig as _savefig
from matplotlib.patches import Wedge

__all__ = ['savefigs', 'panel_label', 'implanet']

def savefigs(fname, outdir='fig', formats=['jpg', 'png']):
    """
    Save thumbnail figure in the current directory of hires in `outdir`.

    :param fname: File name without extension.
    :param outdir: Output directory of the hi-res figures.
    :param formats: Requested formats (jpg|png|eps).
    """
    formats = formats if isinstance(formats, list) else [formats]
    for fdir in [outdir + '/' + fmt for fmt in formats]:
        if os.path.exists(fdir):
            continue

        print 'Creating figure directory: ' + fdir
        os.makedirs(fdir)

    # Thumbnail goes into outdir only
    _savefig(outdir + '/' + fname + '.jpg', dpi=100)

    # Example: fig/jpg/eqplane.jpg
    oname = lambda fmt: outdir + '/' + fmt + '/' + fname + '.' + fmt

    if 'eps' in formats:
        _savefig(oname('eps'), bbox_inches='tight', pad_inches=0)
    if 'jpg' in formats or 'jpeg' in formats:
        _savefig(oname('jpg'), dpi=300, bbox_inches='tight', pad_inches=0)
    if 'png' in formats:
        _savefig(oname('png'), dpi=300, bbox_inches='tight', pad_inches=0, \
                 transparent=True)

def panel_label(ax, label, loc='ul', fsize=14):
    """
    Draw a small box in the corner with a panel label (i.e. `a)`).

    :param ax:    Axes handle.
    :param label: Text of the label (i.e. `a)`).
    :param loc:   Location of the label: `ul` (default), `ur`, `ll`, `lr`.
    :param fsize: Font size (default 14).
    """
    locations={'ul' : (0.01,0.95), 'ur' : (0.99,0.95), \
               'll' : (0.01,0.05), 'lr' : (0.99,0.05)}
    assert loc in locations.keys()

    valign = 'top' if loc[0] == 'u' else 'bottom'
    halign = 'left' if loc[1] == 'l' else 'right'

    ax.text(locations[loc][0], locations[loc][1], label,
            verticalalignment=valign, horizontalalignment=halign,
            transform=ax.transAxes, fontsize=fsize,
            bbox=dict(facecolor='white', alpha=0.6))

def implanet(data, clim=None, dx=1., rx=0.0, radius=-1., scale='none', \
             cmap=None):
    nx = data.shape[::-1] # 2D
    lx = array(nx) * array(dx)    # 2D
    # # lx -= lx * array(rx) # Shift to planet centre
    extent = array([0.,lx[0], 0., lx[1]])
    extent[ix_((0,2))] -= lx * array(rx)
    extent[ix_((1,3))] -= lx * array(rx)
    extent *= -1.
    if scale == 'planet' and radius > 0:
        extent /= radius
        radius = 1.

    # Alternatives are `pcolor` and faster `pcolormesh`
    imshow(data, origin='lower', interpolation='none', aspect='equal', \
           extent=extent, cmap=cmap)
    if radius > 0.:
        dual_half_circle(radius, angle=-90)

    colorbar()
    if clim != None:
        mclim(clim)

def dual_half_circle(radius, angle=0, ax=None, colors=('w','k'), **kwargs):
    """
    Add two half circles to the axes *ax* (or the current axes) at the lower
    left corner of the axes with the specified facecolors *colors* rotated at
    *angle* (in degrees).
    """
    if ax is None:
        ax = gca()
    # kwargs.update(transform=ax.transAxes, clip_on=False)
    center = (0, 0)
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]
