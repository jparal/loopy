
import os
from matplotlib.pyplot import savefig as _savefig

__all__ = ['savefigs', 'panel_label']

def savefigs(fname, outdir='fig', formats=['jpg', 'png']):
    """
    Save thumbnail figure in the current directory of hires in `outdir`.

    :param fname: File name without extension.
    :param outdir: Output directory of the hi-res figures.
    :param formats: Requested formats (jpg|png|eps).
    """
    if not os.path.exists(outdir):
        print 'Creating figure directory: ' + outdir
        os.makedirs(outdir)

    # Thumbnail
    _savefig(fname + '.jpg', dpi=100)

    oname = 'fig/' + fname

    if 'eps' in formats:
        _savefig(oname + '.eps', bbox_inches='tight', pad_inches=0)
    if 'jpg' in formats or 'jpeg' in formats:
        _savefig(oname + '.jpg', dpi=300, bbox_inches='tight', pad_inches=0)
    if 'png' in formats:
        _savefig(oname + '.png', dpi=300, bbox_inches='tight', pad_inches=0, \
                 transparent=True)

def panel_label(ax, label):
    #ydiff = ax.transAxes.transform((0,1))-ax.transAxes.transform((0,0))
    ax.text(0.01, 0.96, label, horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white'))
