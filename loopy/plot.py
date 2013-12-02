
import os
from matplotlib.pyplot import savefig as _savefig

__all__ = ['savefigs', 'panel_label']

def savefigs(fname, dname='fig'):

    if not os.path.exists(dname):
        print 'Creating figure directory: ' + dname
        os.makedirs(dname)

    _savefig(fname + '.jpg', dpi=100)
    oname = 'fig/' + fname

    _savefig(oname + '.eps', bbox_inches='tight', pad_inches=0)
    _savefig(oname + '.jpg', dpi=300, bbox_inches='tight', pad_inches=0)
    _savefig(oname + '.png', dpi=300, bbox_inches='tight', pad_inches=0, \
            transparent=True)

def panel_label(ax, label):
    #ydiff = ax.transAxes.transform((0,1))-ax.transAxes.transform((0,0))
    ax.text(0.01, 0.96, label, horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white'))
