
import os
from matplotlib.pyplot import savefig as _savefig

__all__ = ['savefigs']

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
    
