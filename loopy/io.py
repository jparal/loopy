
import tables as pt
import numpy as np
import loopy as lpy
import shutil as sh # move
import os.path as pth # exists

def readhdf5(fname, path='/'):
    """
    .. py:function:: writehdf5(fname, path='/')

    The function traverse HDF5 files and creates structured dictionary.

    :param fname: File name to read.
    :param path: Root path from where to start reading.
    :rtype: loopy.struct (i.e. dictionary) or variable
    """
    def _traverse_tree(h5f, path):
        # Remove double slashes and the last one
        path = '/'+'/'.join(filter(None, path.split('/')))
        gloc = ''.join(path.rpartition('/')[0:2])
        name = path.rpartition('/')[2]

        # We want to read a single variable
        groups = h5f.listNodes(where=gloc, classname='Group')
        nodes = h5f.listNodes(where=gloc)
        leafs = [n for n in nodes if n not in groups]
        leaf = [n for n in leafs if n.name == name]
        if len(leaf) == 1:
            return leaf[0].read()

        dat = lpy.struct()
        for node in h5f.listNodes(where=path):
            name = node._v_name
            dat[name] = _traverse_tree(h5f, path+'/'+name)

        return dat

    h5f = pt.File(fname, 'r')
    dat = _traverse_tree(h5f, path)
    h5f.close()
    return dat

def writehdf5(fname, data, path='/', append=False, backup=False):
    """
    .. py:function:: writehdf5(fname, data, path='/', append=False)

    The function writes HDF5 file using PyTables and CArray.
    This is high level function which shoud handle the most common scenarios.

    :param fname: name of the HDF5 file
    :param path: location inside of HDF5 file (e.g. /new/Bz)
    :param data: the actual data to be stored
    :type data: dict or ndarray otherwise will be converted into ndarray
    :param append: Should the data be appended to an existing file?
    :param backup: This argument if True rename the file to .bak instead of
    overwriting the file.

    :rtype: none
    """
    if backup and pth.exists(fname):
        sh.move(fname, fname+'.bak')

    mode = 'a' if append else 'w'
    filters = pt.Filters(complevel=6)
    h5f = pt.File(fname, mode)

    # Remove double slashes and the last one
    path = '/'+'/'.join(filter(None, path.split('/')))
    dloc = path.rsplit('/',1)
    root = dloc[0] if np.size(dloc) > 1 else '/'
    root = root if root.startswith('/') else '/' + root
    name = path if np.size(dloc) == 1 else dloc[1]

    if isinstance(data, dict):
        h5f.close()
        for key in data.keys():
            writehdf5(fname, data[key], path=path+'/'+key, append=True)
        return

    if not isinstance(data, np.ndarray):
        data = np.array(data, ndmin=1)

    atm = pt.Atom.from_dtype(data.dtype)
    arr = h5f.createCArray(root, name, atm, data.shape, \
                     createparents=True, filters=filters)
    arr[:] = data
    h5f.close()
    return


from pyhdf.SD import SD, SDC

def loadhdf4(fname,variable):

    data_set = SD(fname, SDC.READ)
    return data_set.select(variable)[:]
