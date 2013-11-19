
import h5py
import tables as pt
import numpy as np

def loadhdf5(fname,variable):

    h5file=h5py.File(fname,'r')
    data = h5file[variable][...]
    h5file.close()
    return data

def writehdf5(fname, loc, data, new=False):
    """
    .. py:function:: writehdf5(fname, name, data, new=False)

    The function writes HDF5 file using PyTables with compression.

    :param fname: name of the HDF5 file
    :param loc: location inside of HDF5 file (e.g. /new/Bz)
    :param data: the actual data to be stored
    :type data: dict or ndarray otherwise will be converted into ndarray
    :rtype: none
    """
    mode = 'w' if new else 'a'
    filters = pt.Filters(complevel=6)
    h5f = pt.File(fname, mode)

    if isinstance(data, dict):
        h5f.close()
        for key in data.keys():
            writehdf5(fname, name+'/'+key, data[key], new=False)
        return

    path = name.rsplit('/',1)
    root = path[0] if np.size(path) > 1 else '/'
    root = root if root.startswith('/') else '/' + root
    name = name if np.size(path) == 1 else path[1]

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
