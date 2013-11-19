
from __future__ import division, print_function, absolute_import

from os.path import exists, realpath, expanduser
import numpy as np
import datetime as dtime
from .mhd import MHDFile

__all__ = ['FileIter']

class FileIter(object):
    def __init__(self, prefix, start, stop, step_sec, file_obj=MHDFile):
        self.prefix = prefix
        start = np.array(start).copy()
        start.resize(6)
        stop = np.array(stop).copy()
        stop.resize(6)

        self.start = dtime.datetime(*start)
        self.stop = dtime.datetime(*stop)
        self.step_sec = step_sec
        self.file_obj = file_obj
        self._iter = 0

    def __iter__(self):
        self._iter = 0
        return self

    def next(self):
        dt=dtime.timedelta(seconds=self.step_sec) * self._iter
        date = self.start + dt
        if date > self.stop:
            raise StopIteration

        self._iter += 1
        # Example: baker_mhd_2012-10-07T00-09-00Z.hdf
        date_str='{:%Y-%m-%dT%H-%M-%SZ}'.format(date)
        fname = self.file_obj.get_fname_format()
        fname = fname.format(prefix=self.prefix, date=date_str)
        fname = realpath(expanduser(fname))

        if not exists(fname):
            raise IOError('File not found: %s' % fname)
        else:
            return self.file_obj(fname)
