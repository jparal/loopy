
import datetime as dt
import numpy as np
import matplotlib.dates as md

__all__ = ['datespace']

def datespace(dbeg, dend, ndate):
    dbeg = dt.datetime(*dbeg)
    dend = dt.datetime(*dend)
    dspace = np.linspace(md.date2num(dbeg), md.date2num(dend), ndate)
    dspace = md.num2date(dspace)
    return dspace
