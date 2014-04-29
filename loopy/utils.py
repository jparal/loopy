
import datetime as dt
import numpy as np
import matplotlib.dates as md

__all__ = ['datespace']

def datespace(dbeg, dend, ndate):
    if not isinstance(dbeg, dt.datetime):
        dbeg = dt.datetime(*dbeg)
    if not isinstance(dend, dt.datetime):
        dend = dt.datetime(*dend)
    dspace = np.linspace(md.date2num(dbeg), md.date2num(dend), ndate)
    dspace = md.num2date(dspace)
    dspace = [dd.replace(tzinfo = None) for dd in dspace]
    return np.array(dspace)
