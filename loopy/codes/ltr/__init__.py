"""
==
LTR: LFM-TIEGCM-RCM (LTR)
==

.. currentmodule:: loopy.codes.ltr

The LFM is a global MHD model written by John Lyon, Joel Fedder and Clark Mobarry. It is primarily used to better understand Earth's Magnetosphere. The LFM can be coupled with a variety of models to aide in understanding of the physical processes in Geospace , such as:

Magnetosphere Ionosphere Coupler/Solver (MIX)
Thermosphere Ionosphere Electrodynamic General Circulation Model (TIEGCM)
Rice Convection Model (RCM)
The coupled combination of models is called the LTR (LFM-TIEGCM-RCM).
The documentation on this wiki contains limited scientific descriptions of the models.  It is intended as a guide to compile & execute the code(s) and documents some of the post-processing tools which have been developed for LTR.
"""

from __future__ import division, print_function, absolute_import

from .mhd import *
from .file import *

__all__ = [s for s in dir() if not s.startswith('_')]
