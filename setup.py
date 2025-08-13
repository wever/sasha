from distutils.core import setup
from distutils.filelist import findall
import os
import py2exe
from numpy import array,zeros,empty,exp,sqrt,pi,set_printoptions,arange,copy,dot,transpose,diag,fill_diagonal,shape,linspace,meshgrid,double,intc,asarray
from numpy.linalg import eigh
from math import factorial
from time import clock
import numpy.ctypeslib as npc
import ctypes as ct
from itertools import chain
import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import QThread, SIGNAL

from matplotlib import use,get_py2exe_datafiles
use('Qt4Agg')
from matplotlib.pyplot import close,rcParams,figure,legend,xlabel,ylabel

from logging import Handler,getLogger,Formatter,DEBUG
import sys
import time
from HF import hartree_fock
import design
from get_basis import inp
import zmq.libzmq
setup(console=['sasha.py'],
    options={
        'py2exe': {
            'includes': ['zmq.backend.cython'],
            'excludes': ['zmq.libzmq'],
            'dll_excludes': ['libzmq.pyd','libiomp5md.pyd'],
            'packages': ['matplotlib','pytz'],
        }
    },
    data_files=get_py2exe_datafiles(),
)

