import numpy as np
import os
import shutil
from sys import argv, exit
import datetime
import glob
# local
# from model_atm_interpolation import get_all_ma_parameters, NDinterpolateGrid,preInterpolationTests
# from read_nlte import read_fullNLTE_grid, find_distance_to_point
# from atmos_package import model_atmosphere
# from run_ts import write_departures_forTS
# import cProfile
# import pstats
#
def atomicZ(el):
    if os.path.isfile('./atomic_numbers.dat'):
        el_z = np.loadtxt('./atomic_numbers.dat', usecols=(0))
        el_id = np.loadtxt('./atomic_numbers.dat', usecols=(1), dtype=str)
    else:
        print("Can not find './atomic_numbers.dat' file. Stopped.")
        exit(1)
    for i in range(len(el_id)):
        if el.lower() == el_id[i].lower():
            return el_z[i]


class ChemElement(object):
    def __init__(self, ID = ''):
        self.z = atomicZ(ID)
        self.nlte = False

        if ID.strip().lower() == 'fe' and self.z == 26:
            self.isFe = True
        else: self.isFe = False
