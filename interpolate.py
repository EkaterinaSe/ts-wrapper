#!/opt/anaconda3/bin/python
# external
import os
from sys import argv
import shutil
import subprocess
import datetime
import numpy as np
import glob
import time
# # local
# import read_config
# from atmos_package import read_atmos_marcs, model_atmosphere
# from read_nlte import grid_to_ts


def find_a_cube(par, par_all):
    """ Find 8 models (a cube) around given points in Teff-logg-FeH space """
    par_all = np.array(par_all).T
    par = np.array(par)
    # t = par[0]
    # lg = par[1]
    # feh = par[2]
    #
    # t-all = par-all[0]
    # lg-all = par-all[1]
    # feh-all = par-all[2]

    # normalised difference
    pos = (par_all - par) / par

    pos = np.sum(pos/par, axis=1)
    pos_min = np.sort(np.abs(pos))[:8]

    print(pos_min)

    # for i in range(len(par_all)):
    #     print(par_all[i])
    #     # pos.append( par_all[:,i] - par )
    # print(pos)


    print(pos)
    return


if __name__ == '__main__':
    all_models = np.loadtxt(argv[1], usecols=(0), dtype=str)
    teff, logg, feh = np.loadtxt(argv[1], usecols=(1,2,3), unpack=True)

    par_to_interp = [5550, 4.35, -0.1 ]

    find_a_cube(par_to_interp, [teff, logg, feh])

exit(0)
