# external
import os
from sys import argv
import shutil
import subprocess
import datetime
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import pickle
import glob
import time
import warnings
# local
import convolve
import read_config
from atmos_package import model_atmosphere
from model_atm_interpolation import *



"""
Test how LinearNDInterpolator is doing at interpolating MARCS grid
by excluding N random points from the grid and interpolating to those points
then compare
"""

if __name__ == '__main__':

    atmos_path = '/Users/semenova/phd/projects/ts-wrapper/input/atmos/MARCS/all/'

    "Read all model atmospheres"
    all_parameters = get_all_ma_parameters(atmos_path, \
                                            format = 'marcs', debug=True)

    all_new = {}
    for k in all_parameters:
        all_new.update({ k : all_parameters[k] })

    all_parameters = all_new
    # size of the model grid
    M = len(all_parameters['file'])
    # how many points to remove?
    count = int(argv[1])

    "Generate random indexes to exclude models from the grid"
    ind = ( np.random.uniform(low=0.0, high=M-1.0, size=count) ).astype(int)
    excluded = np.full(M, False)
    excluded[ind] = True

    " Copy the grid of models omitting those random points"
    all_short = {}
    for k in all_parameters:
        all_short.update({ k : all_parameters[k][~excluded] })

    " Create array of input parameters for interpolation "
    interpol_parameters = { 'teff':[], 'logg':[], 'feh':[], 'vturb':[] }
    for i in np.where(excluded)[0]:
        for k in interpol_parameters:
            interpol_parameters[k].append(all_parameters[k][i])

    interp_f, pars_to_interpolate, models_mask = NDinterpolate(interpol_parameters, all_short)

    comparison = {}
    strtuct_keys = ['tau500', 'temp', 'ne', 'vturb']
    for i in range(len(np.where(excluded)[0])):
        name_exc = all_parameters['file'][i]
        # poitns are normalised to the max value of parameter
        # (stored in pars_to_interpolate)
        point = np.array([ all_parameters[k][i] / pars_to_interpolate[k] for k in pars_to_interpolate])

        # ma.depth_scale, ma.temp, ma.ne, ma.vturb
        value = all_parameters['structure'][i]

        # it returns array of interpolated models, so for one take 0th element
        interpolated_structure = interp_f(point)[0]

        # difference in structure of original and interpolated
        diff =  (value - interpolated_structure) / value

        comparison.update( { name_exc : {} } )
        for j in range(len(strtuct_keys)):
            comparison[name_exc].update( { strtuct_keys[j] : diff[j] } )

    with open('./interpolation_MC_test.pkl', 'wb') as f:
        pickle.dump(comparison, f)


    exit(0)
