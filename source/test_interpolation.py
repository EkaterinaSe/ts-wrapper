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

def mkdir(s):
    if os.path.isdir(s):
        shutil.rmtree(s)
    os.mkdir(s)

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
    # size of the model grid
    M = len(all_parameters['file'])
    # how many points to remove?
    count = int(argv[1])

    "Generate random indexes to exclude models from the grid"
    ind = ( np.random.uniform(low=0.0, high=M-1.0, size=count) ).astype(int)
    excluded = np.full(M, False)
    excluded[ind] = True

    " Copy the grid of models omitting a random point"
    all_short = {}
    for k in all_parameters:
        all_short.update({ k : all_parameters[k][~excluded] })

    " Create array of input parameters for interpolation "
    interpol_parameters = { 'teff':None, 'logg':None, 'feh':None, 'vturb':None }
    for k in interpol_parameters:
        interpol_parameters[k] = all_parameters[k][excluded]

    interp_f, pars_to_interpolate, models_mask = NDinterpolate(interpol_parameters, all_short)

    "Create directory to save interpolated models"
    today = datetime.date.today().strftime("%b-%d-%Y")
    save_interpol_dir = f"./interpolated_models_{today}"
    mkdir(save_interpol_dir)

    comparison = {}
    struct_keys = ['tau500', 'temp', 'ne', 'vturb']
    for i in np.where(excluded)[0]:
        # mask = np.full(M, False)
        # mask[i] = True


        name_exc = all_parameters['file'][i].replace('.mod', '')
        # poitns are normalised to the max value of parameter
        # (stored in pars_to_interpolate)
        point = np.array([ all_parameters[k][i] / pars_to_interpolate[k] for k in pars_to_interpolate])

        # ma.depth_scale, ma.temp, ma.ne, ma.vturb
        value = all_parameters['structure'][i]

        # it returns array of interpolated models, so for one take 0th element
        interpolated_structure = interp_f(point)[0]

        # write to file
        ma_interpol = model_atmosphere()
        ma_interpol.id = name_exc
        ma_interpol.header = f"This model has been interpolated using routines from E.Magg. \n* Date: {today}"

        ma_interpol.depth_scale_type = 'TAU500'
        ma_interpol.logg = all_parameters['logg'][i]
        ma_interpol.teff = all_parameters['teff'][i]
        ma_interpol.feh = all_parameters['feh'][i]

        ma_interpol.depth_scale = interpolated_structure[0]
        ma_interpol.ndep = len(ma_interpol.depth_scale)

        ma_interpol.temp = interpolated_structure[1]
        ma_interpol.ne = interpolated_structure[2]
        ma_interpol.vmac = np.zeros(len(ma_interpol.depth_scale))
        ma_interpol.vturb = interpolated_structure[3]

        ma_interpol.write( f"{save_interpol_dir}/atmos.{name_exc}_interpol", format='m1d')


        # write original model to file
        ma_orig = model_atmosphere()
        ma_orig.id = name_exc
        ma_orig.header = f" "

        ma_orig.depth_scale_type = 'TAU500'
        ma_orig.logg = all_parameters['logg'][i]
        ma_orig.teff = all_parameters['teff'][i]
        ma_orig.feh = all_parameters['feh'][i]

        ma_orig.depth_scale = value[0]
        ma_orig.ndep = len(ma_interpol.depth_scale)

        ma_orig.temp = value[1]
        ma_orig.ne = value[2]
        ma_orig.vmac = np.zeros(len(ma_interpol.depth_scale))
        ma_orig.vturb = value[3]

        ma_orig.write( f"{save_interpol_dir}/atmos.{name_exc}_original", format='m1d')


        comparison.update( { name_exc : {} } )
        for j in range(len(struct_keys)):
            comparison[name_exc].update( { struct_keys[j] : {'interpol' : interpolated_structure[j], 'orig' : value[j] } } )
        for k in ['teff', 'logg', 'feh', 'vturb']:
            comparison[name_exc].update( { k :  all_parameters[k][i] } )

    with open('./interpolation_MC_test.pkl', 'wb') as f:
        pickle.dump(comparison, f)


    exit(0)
