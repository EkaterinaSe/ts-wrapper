# external
import os
from sys import argv
import shutil
import subprocess
import datetime
import numpy as np
import pickle
import glob
import time
import warnings
# local
import convolve
import read_config
from atmos_package import model_atmosphere


def get_all_ma_parameters(models_path, format='m1d', debug = False):
    """
    Get a list of all available model atmopsheres and their parameters
    for later interpolation
    If no list is available, create one by scanning through all available models
    """
    list_file = f"{models_path}/all_models.txt"
    params = {
    'teff':[], 'logg':[], 'feh':[], 'path':[], 'file':[]
    }
    if os.path.isfile(list_file) and os.path.getsize(list_file) > 0:
        params['teff'], params['logg'], params['feh'] = np.loadtxt(list_file, usecols=(0,1,2), unpack=True)
        params['file'] = np.loadtxt(list_file, usecols=(3), dtype=str)

    else:
        print(f"Checking all model atmospheres under {models_path}")

        with os.scandir(models_path) as all_files:
            for entry in all_files:
                if not entry.name.startswith('.') and entry.is_file():
                    try:
                        file_path = models_path + entry.name
                        ma = model_atmosphere(file_path, format=format)
                        params['teff'].append(ma.teff)
                        params['logg'].append(ma.logg)
                        params['feh'].append(ma.feh)
                        params['path'].append(file_path)
                        params['file'].append(entry.name)

                    except: # if it's not a model atmosphere file, or format is wrong
                        if debug:
                            print(f"Cound not read model file {entry.name} for model atmosphere")

        for k in params:
            params[k] = np.array(params[k])

        " Check if any model atmosphere was successfully read "
        if len(params['teff']) == 0:
            raise Exception(f"no model atmosphere parameters were retrived from files under {models_path}.\
Try setting debug = 1 in config file. Check that expected format of model atmosphere is set correctly.")

        "Print UserWarnings about any NaN in parameters"
        for k in params:
            try: # check for NaNs in numeric values:
                if np.isnan(params[k]).any():
                    pos = np.where(np.isnan(params[k]))
                    for p in pos:
                        message = f"NaN in parameter {k} from model atmosphere {params['path'][p]}"
                        warnings.warn(message, UserWarning)
            except TypeError: # ignore other [non-numerical] keys, such as path, name, etc
                pass
        # TODO: should I exclude NaNs here or will they be talen care of when interpolating, or searching for interp. cube?
        #

        with open(list_file, 'w') as f:
            count = [len(v) for v in params.values()][0]
            for i in range(count):
                f.write(f"{params['teff'][i]:8.0f} {params['logg'][i]:8.3f} \
        {params['feh'][i]:8.3f} {params['file'][i]} \n")

    return params


def create_cube(input_par, all_par, debug=False):
    """
    Find a cube in the grid of model atmospheres for interpolation
    by mimnimising the [normalised] distance from input parameters
    to the grid points

    Input:
    input_par (dict) -- parameters to which a grid should be interpolated, 1 point
    all_par (dict) -- parameters of the full model atmospheres grid,
                        from call to get_all_ma_parameters()
    """

    N = len(input_par.keys())
    M = [len(v) for v in all_par.values()][0]

    dist = np.zeros( shape=(N,M) )
    i = 0
    for k, value in input_par.items():
        dist[i, : ] = ( all_par[k] - value )/ max(abs(all_par[k])) # normalise the distance
        i += 1

    """
    Check if any parameter has the same distance to all grid points
    This could mean that
    parameter is [exactly] equally far from any grid point,
            e.g. the whole grid has the same temperature at each point
    Interpolating over this parameter is obsolette
    NOTE: this is more of an exception, than a rule

    Save the parameters over which interpolation needs to be done,
    and number of points for the 'cube'
    """
    n_dim = 0
    params_to_interpolate = []
    i = 0
    for k, value in input_par.items():
        if np.max(dist[i, :]) == np.min(dist[i, :]):
            if debug:
                print(f"{k} is {dist[i, :][0] * max(abs(all_par[k])):.2f} far from EVERY point in the grid")
            else:
                pass
        else: # if not, add a dimension, +2 points to the interpolation n-D cube
            n_dim += 1
            params_to_interpolate.append(k)
        i += 1


    " Collapse the dimension of parameters, only consider iteratable parameters "
    tot_dist = np.zeros(M)
    i = 0
    for k, value in input_par.items():
        if k in params_to_interpolate:
            tot_dist[:] = tot_dist[:] + np.sqrt(dist[i, :]**2)
        i += 1

    " Check if any parameters are exactly as the ones at [some] grid point (or within 0.5%) "
    if np.any(tot_dist < 0.005):
        pos = np.argwhere(tot_dist < 0.005)

        if len(pos) > 1: # if more than one grid point matches,
        # maybe they differ in parameter we don't iterate over, e.g. alpha/fe
            message = f"Found more than one point in the grid matching input parameters within 0.5%\n"
            for k in params_to_interpolate:
                message = message + f"{k}={input_par[k]}\t"
            message = message + "\n"
            for p in pos:
                message = message + f"{all_par['file'][p]} \n"
            raise Exception(message)
        else: # congrats, we don't need to iterate
            pos = pos[0][0]
            if debug:
                message = f"{all_par['file'][pos]} exactly matches "
                for k, v in input_par.items():
                    message = message + f"{k}={v}\t"
                print(message)
            model_int_path = all_par['file'][pos]
    # if interpolation is needed, build a cube
    else:
        cube_size = 2**n_dim
        # find N=cube_size points with the smallest distance
        ind = np.argpartition(tot_dist, cube_size)[:cube_size]
        if debug:
            print(f"Interpolation 'cube': ")
            print(f"Interpolating over {len(params_to_interpolate):.0f} parameter(s): {params_to_interpolate}")
            for i in ind:
                message = f""
                for k in input_par:
                    message = message + f"{k}={all_par[k][i]} "
                print(message)


    return


def interpolate_ma_grid(setup):
    all_parameters = get_all_ma_parameters(setup.atmos_path,  \
                        format=setup.atmos_format, debug=setup.debug)

    input_parameters = {
        'teff' : 7500,
        'logg' : 4.0,
        'feh'  : -2.4
    }
    create_cube(input_parameters, all_parameters, debug=setup.debug)
