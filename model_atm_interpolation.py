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


def get_all_ma_parameters(models_path, format='m1d', debug = False):
    """
    Get a list of all available model atmopsheres and their parameters
    for later interpolation
    also upload the whole grid, or create the file for the quick upload
    If no list is available, create one by scanning through all available models
    """
    # list_file = f"{models_path}/all_models_list.txt"
    save_file = f"{models_path}/all_models_save.pkl"

    params = {
    'teff':[], 'logg':[], 'feh':[], 'vturb':[], 'file':[], 'structure':[]
    }

    if os.path.isfile(save_file) and os.path.getsize(save_file) > 0:
        with open(save_file, 'rb') as f:
            params = pickle.load(f)
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
                        params['vturb'].append(ma.vturb[0])
                        params['file'].append(entry.name)
                        params['structure'].append(np.vstack((ma.depth_scale, ma.temp, ma.ne, ma.vturb)))
                    except: # if it's not a model atmosphere file, or format is wrong
                        if debug:
                            print(f"Cound not read model file {entry.name} for model atmosphere")

        for k in params:
            params[k] = np.array(params[k])

        " Check if any model atmosphere was successfully read "
        if len(params['file']) == 0:
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
        "Dump all in one file (only done once)"
        with open(save_file, 'wb') as f:
            pickle.dump(params, f)
    return params


# def create_cube(input_par, all_par, debug=False):
    """
    Find a cube in the grid of model atmospheres for interpolation
    by mimnimising the [normalised] distance from input parameters
    to the grid points

    Input:
    input_par (dict) -- parameters to which a grid should be interpolated, 1 point
    all_par (dict) -- parameters of the full model atmospheres grid,
                        from call to get_all_ma_parameters()
    Return:
    (list of) index of selected cube points in all_par dictionary,
    (list of) file names of selected for interpolation models
    (list of) parameters to interpolate over
    """

    N = len(input_par.keys()) # number of input parameters
    M = [len(v) for v in all_par.values()][0] # number of points in the grid

    "Compute individual (parameter specific) distance to each grid point"
    dist = np.zeros( shape=(N,M) )
    i = 0
    for k, value in input_par.items():
        dist[i, : ] = ( all_par[k] - value )/ max(abs(all_par[k])) # normalise the distance
        i += 1

    """
    Compute total distance, collapsing the dimension of parameters.
    This is our main diagnostic for finding the interpolation 'cube'
    """
    tot_dist = np.zeros(M)
    for i in range(N):
        tot_dist[:] = tot_dist[:] + np.sqrt(dist[i, :]**2)

    " Check if any parameters are exactly as the ones at [some] grid point (or within 0.5%) "
    threshold =  0.005 * N
    if np.any(tot_dist < threshold):
        pos = np.argwhere(tot_dist < threshold)
        if len(pos) > 1: # if more than one grid point matches,
        # maybe they differ in parameter we don't iterate over, e.g. alpha/fe
            message = f"Found more than one point in the grid matching input parameters within 0.5%\n"
            for k in input_par:
                message = message + f"{k}={input_par[k]}\t"
            message = message + "\n"
            for p in pos:
                message = message + f" {all_par['file'][p]} \n"
            raise Exception(message)
        else: # congrats, we don't need to iterate, this is our model
            pos = pos[0][0]
            if debug:
                message = f"{all_par['file'][pos]} exactly matches "
                for k, v in input_par.items():
                    message = message + f"{k}={v}\t"
                print(message)
            model_int_path = all_par['file'][pos]
            return pos, model_int_path, None
    # if interpolation is needed, build a cube
    else:
        degenerate_params = []

        cube_size = 2**N
        ind = np.argpartition(tot_dist, cube_size)[:cube_size]

        i = 0
        for k, v in input_par.items():
            if v >= max(all_par[k]) or v <= min(all_par[k]):
            # if max(all_par[k][ind ]) == min(all_par[k][ind ]):
                print(f"{k}  is degenerate")
                degenerate_params.append(k)
                cube_size = int(cube_size/2)
            i += 1
            # refine
            ind = np.argpartition(tot_dist, cube_size)[:cube_size]
            # print(all_par['file'][ind])

        i = 0
        a = np.full(M, False)
        tot_dist = np.zeros(M)

        for k in input_par:
            if not k in degenerate_params:
                closest_value =  np.partition(abs(dist[i, :]), 2)[:2]
                ind = np.in1d(abs(dist[i, :]), closest_value)
                b = np.full(M, False)
                b[ind] = True
                a = ~(~a*~b)#?????
                tot_dist[a] = tot_dist[a] +np.sqrt(dist[i, a]**2)

                # print(abs(dist[i, :]), abs(dist[j, :]), abs(dist[i, :])* abs(dist[j, :]))
                # ind.extend(np.argpartition(abs(dist[i, :])*abs(dist[j, :]), 2)[:2])
            i += 1

def interpolate_cube(input, all):
    """
    Interpolate N model atmospheres with respect to requested parameteres

    Input:
    (dict) input -- parameters to interpolate over and their values
    (dict) all -- all info about the grid, including paramteres of all points
    """


def NDinterpolate(inp_par, all_par):
    """
    """

    " Exclude degenerate parameters (the same for all grid points) "
    points = []
    for k in inp_par:
        points.append(all_par[k])
    points = np.array(points).T
    # values = points
    values = all_par['structure']
    print(np.shape(points), np.shape(values))
    interp_f = LinearNDInterpolator(points, values)
    # print(interp_f(3.1567890, -2))



def interpolate_ma_grid(atmos_path, atmos_format, debug):
    all_parameters = get_all_ma_parameters(atmos_path,  \
                        format=atmos_format, debug=debug)

    input_parameters = {
        'teff' : 8000,
        'logg' : 3.5,
        'feh'  : -2.2,
        'vturb': 3.0
    }
    NDinterpolate(input_parameters, all_parameters)

    #
    # # ind, names, params_to_interpolate = \
    # create_cube(input_parameters, all_parameters, debug=debug)
    # if params_to_interpolate is not None:
    #     print(f"Ready to interpolate?")
    #     " Clean up the input a bit"
    #     inp, all = {}, {}
    #     for k in params_to_interpolate:
    #         inp.update({ k : input_parameters[k] })
    #     for k in all_parameters:
    #         all.update({ k : [ all_parameters[k][i] for i in ind ] })
    #     interpolate_cube(inp, all)
    #
    # else: # model exists, no interpolation
    #     ma = model_atmosphere(atmos_path + names, format=atmos_format)
    #     # TODO: write to file here?
    #     return ma
