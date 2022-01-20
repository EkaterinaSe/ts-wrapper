# external
import os
from sys import argv
import shutil
import subprocess
import datetime
import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d
from scipy.spatial import Delaunay
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
    save_file = f"{models_path}/all_models_save.pkl"
    params = {
    'teff':[], 'logg':[], 'feh':[], 'vturb':[], 'file':[], 'structure':[], 'mass':[]
    }

    if os.path.isfile(save_file) and os.path.getsize(save_file) > 0:
        with open(save_file, 'rb') as f:
            params = pickle.load(f)
    else:
        print(f"Checking all model atmospheres under {models_path}")
        d_sc_new = np.linspace(-5, 2, 100)

        with os.scandir(models_path) as all_files:
            for entry in all_files:
                if not entry.name.startswith('.') and entry.is_file():
                # try:
                    file_path = models_path + entry.name
                    ma = model_atmosphere()

                    ma.read(file_path, format=format)

                    if ma.mass <= 1.0:

                        params['teff'].append(ma.teff)
                        params['logg'].append(ma.logg)
                        params['feh'].append(ma.feh)
                        params['vturb'].append(ma.vturb[0])
                        params['mass'].append(ma.mass)

                        params['file'].append(entry.name)

                        ma.temp = np.log10(ma.temp)
                        ma.ne = np.log10(ma.ne)

                        # bring all values to the same depth_scale (tau500)
                        for par in ['temp', 'ne', 'vturb']:
                            f_int = interp1d(ma.depth_scale, ma.__dict__[par], fill_value='extrapolate')
                            ma.__dict__[par] = f_int(d_sc_new)
                        ma.depth_scale = d_sc_new

                        params['structure'].append(np.vstack((ma.depth_scale, ma.temp, ma.ne, ma.vturb)))


                    # except: # if it's not a model atmosphere file, or format is wrong
                    #         if debug:
                    #             print(f"Cound not read model file {entry.name} for model atmosphere")

        for k in params:
            params[k] = np.array(params[k])

        " Check if any model atmosphere was successfully read "
        if len(params['file']) == 0:
            raise Exception(f"no model atmosphere parameters were retrived from files under {models_path}.\
Try setting debug = 1 in config file. Check that expected format of model atmosphere is set correctly.")
        else:
            print(f"{len(params['file'])} model atmospheres in the grid")

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


def NDinterpolate(inp_par, all_par):
    """
    # NOTE: can not extrapolate outside of the grid
    """

    N = int(len(inp_par.keys()) )# number of input parameters
    M = int(len( list( inp_par.values())[0] )) # number of models to create

    " Exclude degenerate parameters (aka the same for all grid points) "
    points = []
    # dict of parameters used for interpolation
    # and their values for normalising parameter space
    params_to_interpolate = {}
    for k in inp_par:
        if not max(all_par[k]) == min(all_par[k]):
            points.append(all_par[k] / max(all_par[k]) )
            params_to_interpolate.update( { k :  max(all_par[k])} )
        else:
            print(f"The grid is degenerate in parameter {k}")
    points = np.array(points).T

    # print("Creating a thing...")
    # tri = Delaunay(points)

    " Check for repeatative points in the grid"
    test = []
    for key in all_par:
        if key not in ['structure', 'file']:
            test.append( [all_par[k]] )
    if len(np.unique(test, axis=1)) != len(test):
        raise Warning(f"Grid has repeatative points.")

    " How many models out of the requested list can be created? "
    doable = np.full(M, False)
    for i in range(M):
        "Skip if outside of the grid"
        outside = np.array( [np.logical_or(inp_par[k][i] > max(all_par[k]),inp_par[k][i] < min(all_par[k]) ) \
         for k in params_to_interpolate] )
        if outside.any():
            print(f"{[ [k,inp_par[k][i]]  for k in params_to_interpolate]} \
outside of the grid, skipping interpolation")
        else:
            doable[i] = True
    # print(f" {len(np.where(doable)[0])} models out of requested {M} \
# can be created")

    "Create interpolator function that interpolates model atmospheres structure"
    values = all_par['structure']
    interp_f = LinearNDInterpolator(points, values)

    return interp_f, params_to_interpolate, doable


def NDinterpolate_NLTE_grid(interpol_parameters, nlte_data):
    N = int(len(interpol_parameters.keys()) ) # number of input parameters
    M = int(len( list( interpol_parameters.values())[0] )) # number of models to create


    points = []
    # dict of parameters used for interpolation
    # and their values for normalising parameter space
    params_to_interpolate = {}
    for k in interpol_parameters:
        if not max(nlte_data[k]) == min(nlte_data[k]):
            points.append(nlte_data[k] / max(nlte_data[k]) )
            params_to_interpolate.update( { k :  max(nlte_data[k])} )
        else:
            print(f"The grid is degenerate in parameter {k}")
    points = np.array(points).T

    values = nlte_data['depart']
    interp_f = LinearNDInterpolator(points, values)

    return interp_f, params_to_interpolate


def interpolate_ma_grid(atmos_path, atmos_format, debug):
    all_parameters = get_all_ma_parameters(atmos_path,  \
                        format=atmos_format, debug=debug)

    input_parameters = {
        'teff' : [8000, 8000],
        'logg' : [3.5, 3.5],
        'feh'  : [-10.2, -2.0],
        'vturb': [3.0, 2.0]    }

    interp_f, pars_to_interpolate, models_mask = NDinterpolate(input_parameters, all_parameters)




if __name__ == '__main__':
    exit(0)
