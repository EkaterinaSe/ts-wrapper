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

    if os.path.isfile(save_file) and os.path.getsize(save_file) > 0:
        with open(save_file, 'rb') as f:
            params = pickle.load(f)
    else:
        print(f"Checking all model atmospheres under {models_path}")
        d_sc_new = np.linspace(-5, 2, 100)

        params = {
        'teff':[], 'logg':[], 'feh':[], 'vturb':[], 'file':[], 'structure':[], 'structure_keys':[], 'mass':[]\
        }

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

                        params['structure'].append( np.vstack( (ma.depth_scale, ma.temp, ma.ne, ma.vturb )  ) )
                        params['structure_keys'].append( ['tau500', 'temp', 'ne', 'vturb'])

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

def preInterpolationTests(data, interpol_coords, dataLabel = 'default'):
    """
    Run multiple tests to catch possible exceptions
    that could affect the performance of the underlying
    Qnull math engine during Delaunay triangulation
    """

    " Check for degenerate parameters (aka the same for all grid points) "
    for k in interpol_coords:
        if max(data[k]) == min(data[k]):
            print(f"Grid {dataLabel} is degenerate in parameter {k}")
            exit()

    " Check for repetitive points within the requested coordinates "
    test = [ data[k] for k in interpol_coords]
    if len(np.unique(test, axis=1)) != len(test):
        print(f"Grid {dataLabel} with coordinates {interpol_coords} \
has repetitive points")
        exit()


    "Any coordinates correspond to the same value? e.g. [Fe/H] and A(Fe) "
    for k in interpol_coords:
        for k1 in interpol_coords:
            if k != k1:
                diff = 100 * ( np.abs( data[k] - data[k1]) ) / data[k]
                if np.max(diff) < 5:
                    print(f"Grid {dataLabel} is only {np.max(diff)} % different \
in parameters {k} and {k1}")
                    exit()

    return


def NDinterpolate_MA(all_par, interpol_par):

    preInterpolationTests(all_par, interpol_par, dataLabel='model_atm')

    " Normalise the coordinates of the grid "
    points = []
    norm_coord = {}
    for k in interpol_par:
            points.append(all_par[k] / max(all_par[k]) )
            norm_coord.update( { k :  max(all_par[k])} )
    points = np.array(points).T

    "Create the function that interpolates model atmospheres structure"
    values = all_par['structure']

    interp_f = LinearNDInterpolator(points, values)

    return interp_f, norm_coord


def NDinterpolate_NLTE_grid(nlte_data, interpol_coords):

    preInterpolationTests(nlte_data, interpol_coords, dataLabel='NLTE')

    " Normalise the coordinates of the grid "
    points = []
    norm_coord = {}
    for k in interpol_coords:
        points.append(nlte_data[k] / max(nlte_data[k]) )
        norm_coord.update( { k :  max(nlte_data[k])} )
    points = np.array(points).T

    values = nlte_data['depart']
    interp_f = LinearNDInterpolator(points, values)

    return interp_f, norm_coord


if __name__ == '__main__':
    exit(0)
