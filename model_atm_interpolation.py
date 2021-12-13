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


def get_all_models(models_path, format='m1d', debug = False):
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
        pass
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

def interpolate_ma_grid(setup):
    get_all_models(setup.atmos_path, format=setup.atmos_format, debug=setup.debug)
