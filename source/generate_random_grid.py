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
from multiprocessing import Pool
# local
from configure_setup import setup
from atmos_package import read_atmos_marcs, model_atmosphere
from read_nlte import grid_to_ts
from run_ts import parallel_worker



if __name__ == '__main__':
    if len(argv) > 1:
        conf_file = argv[1]
    else:
        print("Usage: ./run_ts.py ./configFile.txt")
        exit()
    set = setup(file = conf_file)
<<<<<<< HEAD
    if set.ncpu > set.inputParams['count']:
        set.ncpu = set.inputParams['count']
        print(f"Requested more CPUs than jobs. Will use {set.ncpu}")

    ind = np.arange(set.inputParams['count'])
    args = []
    args = [ [set, ind[i::ncpu]] for i in range(set.ncpu)]
    print(args)

    with Pool(processes=set.ncpu) as pool:
            pool.map( parallel_worker, args )

=======
>>>>>>> 6e32fcf3e867796e043f6e2b65ce179e880a9bc6
    ind = [0]
    set = parallel_worker(set, ind)
