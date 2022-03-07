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
from multiprocessing import Pool, Manager
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

    if set.nnode * set.ncpu > set.inputParams['count']:
        set.nnode = 1
        set.ncpu = set.inputParams['count']
        print(f"Requested more CPUs than jobs. \
Will use {set.nnode} node and {set.ncpu} CPUs")

    

    ind = np.arange(set.inputParams['count'])
    args = [ [set, ind[i::set.ncpu]] for i in range(set.ncpu)]

    with Pool(processes=set.ncpu) as pool:
        pool.map(parallel_worker, args )
