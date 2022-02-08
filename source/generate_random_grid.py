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
# local
import convolve
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
    ind = [0,1,2]
    set = parallel_worker(set, ind)
