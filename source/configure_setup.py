import numpy as np
import os
# local

"""
Reading the config file and preparing for the computations
"""

def read_random_input_parameters(file):
    params = {}

    return params

class setup(object):
    def __init__(self, cnfg_file='./config.txt'):
        self.cwd = os.getcwd()
        self.debug = 0

        "Read all the keys from the config file"
        for line in open(cnfg_file, 'r').readlines():
            line = line.strip()
            if not line.startswith('#') and len(line)>0:
                k, value = line.split('=')
                k, value = k.strip(), value.strip()
                if val.startswith("'") or val.startswith('"'):
                    self.__dict__[k] = val[1:-1]
                elif k == 'nlte_config':
                    val = val.replace('[', '{').replace(']','}')
                    self.__dict__[k] = eval('dict(' + val + ')')
                elif val.startswith("["):
                    self.__dict__[k] = eval('np.array(' + val + ')')
                elif '.' in val:
                    self.__dict__[k] = float(val)
                else:
                    self.__dict__[k] = int(val)
