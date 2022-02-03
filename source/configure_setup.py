import numpy as np
import os
from sys import argv, exit
# local

"""
Reading the config file and preparing for the computations
"""


def read_random_input_parameters(file):
    """
    Read input file listing requested labels

    First four columns are Teff, logg, Vturb, Fe, then all the elements
    Example:
    -----
    Teff logg Vturb Fe H O Mg Ca Mn Ni # Ba
    5535  2.91  0.0  -1.036   0.0   0.810 -0.369  0.663 -0.854 -0.219 # 0.02
    7245  5.74  0.0  -0.493   0.0  -0.589  0.303  0.440  0.595  0.243 # -4.2
    ....
    -----
    """
    data =[ l.split('#')[0] for l in open(file, 'r').readlines() \
                            if not (l.startswith('#') or l.strip()=='') ]
    elements = data[0].split()[3:]

    values =  [ l.split() for l in data[1:] ]
    values = np.array(values).astype(float)

    # print(values)
    input_par = {'teff':values[:, 0], 'logg':values[:, 1], 'vturb':values[:, 2], \
                'elements' : {
                            elements[i].capitalize() : {'abund': values[:, i], 'nlte':False} \
                                                for i in range(len(elements))
                                }
                }
    if 'Fe' not in  input_par['elements']:
        raise Warning('Has to include Fe in the list of input parameters,\
 (defines [Fe/H])')
    else:
        input_par.update( { 'feh' : input_par['elements']['Fe']['abund'] } )


    return input_par

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

        if 'input_params_file' in self.__dict__:
            self.input_params = read_random_input_parameters(self.input_params_file)
            for el in self.input_params['elements']:

                if el in self.nlte_config:
                    self.input_params['elements'][el]['nlte'] = True
                    self.input_params['elements'][el].update({
                            'nlteGrid' : self.nlte_config['nlteGrid'],
                            'nlteAux' : self.nlte_config['nlteAux'],
                            'modelAtom' : self.nlte_config['modelAtom']
                                                            })
