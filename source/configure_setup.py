import numpy as np
import os
# local

"""
Reading the config file and preparing for the computations
"""
def read_elemental_setup(file):
    """
    Read info concerning what chemical elements should be included and how
    """
    elements = {}

    data = np.genfromtxt(file, \
    dtype=[('element', 's'),
            ('nlteFlag', 'i4'), ('modelAtom', 's'), \
            ('nlteBinGrid', 's'), ('nlteAuxFile', 's') ])


    # read which files to use in case of NLTE
    for i in range(len(data['element']))
        el = data['element'][i]
        if data['nlteFlag'][i]:
            elements[el]['nlte'] = True
            elements[el].update( { 'nlteBinGrid':data['nlteBinGrid'][i],\
 'modelAtom':data['modelAtom'][i], 'nlteAuxFile':data['nlteAuxFile'][i] } )
    return elements


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
                elif val.startswith("["):
                    self.__dict__[key] = eval('np.array(' + val + ')')
                elif '.' in val:
                    self.__dict__[key] = float(val)
                else:
                    self.__dict__[key] = int(val)


    if 'elements_input' not in self.__dict__:
        print("Missing input file with elemental abundances, 'elements_input'=''")
        exit()
    else:
        self.elements = read_elemental_setup(self.elements_input)
