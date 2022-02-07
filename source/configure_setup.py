import numpy as np
import os
from sys import argv, exit
from model_atm_interpolation import get_all_ma_parameters, NDinterpolate_MA, NDinterpolate_NLTE_grid
from read_nlte import read_full_grid
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
    elements = data[0].replace("'","").split()[3:]

    values =  [ l.split() for l in data[1:] ]
    values = np.array(values).astype(float)

    # print(values)
    input_par = {'teff':values[:, 0], 'logg':values[:, 1], 'vturb':values[:, 2], \
                'elements' : {
                            elements[i].capitalize() : {'abund': values[:, i+3], 'nlte':False} \
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
    def __init__(self, file='./config.txt'):
        self.cwd = os.getcwd()
        self.debug = 0

        "Read all the keys from the config file"
        for line in open(file, 'r').readlines():
            line = line.strip()
            if not line.startswith('#') and len(line)>0:
                if not '+=' in line:
                    k, val = line.split('=')
                    k, val = k.strip(), val.strip()
                    if val.startswith("'") or val.startswith('"'):
                        self.__dict__[k] = val[1:-1]
                    elif val.startswith("["):
                        if '[' in val[1:]:
                            if not k in self.__dict__ or len(self.__dict__[k]) == 0:
                                self.__dict__[k] = []
                            self.__dict__[k].append(val)
                        else:
                            self.__dict__[k] = eval('np.array(' + val + ')')
                    elif '.' in val:
                        self.__dict__[k] = float(val)
                    else:
                        self.__dict__[k] = int(val)
                elif '+=' in line:
                    k, val = line.split('+=')
                    k, val = k.strip(), val.strip()
                    if len(self.__dict__[k]) == 0:
                        self.__dict__[k] = []
                    self.__dict__[k].append(val)

        # transfer list of lines in nlte_config to a dictionary
        d = {}
        for l in self.nlte_config:
            l = l.replace('[','').replace(']','').replace("'","")
            el, files = l.split(':')[0].strip(), l.split(':')[-1].strip().split(',')
            d.update({el.capitalize() : {'nlteGrid' : files[0], 'nlteAux' : files[1], 'modelAtom' : files[2] }})
        self.nlte_config = d


        if not 'nlte_config' in self.__dict__ or len(self.nlte_config) == 0:
            print(f"{50*'*'}\n Warning: all elements will be computed in LTE!\n To set up NLTE, use 'nlte_config' flag\n {50*'*'}")
        if 'inputParams_file' in self.__dict__:
            self.inputParams = read_random_input_parameters(self.inputParams_file)
            for el in self.inputParams['elements']:
                if el in self.nlte_config:
                    self.inputParams['elements'][el]['nlte'] = True
                    for k in self.nlte_config[el]:
                        self.inputParams['elements'][el].update({
                                k : self.nlte_config[el][k]
                                                            })
            self.prepInterpolation()



    def prepInterpolation(self):
        """
        Read grid of model atmospheres and NLTE grids of departures
        and prepare interpolating functions
        Store for future use
        """
        self.interpolator = {
            'modelAtm' : None,
            'NLTE' : {}
        }

        " Over which parameters (aka coordinates) to interpolate?"
        interpolCoords = ['teff', 'logg', 'feh']
        if 'vturb' in self.inputParams:
            interpolCoords.append('vturb')

        "Model atmosphere grid"
#        if self.debug:
#            print("preparing model atmosphere interpolator...")
#        modelAtmGrid= get_all_ma_parameters(self.atmos_path, \
#                                        format = self.atmos_format, debug=self.debug)
#
#
#        interpFunction, normalisedCoord = NDinterpolate_MA(modelAtmGrid, interpolCoords )
#
#        self.interpolator['modelAtm'] = {'interpFunction' : interpFunction, \
#                                                'normCoord' : normalisedCoord}

        "NLTE grids"
        interpolCoords.append('abund')
        for el in self.inputParams['elements']:
            if self.inputParams['elements'][el]['nlte']:
                if self.debug:
                    print(f"preparing interpolator for {el}")

                nlteData = read_full_grid( self.inputParams['elements'][el]['nlteGrid'], \
                                            self.inputParams['elements'][el]['nlteAux'] )
                interpFunction, normalisedCoord  = NDinterpolate_NLTE_grid(nlteData, interpolCoords)
                self.interpolator['NLTE'].update( { el: {
                                                    'interpFunction' : interpFunction, \
                                                     'normCoord' : normalisedCoord}
                                                 } )
