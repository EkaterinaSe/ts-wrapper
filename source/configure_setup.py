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
    dtype=[('element', 's'), ('abFlag', 'i4'), \
            ('abundStart', 'f4'), ('abundEnd', 'f4'), ('abundStep', 'f4'), \
            ('nlteFlag', 'i4'), ('modelAtom', 's'), \
            ('nlteBinGrid', 's'), ('nlteAuxFile', 's') ])

    # read what abundances are requested
    for i in range(len(data['element']))
        el = data['element'][i]
        elements.update( el : { 'abund' : None, 'nlte':False } )

        abFlag = data['abFlag'][i]
        if abFlag == 0: # solar
        # TODO: not sure how to implement this without hard-coding solar abundances
            print(f"solar abundances not supported, change for {el} in {file}")
            exit()
        elif abFlag == 1: # use new abundance
            abundNew = data['abundStart'][i]
            elements[el]['abund'] = [abundNew]
        elif abFlag == 2: # use an array of abundances
            abStart = data['abStart'][i]
            abEnd = data['abEnd'][i]
            abStep = data['abStep'][i]

            abundArray = np.arange( abStart, abEnd, abStep )
            if abEnd not in abundArray:
                abundArray = np.hstack( [ abundArray, [abEnd] ] )
            elements[el]['abund'] = abundArray

    # check that requested number of abundances is compatible among elements
    abundCount = 1
    for el in elements:
        if len(elements[el]['abund']) != abundCount:
            abundCount = len(elements[el]['abund'])
            for otherEl in elements:
                if otherEl != el:
                    if len(elements[otherEl]['abund']) == 1:
                        elements[otherEl]['abund'] = np.full(abundCount, elements[otherEl]['abund'][0])
                    else:
                        print(f"Requested {abundCount} abundances for {el}, \
which does not match with other elements. Change {file}.")
                        exit()


    # read which files to use in case of NLTE
    for i in range(len(data['element']))
        el = data['element'][i]
        if data['nlteFlag'][i]:
            elements[el]['nlte'] = True
            elements[el].update( { 'nlteBinGrid':data['nlteBinGrid'][i],\
 'modelAtom':data['modelAtom'][i], 'nlteAuxFile':data['nlteAuxFile'][i] } )
    return elements



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
