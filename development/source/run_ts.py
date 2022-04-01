# external
import os
from sys import argv
import shutil
import subprocess
import numpy as np
import pickle
import glob
import time
# local
import convolve
from atmos_package import read_atmos_marcs, model_atmosphere
from read_nlte import grid_to_ts, write_departures_forTS

def mkdir(s):
    if os.path.isdir(s):
        shutil.rmtree(s)
    os.mkdir(s)

def compute_babsma(set, atmos):
    """
    input:
    ts_in (dictionary): contains some of the input values for bsyn and babsma config files
                        property of 'setup' object
    atmos (object):     object of model_atmosphere class, init before passing to this function
    """

    modelOpacFile = set.ts_root + F"/opac{set.ts_input['LAMBDA_MIN']}_{set.ts_input['LAMBDA_MAX']}_AA_{atmos.id}"

    babsma_conf = F""" \
'LAMBDA_MIN:'    '{set.ts_input['LAMBDA_MIN']:.3f}'
'LAMBDA_MAX:'    '{set.ts_input['LAMBDA_MAX']:.3f}'
'LAMBDA_STEP:'   '{set.ts_input['LAMBDA_STEP']:.3f}'
'MODELINPUT:'    '{atmos.path}'
'MARCS-FILE:' '{set.ts_input['MARCS-FILE']}'
'MODELOPAC:' '{modelOpacFile}'
'METALLICITY:'    '{atmos.feh:.3f}'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
    """

    """ Run babsma """
    time0 = time.time()

    os.chdir(set.ts_root)
    pr = subprocess.Popen(['./exec/babsma_lu'], stdin=subprocess.PIPE, \
        stdout=open(set.cwd + '/babsma.log', 'w'), stderr=subprocess.STDOUT )
    pr.stdin.write(bytes(babsma_conf, 'utf-8'))
    pr.communicate()
    pr.wait()
    os.chdir(set.cwd)
    if set.debug:
        print(F"babsma: {time.time()-time0} seconds")

    return modelOpacFile

def compute_bsyn(set, ind, atmos, modelOpacFile, specResultFile, nlteInfoFile=None):
    """
    input:
    atmos (object):     object of model_atmosphere class, init before passing to this function
    """
    bsyn_config = F""" \
'NLTE :'          '{set.ts_input['NLTE']}'
'LAMBDA_MIN:'    '{set.ts_input['LAMBDA_MIN']:.3f}'
'LAMBDA_MAX:'    '{set.ts_input['LAMBDA_MAX']:.3f}'
'LAMBDA_STEP:'   '{set.ts_input['LAMBDA_STEP']:.3f}'
'INTENSITY/FLUX:' 'Flux'
'MARCS-FILE:' '{set.ts_input['MARCS-FILE']}'
'MODELOPAC:'        '{modelOpacFile}'
'RESULTFILE :'    '{specResultFile}'
'HELIUM     :'    '0.00'
'NFILES   :' '{set.ts_input['NFILES']}'
{set.ts_input['LINELIST']}
"""
    if atmos.spherical:
        bsyn_config = bsyn_config + f"""\
'SPHERICAL:'  '.true.'
"""
    else:
        bsyn_config = bsyn_config + f"""\
'SPHERICAL:'  '.false.'
  30
  300.00
  15
  1.30
"""
    if not isinstance(nlteInfoFile,  type(None)):
        bsyn_config = bsyn_config + f"'NLTEINFOFILE:' '{nlteInfoFile}' \n"

    bsyn_config = bsyn_config +\
            f"'INDIVIDUAL ABUNDANCES:'   '{len(set.inputParams['elements'])}' \n"
    for el in set.inputParams['elements']:
        bsyn_config = bsyn_config + f" {set.inputParams['elements'][el]['Z']:.0f} {set.inputParams['elements'][el]['abund'][ind]:5.3f} \n"
# TODO: spherical models???

    """ Run bsyn """
    time0 = time.time()
    os.chdir(set.ts_root)
    pr = subprocess.Popen(['./exec/bsyn_lu'], stdin=subprocess.PIPE, \
        stdout=open(set.cwd + '/bsyn.log', 'w'), stderr=subprocess.STDOUT )
    pr.stdin.write(bytes(bsyn_config, 'utf-8'))
    pr.communicate()
    pr.wait()
    os.chdir(set.cwd)
    if set.debug:
        print(F"bsyn: {time.time()-time0} seconds")

def create_NlteInfoFile(filePath, set, i):
    with open(filePath, 'w') as nlte_info_file:
        nlte_info_file.write('# created on \n')
        nlte_info_file.write('# path for model atom files ! this comment line has to be here !\n')
        nlte_info_file.write(F"{set.modelAtomsPath} \n")

        nlte_info_file.write('# path for departure files ! this comment line has to be here !\n')
        # input NLTE departure file will be written in the cwd, so set path to:
        nlte_info_file.write(F" \n")
        nlte_info_file.write('# atomic (non)LTE setup \n')
        for el in set.inputParams['elements'].values():
            if el.nlte:
                if not isinstance(depart_file, type(None)):
                    model_atom_id = el.modelAtom.split('/')[-1]
                    nlte_info_file.write(F"{el.Z}  '{el.ID}'  'nlte' '{model_atom_id}'  '{el.departFiles[i]}' 'ascii' \n")
                else:
                    if set.debug:
                        print(f"departure file doesn't exist for {el.ID} \
at A({el.ID}) = {el.abund[i]} (i = {i}). No spectrum will be computed.")
                    return False
            else:
                nlte_info_file.write(F"{el.Z}  '{el.ID}'  'lte' ' '  ' ' 'ascii' \n")
        return True

def parallel_worker(arg):
    """
    Run TS on a subset of input parameters (== ind)
    """
    set, ind = arg
    tempDir = f"{set.cwd}/job{min(ind)}_{max(ind)}/"
    mkdir(tempDir)
    print(f"Working on input parameters {ind}...")
    for i in ind:
        # create model atmosphere and run babsma on it
        atmos = model_atmosphere()
        if not isinstance(set.inputParams['modelAtmInterpol'][i], type(None)):
            atmos.depth_scale, atmos.temp, atmos.ne, atmos.vturb = \
                set.inputParams['modelAtmInterpol'][i]
            set.inputParams['modelAtmInterpol'][i] = None
            atmos.temp, atmos.ne = 10**(atmos.temp), 10**(atmos.ne)
            atmos.depth_scale_type = 'TAU500'
            atmos.feh, atmos.logg = set.inputParams['feh'][i], set.inputParams['logg'][i]
            atmos.spherical = False
            atmos.id = f"interpol_{i:05d}"
            atmos.path = f"{tempDir}/atmos.{atmos.id}"
            atmos.write(atmos.path, format = 'ts')

            """ Compute model atmosphere opacity """
            modelOpacFile = compute_babsma(set, atmos)

            """ Compute the spectrum """
            specResultFile = f"{tempDir}/spec"
            for el in set.inputParams['elements'].values():
                specResultFile = specResultFile + f"_{el.ID}{el.abund[i]}"
            if set.nlte:
                specResultFile = specResultFile + '_NLTE'
            else:
                specResultFile = specResultFile + '_LTE'

            if set.nlte:
                nlteInfoFile   = f"{tempDir}/NLTEinfoFile.txt"
                departFilesExist = create_NlteInfoFile(nlteInfoFile, set, i)
                if departFilesExist:
                    compute_bsyn(set, i, atmos, modelOpacFile, specResultFile, nlteInfoFile)
                    if os.path.isfile(specResultFile):
                        shutil.move(specResultFile, f"{set.spectraDir}/{specResultFile.split('/')[-1]}" )

            os.remove(atmos.path)
            os.remove(modelOpacFile)
            os.remove(modelOpacFile+'.mod')
    return set

if __name__ == '__main__':
    if len(argv) > 1:
        conf_file = argv[1]
    else:
        print("Usage: ./run_ts.py ./configFile.txt")
        exit()
    #set = setup(file = conf_file)
    exit(0)
