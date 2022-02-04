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

def print_sep():
    print(50*'-')

def mkdir(s):
    if os.path.isdir(s):
        shutil.rmtree(s)
    os.mkdir(s)

def prep_babsma_input(ts_in, atmos, Z, abund):
    """
    Prepares input config files for babsma
    input:
    ts_in (dictionary): contains some of the input values for bsyn and babsma config files
                        property of 'setup' object
    atmos (object):     object of model_atmosphere class, init before passing to this function
    """

    if 'SEGMENTSFILE' in ts_in.keys():
        babsma_conf =  F""" \
'SEGMENTSFILE:' '{ts_in['SEGMENTSFILE']}'
"""
    else:
        babsma_conf = ""

    babsma_conf = babsma_conf + F""" \
'PURE-LTE  :'  '{ts_in['PURE-LTE']}'
'LAMBDA_MIN:'    '{ts_in['LAMBDA_MIN']:.3f}'
'LAMBDA_MAX:'    '{ts_in['LAMBDA_MAX']:.3f}'
'LAMBDA_STEP:'   '{ts_in['LAMBDA_STEP']:.3f}'
'MODELINPUT:'    '{atmos.path}'
'MARCS-FILE:' '{ts_in['MARCS-FILE']}'
'MODELOPAC:' '{ts_in['MODELOPAC']}'
'NLTEINFOFILE:' '{ts_in['NLTEINFOFILE']}'
'RESULTFILE :'    '{ts_in['RESULTFILE']}'
'METALLICITY:'    '{atmos.feh:.3f}'
'ALPHA/Fe   :'    '{atmos.alpha}'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'XIFIX:' 'T'
1.00
    """
    return babsma_conf

def prep_bsyn_input(ts_in, atmos, Z, abund):
    """
    Prepares input config files for babsma and bsyn TS routines
    input:
    ts_in (dictionary): contains some of the input values for bsyn and babsma config files
                        property of 'setup' object
    atmos (object):     object of model_atmosphere class, init before passing to this function
    """
    if 'SEGMENTSFILE' in ts_in.keys():
        bsyn_config =  F""" \
'SEGMENTSFILE:' '{ts_in['SEGMENTSFILE']}'
"""
    else:
        bsyn_config = ""

    bsyn_config = bsyn_config + F""" \
'PURE-LTE  :'  '{ts_in['PURE-LTE']}'
'NLTE :'          '{ts_in['NLTE']}'
'LAMBDA_MIN:'    '{ts_in['LAMBDA_MIN']:.3f}'
'LAMBDA_MAX:'    '{ts_in['LAMBDA_MAX']:.3f}'
'LAMBDA_STEP:'   '{ts_in['LAMBDA_STEP']:.3f}'
'INTENSITY/FLUX:' 'Flux'
'COS(THETA)    :' '1.00'
'MODELINPUT:'    '{atmos.path}'
'MARCS-FILE:' '{ts_in['MARCS-FILE']}'
'NLTEINFOFILE:' '{ts_in['NLTEINFOFILE']}'
'MODELOPAC:'        '{ts_in['MODELOPAC']}'
'RESULTFILE :'    '{ts_in['RESULTFILE']}'
'METALLICITY:'    '{atmos.feh:.3f}'
'ALPHA/Fe   :'    '{atmos.alpha}'
'HELIUM     :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
{abund:5.3f} {Z}
'ISOTOPES : ' '0'
'NFILES   :' '{ts_in['NFILES']}'
{ts_in['LINELIST']}
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
"""

    return bsyn_config

def prep_nlte_input(set, atmos, abund):
    """ Prepare NLTE file for TS """
    depart_file = grid_to_ts(set.depart_grid, set.aux_file, atmos, abund)

    """ Make an input file for TS, NLTE mode """
    with open(set.cwd + '/LTE_NLTE.dat', 'w') as nlte_info_file:
        set.ts_input['NLTEINFOFILE'] = set.cwd + '/LTE_NLTE.dat'
        nlte_info_file.write('# created on \n')
        nlte_info_file.write('# path for model atom files ! this comment line has to be here !\n')
        nlte_info_file.write(F"{set.model_atom_path} \n")

        nlte_info_file.write('# path for departure files ! this comment line has to be here !\n')
        # input NLTE departure file will be written in the cwd, so set path to:
        nlte_info_file.write(F"{set.cwd}/ \n")
        nlte_info_file.write('# atomic (non)LTE setup \n')
        nlte_info_file.write(F"{set.element_z}  '{set.element}'  'nlte' 'atom.{set.model_atom_id}'  '{depart_file}' 'ascii' ")

    return set


def compute_with_ts(set, atmos, abund, routine='clean'):

    """
    very quick and very dirty solution below, need to be changed
    ## TODO: !!!

    Run TS (bsyn, babsma) with prepared ts_input
    Prepare NLTE if needed
    Read spectrum, return
    """

    if routine.lower() == 'babsma':
        set.ts_input['MODELOPAC'] = set.ts_root + F"/opac{set.ts_input['LAMBDA_MIN']}_{set.ts_input['LAMBDA_MAX']}_AA_{atmos.id}_A_{set.element}"

        """ Make input files for TS rouines """
        babsma_conf = prep_babsma_input(set.ts_input, atmos, set.element_z, abund)

        os.chdir(set.ts_root)
        """ Run babsma """
        if set.debug:
            print("babsma_lu:")
            time0 = time.time()
        pr = subprocess.Popen(['./exec/babsma_lu'], stdin=subprocess.PIPE, \
            stdout=open(set.cwd + '/babsma.log', 'w'), stderr=subprocess.STDOUT )
        pr.stdin.write(bytes(babsma_conf, 'utf-8'))
        pr.communicate()
        pr.wait()
        if set.debug:
            t = time.time()-time0
            set.scaling_outputFile.write(F"{atmos.id} {abund} babsma {t:.2f} \n")
            print(F"{t} seconds")
        os.chdir(set.cwd)
        # return


    if routine.lower() == 'bsyn':
        """ How should we name output spectrum? """
        set.ts_input['RESULTFILE'] = set.ts_root + F"/spec_{set.ts_input['LAMBDA_MIN']}_{set.ts_input['LAMBDA_MAX']}_AA_{atmos.id}_A_{set.element}_{abund:4.4f}"
        # without a path
        set.ts_input['RESULTFILE_filename'] = F"/spec_{set.ts_input['LAMBDA_MIN']}_{set.ts_input['LAMBDA_MAX']}_AA_{atmos.id}_A_{set.element}_{abund:4.4f}"

        bsyn_conf = prep_bsyn_input(set.ts_input, atmos, abund, set.element_z)

        os.chdir(set.ts_root)

        """ Run bsyn """
        if set.debug:
            print("bsyn_lu:")
            time0 = time.time()
        pr = subprocess.Popen(['./exec/bsyn_lu'], stdin=subprocess.PIPE, \
            stdout=open(set.cwd + '/bsyn.log', 'w'), stderr=subprocess.STDOUT )
        pr.stdin.write(bytes(bsyn_conf, 'utf-8'))
        pr.communicate()
        pr.wait()
        if set.debug:
            t = time.time()-time0
            set.scaling_outputFile.write(F"{atmos.id} {abund} bsyn {t:.2f} \n")
            print(F"{t} seconds")

        os.chdir(set.cwd)

        """ Move output spectra to the common folder"""
        shutil.move(set.ts_input['RESULTFILE'], \
            set.output_dir  + '/' + set.ts_input['RESULTFILE_filename'] )

        """ Additionally to saving the spectrum file, read it and add to common dictionary """
        w, f, abs_f = np.loadtxt(set.output_dir + '/' + set.ts_input['RESULTFILE_filename'],  unpack=True)


        # """ EXCLUDE REPEATING WAVELENGTH POINTS (can happen 'cause of segments')"""
        # _, un_ind = np.unique(w, return_index=True)
        # w, f, abs_f =  zip(*sorted(zip(w[un_ind], f[un_ind], abs_f[un_ind])))
        # w, f, abs_f = np.array(w), np.array(f), np.array(abs_f)



        return set, w, f, abs_f

    if routine.lower() == 'clean':
        """ Remove temporary file created by babsma, to avoid their accidental misusage"""
        os.remove(set.ts_input['MODELOPAC'])
        set.ts_input['RESULTFILE'] = ''
        set.ts_input['RESULTFILE_filename'] = ''
        set.ts_input['MODELOPAC'] = ''
        # return


"""
A mini wrapper to operate TurboSpectrum
Started within PLATO Solar project
"""

if __name__ == '__main__':
    if len(argv) > 1:
        conf_file = argv[1]
    else:
        conf_file = './config.txt'
    set = setup(file = conf_file)
    exit(0)

    """ Make directory to save output spectra """
    today = datetime.date.today().strftime("%b-%d-%Y")
    set.output_dir = set.cwd + F"/spectra-{set.element}-{today}/"
    mkdir(set.output_dir)
    """ In addition dump them all to a dictionary"""
    spectra_all = { }

    """ Print scaling times"""
    if set.debug:
        set.scaling_outputFile = open(set.cwd + '/scaling.dat', 'w')
        set.scaling_outputFile.write("# atmosphere, A(X), TS routine, run time \n")

    """ Loop over requested model atmospheres """
    for atm_file in set.atmos_list:

        """ First read model atmosphere for parameteres """
        atmos = model_atmosphere()
        atmos.read(file = atm_file, format=set.atmos_format)

        atmos.path = atm_file

        spectra_all.update( { atmos.id :  { 'atmos_obj':atmos, 'abund':[], 'spec':[] } } )

        print_sep()
        print(F"Model atmosphere: {atmos.id}")

        """ Compute model atmosphere opacity """
        compute_with_ts(set, atmos, np.nan, routine='babsma')


        for abund in set.abund_list:
            """ Scale abundance wuth [Fe/H]"""
            abund = abund + atmos.feh

            print(F"A({set.element})={abund:4.2f}")

            """ Prepare NLTE input if requested """
            if set.nlte:
                set = prep_nlte_input(set, atmos, abund)

            """ Run line synthesis """
            set, w, f, abs_f = compute_with_ts(set, atmos, abund, routine='bsyn')

            """ Save to common dictionary """
            spectra_all[atmos.id]['abund'].append(abund)
            spectra_all[atmos.id]['spec'].append( [w, f, abs_f] )


    compute_with_ts(set, atmos, np.nan, routine='clean')


    with open(set.output_dir + F"/spectra_TS_{today}.pkl", 'wb') as f:
            pickle.dump(spectra_all, f)

    """Remove intermediate TS files"""
    mods = glob.glob(set.ts_root + '/*.mod')
    for fle in mods:
        os.remove(fle)

    if set.debug:
        set.scaling_outputFile.close()




    exit(0)
