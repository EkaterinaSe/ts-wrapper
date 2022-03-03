import numpy as np
import os
from scipy.interpolate import interp1d
from sys import exit

def grid_to_ts(grid_file, aux_file, atmos, abund):

    """ Find model atmosphere in the auxiliarly file """
    aux_atmos = np.loadtxt(aux_file, usecols=(0), dtype=str)
    aux_atmos = np.array([x.replace("'","").replace('"',"") for x in aux_atmos])
    aux_abund, aux_pointer = np.loadtxt(aux_file, usecols=(7,8), unpack=True)

    if atmos.id in aux_atmos:
        pos = np.where(aux_atmos == atmos.id)[0]
        """ It will always pick the first record, that might need to be changed later
        if e.g. Vmic is included """
        pos1 = np.where(np.abs(aux_abund[pos]-abund) == min(np.abs(aux_abund[pos]-abund)))[0][0]
        pointer = int(aux_pointer[pos][pos1])


    else:
        print("Can't find atmos in aux file. Stopped")
        exit(1)

    if np.abs(aux_abund[pos][pos1]-abund) > 0.1:
        print(F"Warning: abundance for NLTE departure coefficients: {aux_abund[pos][pos1]:.3f}")
        print(F"Warning: {aux_abund[pos][pos1]-abund:.3f} dex different from target abundance: {abund:.3f}")

    ndep, nk, depart, tau = read_binary_grid(grid_file, pointer)

    if atmos.depth_scale_type.startswith('T'):
        if atmos.ndep == len(tau):
            diff = np.sqrt(np.sum( (tau - atmos.depth_scale)**2 ) )
            if diff > 0.05:
                print("Seems like depth scale in model atmosphere and departure file differ:")
                print(F"sqrt(sum(d_atm - d_dep)**2) = {diff}. Stopped")
                exit(1)
            else:
                # because TS has rather corrupt way of checking it..
                tau = atmos.depth_scale.copy()
        else:
            print("NDEP in departure grid != NDEP in model atmosphere. Stopped.")
            exit(1)



    depFile = './nlte_dep.dat'
    depart = depart.T
    with open(depFile, 'w') as f:
        for i in range(8):
            f.write('# parameter 1.0 1.0\n')
        f.write(F"{abund:.3f}\n")
        f.write(F"{ndep:.0f}\n")
        f.write(F"{nk:.0f}\n")
        for t in tau:
            f.write(F"{t}\n")
        for i in range(ndep):
            for j in range(nk):
                f.write(F"{depart[i,j]}" + ' ') # separate numbers in lines with blank space
            f.write('\n') # separate lines


    return depFile

def write_departures_forTS(filePath, tau, depart, abund):
    ndep = len(tau)
    nk = len(depart)
    with open(filePath, 'w') as f:
        # these comment lines have to be here for TS
        # I can not help it, someone help me
        for i in range(8):
            f.write('# parameter 1.0 1.0\n')

        f.write(F"{abund:.3f}\n")
        f.write(F"{ndep:.0f}\n")
        f.write(F"{nk:.0f}\n")
        for t in tau:
            f.write(F"{t:15.8E}\n")

        for i in range(ndep):
            f.write( f"{'  '.join(str(depart[j,i]) for j in range(nk))} \n" )


def read_binary_grid(grid_file, pointer=1):
    """ Read a record from NLTE grid """
    with open(grid_file, 'rb') as f:
        # -1 since Python stars at 0
        pointer = pointer - 1

        f.seek(pointer)
        atmosStr = f.readline(500)#.decode('utf-8', 'ignore').strip()
        ndep = int.from_bytes(f.read(4), byteorder='little')
        nk = int.from_bytes(f.read(4), byteorder='little')
        tau  = np.log10(np.fromfile(f, count = ndep, dtype='f8'))
        depart = np.fromfile(f, dtype='f8', count=ndep*nk).reshape(nk, ndep)

    return ndep, nk, depart, tau


def read_fullNLTE_grid(bin_file, aux_file, rescale=False, depthScale=None):
    if rescale and isinstance(depthScale, type(None)):
            print(f"to re-scale NLTE departure coefficient, please supply new depth scale to read_fullNLTE_grid() ")
            exit()
    aux = np.genfromtxt(aux_file, \
    dtype = [('atmos_id', 'str'), ('teff','f8'), ('logg','f8'), ('feh', 'f8'),\
             ('alpha', 'f8'), ('mass', 'f8'), ('vturb', 'f8'), ('abund', 'f8'), \
             ('pointer', 'i8')])

    data = {}
    for k in aux.dtype.names:
        data.update( { k : aux[k] } )

    # read and save each record from the binary file
    data.update( { 'depart' : [] } )
    for p in data['pointer']:
        #print(f"{100*p/data['pointer'][-1]:.2f}")
        ndep, nk, depart, tau = read_binary_grid(bin_file, pointer=p)
        if rescale:
            depart_new = np.full(  shape=(len(depthScale), nk), fill_value=np.nan)
            f_int = interp1d(tau, depart, fill_value='extrapolate')
            depart_new = f_int(depthScale)
        data['depart'].append( np.vstack([depthScale, depart_new]) )



    data['depart'] = np.array( data['depart'] )
    return data

    # TODO: make sure that all departure coefficient are at the same tau scale
    # also, same as model atmospheres... i guess?
