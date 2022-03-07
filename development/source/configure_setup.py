import numpy as np
import os
import shutil
from sys import argv, exit
import datetime
import glob
from scipy.spatial import Delaunay
# local
from model_atm_interpolation import get_all_ma_parameters, NDinterpolateGrid, tensorInterpN
from read_nlte import read_fullNLTE_grid
from atmos_package import model_atmosphere
from run_ts import write_departures_forTS
import cProfile
import pstats

def in_hull(p, hull):
   return hull.find_simplex(p) >= 0


"""
Reading the config file and preparing for the computations
"""
def mkdir(s):
    if os.path.isdir(s):
        shutil.rmtree(s)
    os.mkdir(s)

def atomicZ(el):
    if os.path.isfile('./atomic_numbers.dat'):
        el_z = np.loadtxt('./atomic_numbers.dat', usecols=(0))
        el_id = np.loadtxt('./atomic_numbers.dat', usecols=(1), dtype=str)
    else:
        print("Can not find './atomic_numbers.dat' file. Stopped.")
        exit(1)
    for i in range(len(el_id)):
        if el.lower() == el_id[i].lower():
            return el_z[i]

    self.MAhull = np.array([ modelAtmGrid[k] for k in interpolCoords ]).T


def read_random_input_parameters(file):
    """
    Read input file listing requested labels

    First four columns are Teff, logg, Vturb, Fe, then all the elements
    Example:
    -----
    Teff logg Vturb FeH Fe H O # Ba
    5535  2.91  0.0  -1.03  6.470  12.0   9.610 # 2.24
    7245  5.74  0.0  -0.50  7.000  12.0   8.009 # -2.2
    ....
    -----
    """
    data =[ l.split('#')[0] for l in open(file, 'r').readlines() \
                            if not (l.startswith('#') or l.strip()=='') ]
    elements = data[0].replace("'","").split()[4:]

    values =  [ l.split() for l in data[1:] ]
    values = np.array(values).astype(float)

    input_par = {'teff':values[:, 0], 'logg':values[:, 1], 'vturb':values[:, 2], 'feh':values[:,3], \
                'elements' : {
                            elements[i].capitalize() : {'abund': values[:, i+4], 'nlte':False, 'Z' : atomicZ(elements[i])} \
                                                for i in range(len(elements))
                                }
                }
    input_par.update({'count' : len(input_par['teff'])})
    if 'Fe' not in input_par['elements']:
        print(f"Warning: input contains [Fe/H], but no A(Fe)")
    absAbundCheck = np.array([ input_par['elements'][el]['abund'] / 12. for el in input_par['elements'] ])
    if (absAbundCheck < 0.1).any():
        print(f"Warning: abundances must be supplied relative to H, on log12 scale. Please double check input file '{file}'")



    return input_par

class setup(object):
    def __init__(self, file='./config.txt'):
        if 'cwd' not in self.__dict__.keys():
            self.cwd = f"{os.getcwd()}/"
        self.debug = 0
        self.ncpu  = 1
        self.nnode  = 1


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
        if 'nlte_config' in self.__dict__:
            d = {}
            for l in self.nlte_config:
                l = l.replace('[','').replace(']','').replace("'","")
                el, files = l.split(':')[0].strip(), l.split(':')[-1].split(',')
                files = [ f.strip() for f in files ]
                files = [ f.replace('./', self.cwd) if f.startswith('./') else f  for f in files]
                d.update({el.capitalize() : {'nlteGrid' : files[0], 'nlteAux' : files[1], 'modelAtom' : files[2] }})
            self.nlte_config = d

            "TS needs to access model atoms from the same path for all elements"
            if 'modelAtomsPath' not in self.__dict__.keys():
                self.modelAtomsPath = f"{self.cwd}/modelAtoms_links/"
                mkdir(self.modelAtomsPath)
                for el in self.nlte_config:
                    dst = self.modelAtomsPath + self.nlte_config[el]['modelAtom'].split('/')[-1]
                    os.symlink(self.nlte_config[el]['modelAtom'], dst )

        if 'inputParams_file' in self.__dict__:
            self.inputParams = read_random_input_parameters(self.inputParams_file)
        else:
            print("Missing file with input parameters: inputParams_file")
            exit()

        if 'nlte_config' not in self.__dict__ or len(self.nlte_config) == 0:
            print(f"{50*'*'}\n Note: all elements will be computed in LTE!\n \
To set up NLTE, use 'nlte_config' flag\n {50*'*'}")
        else:
            for el in self.nlte_config:
                self.inputParams['elements'][el]['nlte'] = True
                for k in self.nlte_config[el]:
                    self.inputParams['elements'][el].update({
                                                    k : self.nlte_config[el][k]
                                                            })

        self.nlte = False
        for el in self.inputParams['elements']:
            if self.inputParams['elements'][el]['nlte']:
                self.nlte = True
                break

        today = datetime.date.today().strftime("%b-%d-%Y")
        self.spectraDir = self.cwd + f"/spectra-{today}/"
        if not os.path.isdir(self.spectraDir):
            os.mkdir(self.spectraDir)


        "Temporary directories for NLTE files"
        for el in self.inputParams['elements']:
            if self.inputParams['elements'][el]['nlte']:
                path = self.cwd + f"/{el}_nlteDepFiles/"
                mkdir(path)
                self.inputParams['elements'][el].update({'dirNLTE':path})

        atmDepthScale, interpolCoords = self.prepInterpolation_MA()
        # self.interpolateAllPoints_MA()
        exit()
        self.interpolator['modelAtm'] = None
        # have to work with one nlte grid at a time to avoid memory overflow
        if self.nlte:
            self.interpolator.update({'NLTE':{}})
            for el in self.inputParams['elements']:
                if self.inputParams['elements'][el]['nlte']:
                    self.prepInterpolation_NLTE(el,interpolCoords, rescale = True, depthScale = atmDepthScale)
                    self.interpolateAllPoints_NLTE(el)
                    self.interpolator['NLTE'][el] = None

        self.createTSinputFlags()



    def prepInterpolation_MA(self):
        """
        Read grid of model atmospheres and NLTE grids of departures
        and prepare interpolating functions
        Store for future use
        """
        self.interpolator = {'modelAtm' : None}

        " Over which parameters (== coordinates) to interpolate?"
        interpolCoords = ['teff', 'logg', 'feh'] # order should match input file!
        if 'vturb' in self.inputParams:
            interpolCoords.append('vturb')

        "Model atmosphere grid"
        if self.debug: print("preparing model atmosphere interpolator...")
        modelAtmGrid, atmDepthScale = get_all_ma_parameters(self.atmos_path, \
                                        format = self.atmos_format, debug=self.debug)

        for i in range(self.inputParams['count']):
            q = (c in modelAtmGrid[k] for k in interpolCoords )
            print(f"q={q}")
            qi = (c in self.inputParams[k][i] for k in interpolCoords )
            v = modelAtmGrid['structure']
            value = tensorInterpN(q, v, qi, method='cubic')
            print(value)

        # interpFunction, normalisedCoord = NDinterpolateGrid(modelAtmGrid, interpolCoords, \
                                        # valueKey='structure', dataLabel = 'model atmosphere grid' )



        # """
        # Create hull object to test whether which of the requested points
        # are within the original grid
        # Interpolation outside of hull returns NaNs, therefore skip those points
        # """
        # hull = Delaunay(np.array([ modelAtmGrid[k] / normalisedCoord[k] for k in interpolCoords ]).T)
        #
        # self.interpolator['modelAtm'] = {'interpFunction' : interpFunction, \
        #                                 'normCoord' : normalisedCoord, \
        #                                 'hull': hull}
        # return atmDepthScale, interpolCoords

    def prepInterpolation_NLTE(self, el, interpolCoords, rescale = True, depthScale = None):
        profiler = cProfile.Profile()
        profiler.enable()
        "NLTE grids"
        if self.debug: print(f"reading grid {self.inputParams['elements'][el]['nlteGrid']}")
        " 0th element is tau, 1th-Nth are departures for N levels "
        nlteData = read_fullNLTE_grid( self.inputParams['elements'][el]['nlteGrid'], \
                                    self.inputParams['elements'][el]['nlteAux'], \
                                    rescale=True, depthScale = depthScale )
        if self.debug: print('done reading')
        " interpolate over abundance? "
        interpolCoords_el = interpolCoords.copy()
        if min(nlteData['abund']) == max(nlteData['abund']):
            pass
        elif len(np.unique(nlteData['feh'])) == len(np.unique(nlteData['abund'])): # it is Fe
            coords = [c for c in interpolCoords_el if c!='feh']
            interpolCoords_el = coords
            interpolCoords_el.append('abund')
            if self.debug: print(f"for Fe grid interpolating over abundance rather than feh")
        else:
            if self.debug: print(f"included interpolation over abundance of {el}")
            interpolCoords_el.append('abund')

        #if self.debug: print("starting interpolation")
        #interpFunction, normalisedCoord  = NDinterpolateGrid(nlteData, interpolCoords_el,\
        #                                        valueKey='depart', dataLabel=f"NLTE grid {el}")
        if self.debug: print(f"Interpolated, creating hull")
        hull = Delaunay(np.array([ nlteData[k] /normalisedCoord[k] for k in interpolCoords_el ]).T)
        self.interpolator['NLTE'].update( { el: {
                                            'interpFunction' : interpFunction, \
                                             'normCoord' : normalisedCoord, 'hull':hull }
                                         } )
        nlteData = None
        if self.debug: print(f"{el} done")
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumulative')
        stats.print_stats()


    def interpolateAllPoints_MA(self):
        """
        Python parallelisation libraries can not send more than X Gb of data between processes
        To avoid that, interpolation at each requested point is done before the start of computations
        """
        if self.debug: print(f"Interpolating to each of {self.inputParams['count']} requested points...")

        "Model atmosphere grid"
        self.inputParams.update({'modelAtmInterpol' : np.full(self.inputParams['count'], None) })

        countOutsideHull = 0
        for i in range(self.inputParams['count']):
            point = [ self.inputParams[k][i] / self.interpolator['modelAtm']['normCoord'][k] \
                    for k in self.interpolator['modelAtm']['normCoord'] ]


            if not in_hull(np.array(point).T, self.interpolator['modelAtm']['hull']):
                countOutsideHull += 1
            else:
                values =  self.interpolator['modelAtm']['interpFunction'](point)[0]
                self.inputParams['modelAtmInterpol'][i] = values
        if countOutsideHull > 0 and self.debug:
            print(f"{countOutsideHull}/{self.inputParams['count']}requested \
points are outside of the model atmosphere grid. No computations will be done")

    def interpolateAllPoints_NLTE(self, el):
        " NLTE grids "
        self.inputParams['elements'][el].update({
                        'departFile' : np.full(self.inputParams['count'], None)
                                                })
        profiler = cProfile.Profile()
        profiler.enable()
        for i in range(self.inputParams['count']):
            depart = None
            if not isinstance(self.inputParams['modelAtmInterpol'][i], type(None)):
                point = [ self.inputParams[k][i] / self.interpolator['NLTE'][el]['normCoord'][k] \
                        for k in self.interpolator['NLTE'][el]['normCoord'] if k !='abund']
                if 'abund' in self.interpolator['NLTE'][el]['normCoord']:
                    abund = self.inputParams['elements'][el]['abund'][i]
                    point.append(abund)

                fallsIntoHull = in_hull(np.array(point).T, self.interpolator['NLTE'][el]['hull'])
                if not fallsIntoHull:
                    " Trying to find a closer abundance value "
                    if 'abund' in self.interpolator['NLTE'][el]['normCoord']:
                        shift = 0.0
                        sign = [-1 if abund > np.mean(self.inputParams['elements'][el]['abund']) else 1][0]

                        while not fallsIntoHull and shift <= 0.5: # dont walk futher than 0.5 dex from the requested abundance
                            abund = abund + sign * 0.1
                            point = [ self.inputParams[k][i] / self.interpolator['NLTE'][el]['normCoord'][k] \
                                    for k in self.interpolator['NLTE'][el]['normCoord'] if k != 'abund' ]
                            point.append(abund)
                            fallsIntoHull = in_hull(np.array(point).T, self.interpolator['NLTE'][el]['hull'])
                            shift += 0.1
                            if self.debug: print(f"{el} at abund={abund:.2f} is outside of hull. trying abund={abund + sign * shift:.2f}")
                        if fallsIntoHull:
                            depart = self.interpolator['NLTE'][el]['interpFunction'](point)[0]

                    else: print(f"WARNING: point outside of hull, but abundance is not a free parameter for {el}")

                else:
                    depart = self.interpolator['NLTE'][el]['interpFunction'](point)[0]

                if not isinstance(depart, type(None)):
                    tau = depart[0]
                    depart_coef = depart[1:]
                    departFile = f"{self.inputParams['elements'][el]['dirNLTE']}/depart_{el}_{i}"
                    write_departures_forTS(departFile, tau, depart_coef, abund)
                    self.inputParams['elements'][el]['departFile'][i] = departFile
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumulative')
        stats.print_stats()
        exit(0)

    def createTSinputFlags(self):
        self.ts_input = { 'PURE-LTE':'.false.', 'MARCS-FILE':'.false.', 'NLTE':'.false.',\
        'NLTEINFOFILE':'', 'LAMBDA_MIN':4000, 'LAMBDA_MAX':9000, 'LAMBDA_STEP':0.05,\
         'MODELOPAC':'./OPAC', 'RESULTFILE':'' }


        """ At what wavelenght range to compute a spectrum? """
        self.lam_start, self.lam_end = min(self.lam_end, self.lam_start), \
                                    max(self.lam_end, self.lam_start)
        self.wave_step = np.mean([self.lam_start, self.lam_end]) / self.resolution
        self.ts_input['LAMBDA_MIN'] = self.lam_start
        self.ts_input['LAMBDA_MAX'] = self.lam_end
        self.ts_input['LAMBDA_STEP'] = self.wave_step


        """ Linelists """
        if type(self.linelist) == np.array or type(self.linelist) == np.ndarray:
            pass
        elif type(self.linelist) == str:
            self.linelist = np.array([self.linelist])
        else:
            print(f"Can not understand the 'linelist' flag: {self.linelist}")
            exit(1)
        llFormatted = []
        for path in self.linelist:
            if '*' in path:
                llFormatted.extend( glob.glob(path) )
            else:
                llFormatted.append(path)
        self.linelist = llFormatted
        print(f"Linelist(s) will be read from: {' ; '.join(str(x) for x in self.linelist)}")

        self.ts_input['NFILES'] = len(self.linelist)
        self.ts_input['LINELIST'] = '\n'.join(self.linelist)


        "Any element in NLTE?"
        if self.nlte:
            self.ts_input['NLTE'] = '.true.'
