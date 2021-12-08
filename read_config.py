import numpy as np
import os
# local
from observations import *

"""
Setup that will be used to run TS
"""
class setup(object):
    def __init__(self, file='config.txt'):

        self.cwd = os.getcwd()

        """
        Reads specifications for the TS run from a config file
        """
        "Some default keys:"
        self.debug = 0
        self.interpolate = 0

        print('Reading configuration file %s' %(file ) )
        for line in open(file, 'r'):
            line=line.strip()
            if not line.startswith("#") and line != '':

                    key, val = line.split('=')
                    key, val = key.strip(), val.strip()
                    if val.startswith("'") or val.startswith('"'):
                        self.__dict__[key] = val[1:-1]
                    elif val.startswith("["):
                        self.__dict__[key] = eval('np.array(' + val + ')')
                    elif '.' in val:
                        self.__dict__[key] = float(val)
                    else:
                        self.__dict__[key] = int(val)

        if self.use_abund == 1:
            self.abund_list = [self.new_abund]
        elif self.use_abund == 2:
            self.abund_list = np.arange(self.start_abund, self.end_abund, self.step_abund)
            # [start, end) --> [start, end]
            self.abund_list = np.hstack((self.abund_list,  self.end_abund ))

        print("Reading a list of model atmospheres from %s" %( self.use_atmos ))
        self.atmos_list = np.loadtxt(self.use_atmos, ndmin=1, dtype=str, usecols=(0))
        self.atmos_list = [ self.atmos_path + F"/{x}" for x in self.atmos_list ]

        print('Element: %s' %self.element)
        print(F"Use {len(self.abund_list)} abundances")
        print(F"Use {len(self.atmos_list)} model atmospheres")

        if self.nlte == 1:
            self.nlte = True
            self.lte = False

        elif self.nlte == 0:
            self.lte = True
            self.nlte = False

        else:
            print("'nlte' flag unrecognised")
            exit(1)

        if self.nlte:
            if not os.path.isfile(self.depart_grid):
                print(F"Specified departure grid not found: {self.depart_grid}")
            if not os.path.isfile(self.aux_file):
                print(F"Specified auxliarly file not found: {self.aux_file}")


        self.ts_input = { 'PURE-LTE':'', 'NLTE':'', 'MARCS-FILE':'.true.', 'NLTEINFOFILE':'', 'LAMBDA_MIN':0, 'LAMBDA_MAX':0, 'LAMBDA_STEP':0, 'MODELOPAC':'OPAC', 'RESULTFILE':'' }
        if self.lte:
            self.ts_input['PURE-LTE'] = '.true.'
            self.ts_input['NLTE'] = '.false.'

        if self.nlte:
            self.ts_input['PURE-LTE'] = '.false.'
            self.ts_input['NLTE'] = '.true.'

        if self.atmos_format.lower() != 'marcs':
            self.ts_input['MARCS-FILE'] = '.false.'
            print("NOTE: specific commenting in the model atmosphere is required")


        """ At what wavelenght range to compute a spectrum? """
        if  self.lam_end < self.lam_start:
            tmp = self.lam_start
            self.lam_start = self.lam_end
            self.lam_end = tmp


        """ Compute wavelenght step from resolution """
        self.wave_step = np.mean([self.lam_start, self.lam_end]) / self.resolution

        self.ts_input['LAMBDA_MIN'] = self.lam_start
        self.ts_input['LAMBDA_MAX'] = self.lam_end
        self.ts_input['LAMBDA_STEP'] = self.wave_step

        """ Linelists """
        if type(self.linelist) == np.ndarray:
            pass
        elif type(self.linelist) == str:
            self.linelist = np.array([self.linelist])
        else:
            print("Do not understand linelist argument. Stopped.")
            exit(1)
        print("Linelist(s) will be read from:", self.linelist)

        self.ts_input['NFILES'] = len(self.linelist)
        self.ts_input['LINELIST'] = '\n'.join(self.linelist)

        """
        Find Z (atomic number of the element)
        and check the agreenment with the config file
        """
        if os.path.isfile('./atomic_numbers.dat'):
            el_z = np.loadtxt('./atomic_numbers.dat', usecols=(0))
            el_id = np.loadtxt('./atomic_numbers.dat', usecols=(1), dtype=str)
        else:
            print("Can not find './atomic_numbers.dat' file. Stopped.")
            exit(1)
        for i in range(len(el_id)):
            if self.element.lower() == el_id[i].lower():
                z = el_z[i]
        if not 'element_z' in self.__dict__.keys():
            self.element_z = z
        elif z != self.element_z :
            print(F"config file states Z={self.element_z} for element {self.element}, but according to ./atomic_numbers.dat Z={z}. Stopped.")
            exit(1)


        """ Fitting observations? """
        fitting = False
        if 'fitting' in self.__dict__.keys():
            if self.fitting == 1:
                fitting = True
        self.fitting = fitting

        """ Saving all spectra computed while fitting? """
        save_all_spec = False
        if 'save_all_spec' in self.__dict__.keys():
            if self.save_all_spec == 1:
                save_all_spec = True
        self.save_all_spec = save_all_spec


        if self.fitting:
            if 'observed_path' not in self.__dict__.keys():
                print(F"Specify path to observations. Stopped")
                exit(1)

            print(F"Reading observed spectrum from {self.observed_path}..")
            w_obs, f_obs = read_observations(self.observed_path, self.observed_format)
            self.observed_spectrum = spectrum(w_obs, f_obs, res=self.observed_resolution)

            """ Read linemasks """
            lms = glob.glob( self.line_data + '/*-lmask.txt')
            self.line_masks = {}
            for lm in lms:
                element = lm.split('/')[-1].split('-lmask.txt')[0].lower()
                if self.debug:
                    print(F"Found linemask for {element.capitalize()}")
                w_c, w_start, w_end = np.loadtxt(lm, unpack=True, usecols=(0,1,2), comments=';', ndmin=1)

                if type(w_c) == np.float64:
                    w_c = np.array([w_c])
                    w_start = np.array([w_start])
                    w_end = np.array([w_end])
                w_c, w_start, w_end = zip(*sorted(zip(w_c, w_start, w_end)))
                w_c, w_start, w_end = np.array(w_c), np.array(w_start), np.array(w_end)
                self.line_masks.update({ element : { 'w_center':w_c, 'w_start':w_start, 'w_end':w_end } })

            """
            Segments (windows) are used by default
            Windows are defined by the linemask +- line_offset value (in AA)
            """
            if  'line_offset' not in self.__dict__.keys():
                print("WARNING: seems that you did not set line_offset in the configuration file")
                print("It defines how far to compute the spectrum from the center of the lines in the linemask")
                self.line_offset = 5
                print(F"set to {self.line_offset} AA")


            """ Fitting for macro turbulence? """
            fit_vmac = False
            self.vmac = [np.nan]
            if 'fit_vmac' in self.__dict__.keys():
                if self.fit_vmac == 1:
                    fit_vmac = True
                    self.vmac = np.arange(self.start_vmac, self.end_vmac, self.step_vmac)
                    # [start, end) --> [start, end]
                    if not self.end_vmac in self.vmac:
                        self.vmac = np.hstack((self.vmac, self.end_vmac))
            self.fit_vmac = fit_vmac

            """ Fitting for rotational broadening? """
            fit_vrot = False
            self.vrot = [np.nan]
            if 'fit_vrot' in self.__dict__.keys():
                if self.fit_vrot == 1:
                    fit_vrot = True
                    self.vrot = np.arange(self.start_vrot, self.end_vrot, self.step_vrot)
                    # [start, end) --> [start, end]
                    if not self.end_vrot in self.vrot:
                        self.vrot = np.hstack((self.vrot, self.end_vrot))
            self.fit_vrot = fit_vrot

            if self.debug:
                print(F"Requested Vmac:", self.vmac)
            if self.debug:
                print(F"Requested Vrot:", self.vrot)

            if (self.fit_vmac and not self.line_by_line) or (self.fit_vrot and not self.line_by_line):
                print("fitting for Vmac or Vrot not line by line is not supported yet. Stopped")
                exit(1)
