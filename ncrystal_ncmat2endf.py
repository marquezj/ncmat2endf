#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  Copyright 2024 ESS Spallation Physics Group                               ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

"""Script for creating a set of ENDF-6 thermal scattering files from a .ncmat
file. Parameters for the ENDF-6 file can be defined with command line arguments or
changing the endf_defaults dictionary in the script. 

The script allows to handle multiple temperatures in one ENDF-6 file, but this is not
recommended, because NCrystal computes an optimal (alpha, beta) grid for each material
and temperature.

Ths script uses the endf-parserpy package from IAEA to format and check the syntaxis of
the ENDF-6 file.
"""

import numpy as np
import scipy.interpolate as scint
import os
import sys
from datetime import datetime
import argparse
import json
import warnings

try:
    import NCrystal as NC
except ImportError:
    raise SystemExit('Could not import NCrystal. Check the package was correctly installed, with a version equal or higher than 3.9.7.\nhttps://github.com/mctools/ncrystal/wiki')
assert NC.version_num >=  3009007, "Too old NCrystal found. Version 3.9.7 or above required"

try:
    import endf_parserpy
    from endf_parserpy.interpreter.fortran_utils import read_fort_floats
    from endf_parserpy.interpreter.fortran_utils import write_fort_floats
except ImportError:
    raise SystemExit('Could not import endf_parserpy. Check the package was correctly installed, with a version equal or higher than 0.11.0.\nhttps://endf-parserpy.readthedocs.io/')
version_num = sum([int(x)*10**(3*n) for x,n in zip(endf_parserpy.__version__.split('.'), reversed(range(3)))])
assert version_num >=  1100, "Too old endf-parserpy found. Version 0.11.0 or above required"


available_elastic_modes = ('greater', 'scaled', 'mixed')
mass_neutron = 1.04540751e-4 #  eV*ps^2*Angstrom^-2 
hbar = 0.658211951e-3 # eV*ps 
kT0 = 0.0253 # eV
T0 = 293.6 # K
        
def endf_roundoff(x):
    """Limit the precision of a float to what can be represented in an ENDF-6 file
       
    Parameters
    ----------
    x : Iterable of float
        Array to process

    Returns
    -------
    numpy array
        Processed array
        
    """
    return np.array(read_fort_floats(write_fort_floats(x), n=len(x)))
    
def wrap_string(inp, lim=66):
    """Return string with lines wrapped to a given width.
       
    Parameters
    ----------
    inp : str
        String to wrap

    lim : int
        Maximum line width

    Returns
    -------
    s: str
        Output string with lines wrapped
        
    """
    ll = []
    for s in inp.split("\n"):
        if s == "": 
            ll.append('')
            continue
        w=0 
        l = []
        if len(s.split()) == 1:
            w = s.split()[0]
            n = len(w) // lim
            for i in range(n):
                ll.append(w[i*lim:(i+1)*lim])
            ll.append(w[n*lim:])
        else:
            for d in s.split():
                if w + len(d) + 1 <= lim:
                    l.append(d)
                    w += len(d) + 1 
                else:
                    ll.append(" ".join(l))
                    l = [d] 
                    w = len(d)
            if (len(l)): ll.append(" ".join(l))
    return('\n'.join(ll))

class ElementData():
    r"""Container for nuclear data for a single element or isotope.

    Attributes
    ----------
    alpha : numpy array
        alpha grid
    beta : numpy array
        positive beta grid
    beta_total : numpy array
        asymmetric beta grid
    sab : iterable of numpy array
        symmetric S(alpha, beta) table
    dwi : iterable of float
        Debye-Waller integral
    teff : iterable of float
        Effective temperatures for short collision time approximation
    elastic : string
        Elastic approximation used in the element (coherent, incoherent or mixed)
    awr : float
        Atomic mass in neutron mass units
    za : int
        ZAID (Z*1000 + A)
    zsymam : float
        Text representation of the element or isotope       
    sigma_i : float
        Incoherent bound atom cross section
    sigma_free : float
        Scattering free atom cross section
    """

    def __init__(self, ad):
        r"""
        Parameters
        ----------
        ad : NCrystal AtomData
        """
        self._sigma_i = ad.incoherentXS()
        self._sigma_free = ad.freeScatteringXS()
        self._awr = ad.averageMassAMU() / NC.constants.const_neutron_mass_amu
        self._dwi = []
        self._alpha = np.array([])
        self._beta = np.array([])
        self._beta_total = np.array([])
        self._sab_total = []
        self._teff = []
        self._elastic = None
        self._element_name = ad.displayLabel()
        Z = '{:3d}'.format(ad.Z())
        A = '   ' if ad.isNaturalElement() else '{:3d}'.format(ad.A())
        sym = self._get_symbol(self._element_name).ljust(2)
        self._zsymam = '{:3s}-{:2s}-{:3s}'.format(Z,sym,A)
        self._za = ad.Z()*1000+ad.A()

    def _get_symbol(self, isotope_name):
        """Get element symbol from isotope name. E.g. Be9 -> Be
    
        Parameters
        ==========
        isotope_name: string

        Returns
        ==========
        element_symbol: string

        """

        symbol = ''        
        for c in isotope_name:
            if (ord(c) >= 97 and ord(c) <= 122) or (ord(c) >= 65 and ord(c) <= 90):
                symbol += c
            else:
                break
        return symbol


    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, x):
        self._alpha = x

    @property
    def beta(self):
        return self._beta
    @beta.setter
    def beta(self, x):
        self._beta = x

    @property
    def beta_total(self):
        return self._beta_total
    @beta_total.setter
    def beta_total(self, x):
        self._beta_total = x

    @property
    def sab_total(self):
        return self._sab_total
    @sab_total.setter
    def sab_total(self, x):
        self._sab_total = x

    @property
    def dwi(self):
        return self._dwi
    @dwi.setter
    def dwi(self, x):
        self._dwi = x

    @property
    def teff(self):
        return self._teff
    @teff.setter
    def teff(self, x):
        self._teff = x

    @property
    def elastic(self):
        return self._elastic
    @elastic.setter
    def elastic(self, x):
        self._elastic = x

    @property
    def awr(self):
        return self._awr

    @property
    def za(self):
        return self._za

    @property
    def zsymam(self):
        return self._zsymam

    @property
    def sigma_free(self):
        return self._sigma_free
    @property
    def sigma_i(self):
        return self._sigma_i
    @sigma_i.setter
    def sigma_i(self, x):
        self._sigma_i = x
        
class NuclearData():
    r"""Container for nuclear data for a material.

    Attributes
    ----------
    ncrystal_comments : string
        Comments in the ncmat file
    temperatures : iterable of flotat
        List of temperatures to process
    ncmat_fn : string
        NCrystal ncmat filename to convert
    composition : iterable of tuples (float, NCrystal AtomData)
        Composition of the material
    elements : iterable of ElementData
        Nuclear data for each of the elements of isotopes in the materal
    edges : iterable of numpy array
        Energies for the Bragg edges for each temperature
    sigmaE : iterable of numpy array
        XS*E for the Bragg edges for each temperature
    """

    def __init__(self, ncmat_fn, temperatures, elastic_mode, vdoslux=3, verbosity=1):
        r"""
        Parameters
        ----------
        ncmat_fn : string
            NCrystal ncmat filename to convert
        temperatures : iterable of float
            List of temperatures to process
        elastic_mode : string
            Elastic approximation used in the material (greater, scaled or mixed)
        vdoslux : integer
            Level of luxury to generate data in NCrystal
        verbosity : integer
            Level of verbosity for the output
        """
        self._temperatures = np.sort(np.asarray(temperatures))
        self._ncmat_fn = ncmat_fn
        mat = NC.load(ncmat_fn+f';vdoslux={vdoslux}')
        self._composition = mat.info.composition
        self._elements = {}
        self._ncrystal_comments = None
        self._vdoslux = vdoslux
        self._verbosity = verbosity
        self._combine_temperatures = False # False: use (alpha, beta) grid for lowest temperature; True: combine all temperatures
        if elastic_mode not in available_elastic_modes:
            raise ValueError(f"Elastic mode '{elastic_mode}' not in {available_elastic_modes}")
        for frac, ad in self._composition:
            element_name = ad.displayLabel()
            self._elements[element_name] = ElementData(ad)
        if mat.info.hasAtomPositions():
            self._edges = []
            self._sigmaE = []
        else:
            self._edges = None
            self._sigmaE = None
        self._incoherent_fraction = -1
        if (len(self._composition) > 1) and (elastic_mode == 'scaled'):
            #
            # Find element with minimum incoherent contribution. 
            #
            self._designated_coherent_atom = None
            for frac, ad in self._composition:
                element_name = ad.displayLabel()
                if (frac/(1.0-frac)*ad.incoherentXS() < self._incoherent_fraction or self._incoherent_fraction == -1):
                    self._incoherent_fraction = frac/(1.0-frac)*ad.incoherentXS()
                    self._designated_coherent_atom = element_name
            if self._verbosity > 1:
                print(f'Designated incoherent: {element_name}')
        self._get_alpha_beta_grid()
        self._get_elastic_data(elastic_mode)
        self._get_inelastic_data()
        self._get_ncrystal_comments()

    @property
    def ncrystal_comments(self):
        return self._ncrystal_comments
        
    @property
    def temperatures(self):
        return self._temperatures

    @property
    def ncmat_fn(self):
        return self._ncmat_fn

    @property
    def composition(self):
        return self._composition
    @composition.setter
    def composition(self, x):
        self._composition = x

    @property
    def elements(self):
        return self._elements
    @elements.setter
    def elements(self, x):
        self._elements = x

    @property
    def edges(self):
        return self._edges
    @edges.setter
    def edges(self, x):
        self._edges = x

    @property
    def sigmaE(self):
        return self._sigmaE
    @sigmaE.setter
    def sigmaE(self, x):
        self._sigmaE = x

    def _get_alpha_beta_grid(self):
        T = self._temperatures[0]
        m = NC.load(f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux}')
        kT = kT0*T/T0 # eV
        for di in m.info.dyninfos:
            element_name = di.atomData.displayLabel()
            if type(di) in [NC.core.Info.DI_VDOS, NC.core.Info.DI_VDOSDebye]:
                sctknl = di.loadKernel(vdoslux=self._vdoslux)
                self._elements[element_name].alpha = sctknl['alpha']*kT/kT0
                self._elements[element_name].beta_total = sctknl['beta']*kT/kT0
            else:
                raise NotImplementedError('Conversion supported only for VDOS and VDOSDebye dyninfos')
        if self._combine_temperatures:
            #
            # Combine (alpha, beta) grids from different temperatures into a single grid.
            # This usually results in a huge grid and it is only kept as an option
            # to debug libraries.
            #
            for T in self._temperatures[1:]:
                m = NC.load(f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux}')
                kT = kT0*T/T0 # eV
                for di in m.info.dyninfos:
                    element_name = di.atomData.displayLabel()
                    if type(di) in [NC.core.Info.DI_VDOS, NC.core.Info.DI_VDOSDebye]:
                        sctknl = di.loadKernel(vdoslux=self._vdoslux)
                        self._elements[element_name].alpha = np.unique(np.concatenate((self._elements[element_name].alpha, sctknl['alpha']*kT/kT0)))
                        self._elements[element_name].beta_total = np.unique(np.concatenate((self._elements[element_name].beta_total, sctknl['beta']*kT/kT0)))
        for frac, ad in self._composition:
            element_name = ad.displayLabel()
            #
            # Remove points that cannot be represented in ENDF data
            #
            self._elements[element_name].alpha = np.unique(endf_roundoff(self._elements[element_name].alpha))
            self._elements[element_name].beta_total = np.unique(endf_roundoff(self._elements[element_name].beta_total))
            x = self._elements[element_name].beta_total[np.where(self._elements[element_name].beta_total<=0)] # get negative beta
            self._elements[element_name].beta = -x[::-1] # Invert beta and change sign
            self._elements[element_name].beta[0] = 0.0
            if self._verbosity > 2:
                print(f'>>> alpha points: {len(self._elements[element_name].alpha)}, alpha range: ({np.min(self._elements[element_name].alpha*kT0/kT)}, {np.max(self._elements[element_name].alpha*kT0/kT)})')
                print(f'>>> beta points: {len(self._elements[element_name].beta)}, beta range: ({np.min(self._elements[element_name].beta*kT0/kT)}, {np.max(self._elements[element_name].beta*kT0/kT)})')

    def _get_elastic_data(self, elastic_mode):
        for T in self._temperatures:
            m = NC.load(f'{self._ncmat_fn};temp={T}K;vdoslux={self._vdoslux};comp=bragg;dcutoff=0.1')
            if m.info.hasAtomPositions():
                #
                # Load coherent elastic data
                #
                if T == self._temperatures[0]:
                    # Find unique Bragg edges, as represented in ENDF-6 floats
                    edges = np.array([NC.wl2ekin(2.0*e.dspacing) for e in m.info.hklObjects()])
                    self._edges = np.unique(endf_roundoff(edges))
                    # Coherent scattering XS is evaluated between edges
                    eps = 1e-3
                    self._evalpoints = np.concatenate((self._edges[:-1]**(1-eps)*self._edges[1:]**eps, [self._edges[-1]]))
                sigmaE = endf_roundoff(m.scatter.xsect(self._evalpoints)*self._evalpoints)
                assert np.all(sigmaE[:-1] <= sigmaE[1:]), 'Sigma*E in Bragg edges not cummulative'
                self._sigmaE.append(sigmaE)
            #
            for di in m.info.dyninfos:
                element_name = di.atomData.displayLabel()
                if type(di) in [NC.core.Info.DI_VDOS, NC.core.Info.DI_VDOSDebye]:
                    emin = di.vdosData()[0][0]
                    emax = di.vdosData()[0][1]
                    rho = di.vdosData()[1]
                    res = NC.analyseVDOS(emin, emax, rho, di.temperature, di.atomData.averageMassAMU())
                    #
                    # Load incoherent elastic data
                    #
                    if m.info.stateOfMatter().name == 'Solid':
                        msd = res['msd']
                        self._elements[element_name].dwi.append(msd*2*mass_neutron/hbar**2)
        if self._verbosity > 1:
            print('>> Prepare elastic approximations')
        if elastic_mode == 'scaled' and self._incoherent_fraction < 1e-6:
            elastic_mode = 'greater'
            if self._verbosity>1:
                print(f'>> Scaled elastic mode requested but all elements are coherent.')
        for frac, ad in self._composition:
            element_name = ad.displayLabel()
            if elastic_mode == 'mixed': # iel = 100
                if (self._sigmaE == None):      # mixed elastic requested but only incoherent available
                    if self._verbosity>1:
                        print(f'>> Mixed elastic mode for {element_name} but no Bragg edges found: incoherent approximation')
                    self._elements[element_name].sigma_i = (ad.incoherentXS() + ad.coherentXS())
                    self._elements[element_name].elastic = 'incoherent'
                else:
                    if self._verbosity>1:
                        print(f'>> Mixed elastic mode for {element_name}')
                    self._elements[element_name].elastic = 'mixed'
            if elastic_mode == 'greater': # iel = 98
                if (ad.incoherentXS() > ad.coherentXS())  or (self._sigmaE == None):
                    if self._verbosity>1:
                        print(f'>> Principal elastic mode for {element_name}: incoherent')
                    self._edges = None
                    self._sigmaE = None
                    self._elements[element_name].elastic = 'incoherent'
                else:
                    if self._verbosity>1:
                        print(f'>> Principal elastic mode for {element_name}: coherent')
                    self._elements[element_name].sigma_i =  None
                    self._elements[element_name].dwi =  None
                    self._elements[element_name].elastic = 'coherent'
            elif elastic_mode == 'scaled':
                if len(self._composition) == 1: # iel = 99, single atomic case
                    if (ad.incoherentXS() > ad.coherentXS()) or (self._sigmaE == None):
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for single atom {element_name}: incoherent')
                        self._edges = None
                        self._sigmaE = None
                        self._elements[element_name].sigma_i = (ad.incoherentXS() + ad.coherentXS())
                        self._elements[element_name].elastic = 'incoherent'
                    else:
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for single atom {element_name}: coherent')
                        self._elements[element_name].sigma_i =  None
                        self._elements[element_name].dwi =  None
                        self._elements[element_name].elastic = 'coherent'
                        self._sigmaE = [x*(ad.incoherentXS() + ad.coherentXS())/ad.coherentXS() for x in self._sigmaE]
                elif (self._sigmaE == None):      # iel = 99, incoherent approximation
                    if self._verbosity>1:
                        print(f'>> Scaled elastic mode for {element_name}: incoherent approximation')
                    self._elements[element_name].sigma_i = (ad.incoherentXS() + ad.coherentXS())
                    self._elements[element_name].elastic = 'incoherent'
                else:                              # iel = 99, multi atomic case
                    if element_name == self._designated_coherent_atom:
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for {element_name} in compound: designated coherent atom, dividing by frac^2={frac**2}')
                        self._elements[element_name].sigma_i =  None
                        self._elements[element_name].dwi =  None
                        self._elements[element_name].elastic = 'coherent'
                        self._sigmaE = [x/frac for x in self._sigmaE]
                    else:
                        if self._verbosity>1:
                            print(f'>> Scaled elastic mode for {element_name} in compound: incoherent')
                        self._elements[element_name].elastic = 'incoherent'
                        self._elements[element_name].sigma_i = (1.0+self._incoherent_fraction/ad.incoherentXS())*self._elements[element_name].sigma_i
        
    def _get_inelastic_data(self):
        for T in self._temperatures:
            m = NC.load(f'{self._ncmat_fn};temp={T}K')
            for di in m.info.dyninfos:
                element_name = di.atomData.displayLabel()
                if type(di) in [NC.core.Info.DI_VDOS, NC.core.Info.DI_VDOSDebye]:
                    #
                    # Load incoherent inelastic data
                    #
                    sctknl = di.loadKernel(vdoslux=self._vdoslux)
                    if self._verbosity > 2:
                        print(f'>>> Interpolating T={T}K for {element_name}')
                
                    alpha = sctknl['alpha']
                    beta = sctknl['beta']
                    sab = sctknl['sab']
                    sab.shape = (len(beta), len(alpha))
                    kT = kT0/T0*T # eV
                    alpha_grid, beta_grid = np.meshgrid(self._elements[element_name]._alpha*kT0/kT, self._elements[element_name]._beta_total*kT0/kT)
                    points0 = np.column_stack((alpha_grid.ravel(), beta_grid.ravel()))
                    #
                    # We need to interpolate S(a,b) because the NCrystal grid might contain numbers that cannot be represented
                    # as distinct FORTRAN reals in the ENDF-6 file
                    #
                    sab_int = scint.interpn((alpha, beta), sab.transpose(), points0, bounds_error=False, fill_value=0.0, method='linear')
                    sab_int.shape = np.shape(alpha_grid)
                    self._elements[element_name]._sab_total.append(sab_int) 
                    emin = di.vdosData()[0][0]
                    emax = di.vdosData()[0][1]
                    rho = di.vdosData()[1]
                    res = NC.analyseVDOS(emin, emax, rho, di.temperature, di.atomData.averageMassAMU())
                    self._elements[element_name]._teff.append(res['teff'])
                else:
                    self._elements[element_name]._sab = None
                    self._elements[element_name]._teff = None
    def _get_ncrystal_comments(self):
        ll = []
        for l in NC.createTextData(self._ncmat_fn).rawData.split('\n')[:]:
            if len(l) > 0 and l[0] == '#':
                ll.append(l[1:])
        self._ncrystal_comments = wrap_string("\n".join(ll),66)

class EndfFile():
    r"""Creates thermal ENDF file. 
       using endf-parserpy

    Methods
    ----------
    write(endf_fn)
        Write ENDF file.
    """
    def __init__(self, element_name, data, mat, endf_parameters, include_gif=False, isotopic_expansion=False, verbosity=1):
        r"""
        Parameters
        ----------
        element_name : string
            Element to be output 
    
        data : NuclearData
            Nuclear data for the material
    
        mat: int
            ENDF material number
            
        endf_parameters : EndfParameters
            Parameters for the ENDF-6 file
        
        include_gif: boolean
            Include the generalized information in MF=7/MT=451 in isotopes
        
        isotopic_expansion: boolean
            Expand the information in MF=7/MT=451 in isotopes
        
        verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        self._endf_dict = endf_parserpy.EndfDict()
        self._parser = endf_parserpy.EndfParser(explain_missing_variable=True, cache_dir=False)
        self._element_name = element_name
        self._mat = mat
        assert ((not isotopic_expansion) or include_gif), "Isotopic expansion requires generalized information file, use --gif"
        self._include_gif = include_gif
        self._isotopic_expansion = isotopic_expansion
        self._verbosity = verbosity
        self._endf_dict['0/0'] = {}
        self._endf_dict['0/0']['MAT'] = self._mat
        self._endf_dict['0/0']['TAPEDESCR'] = 'Created with ncmat2endf'
        self._createMF1(data, endf_parameters)
        self._createMF7(data, endf_parameters)
        endf_parserpy.update_directory(self._endf_dict, self._parser)

    def _createMF7(self, data, endf_parameters):
        """Creates MF=7 file of a thermal ENDF file. 
           See ENDF-102, sect. 7.
    
        Parameters
        ----------
        data : NuclearData
            Nuclear data for the material
    
        endf_parameters : EndfParameters
            Parameters for the ENDF-6 file

        verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        if self._verbosity > 1:
            print('> Generate MF7')
        awr = data.elements[self._element_name].awr
        mat = self._mat
        za = data.elements[self._element_name].za
        temperatures = data.temperatures
        #
        # Prepare dictionary for thermal elastic reaction
        #
        elastic = data.elements[self._element_name].elastic            
        if elastic is not None:
            self._endf_dict['7/2'] = {}
            d = self._endf_dict['7/2']
            d['MAT'] = mat
            d['ZA'] = za
            d['AWR'] = awr
            if elastic in ['coherent', 'mixed']:
                edges  = data.edges
                sigmaE = data.sigmaE[0]
                d['T0'] = temperatures[0]
                d['LT'] = len(temperatures)-1
                d['S_T0_table'] = {}
                d['S_T0_table']['NBT'] = [len(edges)]
                d['S_T0_table']['INT'] = [1]
                d['S_T0_table']['Eint'] = edges.tolist()
                d['S_T0_table']['S'] = sigmaE.tolist()
                d['T'] = {k: v for k, v in enumerate(temperatures[1:], start=1)}
                S = {}
                for q, v in enumerate(edges, start=1):
                    S[q] = {}
                    for i,v in enumerate(temperatures[1:], start=1):
                        S[q][i] = data.sigmaE[i][q-1]
                d['S'] = S
                d['LI'] = 2
                d['NP'] = len(edges)
            if elastic in ['incoherent', 'mixed']:
                d['SB'] = data.elements[self._element_name].sigma_i
                d['Wp']  = data.elements[self._element_name].dwi
                d['Tint'] = temperatures.tolist()
                d['INT'] = [2]
                d['NBT'] = [len(temperatures)]
            lthr_values = {'coherent':1, 'incoherent':2, 'mixed':3}
            d['LTHR'] = lthr_values[data.elements[self._element_name].elastic]
        #
        # Prepare dictionary for thermal inelastic reaction
        #
        self._endf_dict['7/4'] = {}
        d = self._endf_dict['7/4']
        d['MAT'] = mat
        d['ZA'] = za
        d['AWR'] = awr
        d['LAT'] =  1 # (alpha, beta) grid written for kT0 = 0.0253 eV
        d['LASYM'] = endf_parameters.lasym # symmetric/asymmetric S(a,b)
        d['LLN'] = 0 # linear S is stored
        d['NI'] = 6
        d['NS'] = 0
        d['B'] = {1:data.elements[self._element_name].sigma_free, 
                  2:endf_parameters.emax, 
                  3:awr, 
                  4:endf_parameters.emax, 
                  5:0,                                # unused
                  6:1                                 # natom
                 }
        alpha = data.elements[self._element_name].alpha
        beta = data.elements[self._element_name].beta
        sab_data = []
        for sab_total, T in zip(data.elements[self._element_name].sab_total, temperatures):
            kT = T/T0*kT0            
            alpha_grid, beta_grid = np.meshgrid(alpha*kT0/kT, data.elements[self._element_name].beta_total*kT0/kT)
            if (endf_parameters.lasym == 0) or (endf_parameters.lasym == 1):
                sym = np.exp(beta_grid/2)
            if endf_parameters.lasym == 3:
                # S(a,b) for all beta
                sab_data.append(sab_total.transpose())
                continue
            if endf_parameters.lasym == 2:
                # S(a,b) for negative beta
                sab_sym2 = sab_total[np.where(beta_grid<=0)] # get negative branch of S(a,b)
                sab_sym2.shape = (len(beta), len(alpha))
                sab_sym3 = sab_sym2[::-1,:]  # Invert S(a,b) for negative beta
                sab_data.append(sab_sym3.transpose())
                continue
            if endf_parameters.lasym == 1:
                # S(a,b)*exp(-b/2) for all beta
                sab_sym = sab_total*sym
                sab_data.append(sab_sym.transpose())
                continue
            if endf_parameters.lasym == 0:
                # S(a,b)*exp(-b/2) for negative beta
                sab_sym = sab_total*sym
                sab_sym2 = sab_sym[np.where(beta_grid<=0)] # get negative branch of S(a,b)
                sab_sym2.shape = (len(beta), len(alpha))
                sab_sym3 = sab_sym2[::-1,:]  # Invert S(a,b) for negative beta
                sab_data.append(sab_sym3.transpose())
                continue
                    
        if (endf_parameters.lasym == 0) or (endf_parameters.lasym == 2):
            # Save S(a,b) or S(a,b)*exp(-b/2) for negative beta
            d['NB'] = len(beta) 
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]
        
            d['T0'] = temperatures[0]
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1 for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < endf_parameters.smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = (alpha/awr).tolist()
                S1[j]['S'] = sab[:,j-1].tolist()
            d['S_table'] = S1
                
            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < endf_parameters.smin:
                                    sval = 0.0
                                sab.append(sval)
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2
        elif (endf_parameters.lasym == 1) or (endf_parameters.lasym == 3):
            # Save S(a,b) or S(a,b)*exp(-b/2) for all beta
            alpha = data.elements[self._element_name].alpha
            beta = data.elements[self._element_name].beta_total
            
            d['NB'] = len(beta) 
            d['beta_interp/NBT'] = [len(beta)]
            d['beta_interp/INT'] = [4]
        
            d['T0'] = temperatures[0]
            d['beta'] = {k:v for k, v in enumerate(beta,start=1)}
            d['LT'] = {k:len(temperatures)-1 for k, v in enumerate(beta,start=1)}
            d['T'] = {k:v for k,v in enumerate(temperatures[1:], start=1)}
            d['LI'] = {k:4 for k,v in enumerate(temperatures[1:], start=1)}
            d['NP'] = len(alpha)
            S1 = {}
            sab = sab_data[0]
            sab[sab < endf_parameters.smin] = 0.0
            for j,v in enumerate(beta, start=1):
                S1[j] = {}
                S1[j]['NBT'] = [len(alpha)]
                S1[j]['INT'] = [4]
                S1[j]['alpha'] = (alpha/awr).tolist()
                S1[j]['S'] = sab[:,j-1].tolist()
            d['S_table'] = S1
            
            S2 = {}
            if len(temperatures) > 1:
                for q, v in enumerate(alpha, start=1):
                    S2[q] = {}
                    for j,v in enumerate(beta, start=1):
                        sab = []
                        for i,v in enumerate(temperatures[1:], start=1):
                                sval = sab_data[i][q-1,j-1]
                                if sval < endf_parameters.smin:
                                    sval = 0.0
                                sab.append(sval)
                        S2[q][j] = {k: v for k,v in enumerate(sab, start =1)}
            d['S'] = S2

        d['teff0_table/Teff0'] = data.elements[self._element_name].teff
        d['teff0_table/Tint'] = temperatures.tolist()
        d['teff0_table/NBT'] = [len(temperatures)]
        d['teff0_table/INT'] = [2]
        if self._include_gif:        
            if self._isotopic_expansion:
                raise NotImplementedError('Isotopic expansion not yet implemented')
            else:
                self._endf_dict['7/451'] = {}
                d = self._endf_dict['7/451']
                d['MAT'] = mat
                d['ZA'] = za
                d['AWR'] = awr
                d['NA'] = 1
                d['NAS'] = 1
                d['NI'] = {1:1}
                d['ZAI'] = {1:{1:za}}
                d['LISI'] = {1:{1:0}}
                d['AFI'] = {1:{1:1.0}}
                d['SFI'] = {1:{1:data.elements[self._element_name].sigma_free}}
                d['AWRI'] = {1:{1:awr}}

    def _createMF1(self, data, endf_parameters):
        """Creates MF=1 file of a thermal ENDF file. 
           See ENDF-102, sect. 1.
    
        Parameters
        ----------
        endf_dict:
            endf-parserpy dictionary
    
        element_name : str
            Element to write
    
        data : dictionary
            Dictionary containing the nuclear data extrated by get_nuclear_data()
    
        self._verbosity : int
            Level of verbosity of the output (0: quiet)
        """
        awr = data.elements[self._element_name].awr
        mat = self._mat
        za = data.elements[self._element_name].za
        zsymam = data.elements[self._element_name].zsymam
        self._endf_dict['1/451'] = {}
        d = self._endf_dict['1/451']
        d['MAT'] = mat;     d['ZA'] = za;     d['AWR'] = awr;     d['LRP'] = -1
        d['LFI'] = 0;       d['NLIB'] = endf_parameters.nlib;    d['NMOD'] = 0;      d['ELIS'] = 0
        d['LIS'] = 0;       d['LISO'] = 0;    d['STA'] = 0;       d['NFOR'] = 6
        d['AWI'] = 1.0;     d['EMAX'] = endf_parameters.emax;  d['LREL'] = endf_parameters.lrel
        d['NSUB'] = 12;     d['NVER'] = endf_parameters.nver;           d['TEMP'] = 0.0
        d['LDRV'] = 0
        d['HSUB/1'] = f'----{endf_parameters.libname:18s}MATERIAL {mat:4d}'.ljust(66)
        d['HSUB/2'] = f'-----THERMAL NEUTRON SCATTERING DATA'.ljust(66)
        d['HSUB/3'] = f'------ENDF-6 FORMAT'.ljust(66)
        d['NXC'] = 1 
        d['ZSYMAM'] = zsymam.ljust(11)
        d['ALAB'] = endf_parameters.alab.ljust(11)
        d['AUTH'] = endf_parameters.auth.ljust(33)
        d['REF'] = endf_parameters.reference.ljust(21)
        now = datetime.now()
        months = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AGO', 'SEP', 'OCT', 'NOV', 'DEC')
        edate= f'EVAL-{months[now.month-1]}{now.strftime("%y")}'
        ddate= f'DIST-{months[now.month-1]}{now.strftime("%y")}'
        rdate= f'REV{endf_parameters.lrel:1d}-{months[now.month-1]}{now.strftime("%y")}'
        d['EDATE'] = edate
        d['DDATE'] = ddate
        d['RDATE'] = rdate
        d['ENDATE'] = endf_parameters.endate.ljust(8)
        description = []
        description.append(''.ljust(66))
        description.append(' Converted from:'.ljust(66))
        description.append(data.ncmat_fn.center(66))
        description.append(' Temperatures:'.ljust(66))
        for T in data.temperatures:
            description.append(f'       {T:.2f} K'.ljust(66))
        description.append(''.ljust(66))
        description.append(66*'*')
        description.append(''.ljust(66))
        description.append('Comments from NCMAT file:'.ljust(66))
        description.append(''.ljust(66))
        for line in data.ncrystal_comments.split('\n'):
            description.append(line.ljust(66))
        # description.append(''.ljust(66))
        description.append(66*'*')
        d['DESCRIPTION'] = {k:v for k, v in enumerate(description, start=1)}
        d['NWD'] = 5+len(description)
        d['MFx/1'] = 1;   d['MTx/1'] = 451;  d['NCx/1'] = 5;  d['MOD/1'] = d['NMOD']

    def write(self, endf_fn):
        if self._verbosity > 0:
            print(f'Write ENDF file {endf_fn}...')
        with open(endf_fn, mode='w') as f:
            f.write('\n'.join(self._parser.write(self._endf_dict, zero_as_blank=True)))

class EndfParameters():
    """Parameters for the ENDF-6 file

    Attributes
    ----------
    alab : string
        Mnemonic for the originating laboratory 

    smin : float
        Minimum value of S(alpha, beta) to be stored in the file

    libname : string
        Name of the library

    nlib : int
        Library identifier (e.g. NLIB= 0 for ENDF/B).

    auth : string
        Author(s) name(s).

    reference : string
        Primary reference for the evaluation.

    emax : float
        Upper limit of the energy range for evaluation (eV).

    lrel : int
        Library release number.

    nver : int
        Library version number.

    endate: string
        Master File entry date in the form YYYYMMDD.
    """

    def __init__(self):
        self._alab = 'MyLAB'
        self._auth = 'NCrystal'
        self._reference = 'REFERENCE'
        self._nver = 1
        self._libname = 'MyLib'
        self._endate = 'YYYYMMDD'
        self._nlib = 0
        self._lrel = 0
        self._smin = 1e-100
        self._emax = 5.0
        self._lasym = 0 # Symmetric S(a,b) as default 

    @property
    def alab(self):
        return self._alab
    @alab.setter
    def alab(self, x):
        self._alab = x

    @property
    def smin(self):
        return self._smin
    @smin.setter
    def smin(self, x):
        self._smin = x

    @property
    def libname(self):
        return self._libname
    @libname.setter
    def libname(self, x):
        self._libname = x

    @property
    def nlib(self):
        return self._nlib
    @nlib.setter
    def nlib(self, x):
        self._nlib = x

    @property
    def auth(self):
        return self._auth
    @auth.setter
    def auth(self, x):
        self._auth = x

    @property
    def reference(self):
        return self._reference

    @property
    def emax(self):
        return self._emax

    @property
    def lrel(self):
        return self._lrel

    @property
    def nver(self):
        return self._nver
        
    @property
    def endate(self):
        return self._endate

    @property
    def lasym(self):
        return self._lasym
    @lasym.setter
    def lasym(self, x):
        self._lasym = x

def ncmat2endf(ncmat_fn, name, endf_parameters, temperatures=(293.6,), mat_numbers=None, elastic_mode='scaled', include_gif=False, isotopic_expansion=False, vdoslux=3, verbosity=1):
    """Generates a set of ENDF-6 formatted files for a given NCMAT file.
       
    Parameters
    ----------
    ncmat_fn : str
        Filename of the ncmat file to convert

    temperatures : float or iterable of float
        Temperatures in Kelvin to generate the nuclear data

    mat_numbers : dict of str to int
        Material number for each element

    elastic_mode : str
        Treatment mode for the elastic component
        "greater" = only the greater ellastic component (coherent or incoherent) is saved
        "mixed"   = both the coherent and incoherent inelastic components are saved
        "scaled"  = for monoatomic scatterers, the major component is saved, scaled to the total bound XS
                    for polyatomic scatterers, coherent scattering for the whole system is assigned to the
                    atom with minimum incoherent cross section, and its incoherent contribution is distributed
                    among the other atoms

    include_gif: boolean
        Include the generalized information in MF=7/MT=451 in isotopes

    isotopic_expansion: boolean
        Expand the information in MF=7/MT=451 in isotopes

    vdoslux : integer
        Level of luxury to generate data in NCrystal

    verbosity : int
        Level of verbosity of the output (0: quiet)

    Returns
    -------
    file_names: list of (str, float)
        List of tuples contanining the ENDF-6 files and their fraction in the composition
        
    """
    if verbosity > 0 and endf_parameters.lasym > 0:
        print(f'Creating non standard S(a,b) with LASYM = {endf_parameters.lasym}')
        
    if type(temperatures) in [int, float]:
        temperatures = (temperatures,)
    
    if len(temperatures) > 1:
        warnings.warn('Multiple temperatures requested. Although this is supported, '
        +'it is not recommended because NCrystal generates a custom (alpha,beta) grid for each temperature. '
        +'The (alpha,beta) grid for first temperature will be used, and S(alpha, beta) for other temperatures will be interpolated.')
    if verbosity > 0:
        print('Get nuclear data...')
        
    data = NuclearData(ncmat_fn, temperatures, elastic_mode, vdoslux, verbosity)

    if mat_numbers is not None:
        n = len(mat_numbers)
        for frac, ad in data.composition:
            if ad.displayLabel() in mat_numbers.keys():
                n = n - 1
        assert n==0, 'Incorrect material number assignement'

    file_names = []
    for frac, ad in data.composition:
        element_name = ad.displayLabel()
        mat = 999 if mat_numbers is None else mat_numbers[element_name]
        endf_fn = f'tsl_{name}.endf' if element_name == name else f'tsl_{element_name}_in_{name}.endf'
        if data.elements[element_name].sab_total is not None:
            endf_file = EndfFile(element_name, data, mat, endf_parameters, include_gif, isotopic_expansion, verbosity)
            endf_file.write(endf_fn)
            file_names.append((endf_fn, frac))
        else:
            if verbosity > 0:
                print(f'Scattering kernel not available for: {endf_fn}')
    return(file_names)

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_args(args=sys.argv[1:]):
    '''Parse arguments.'''
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter)
    endf_defaults = EndfParameters()
    parser.add_argument('input',
                        help='NCMAT file to convert')
    parser.add_argument('-n', '--name',
                        help='Name of the compound in the NCMAT file')
    parser.add_argument('-t', '--temperatures',
                        nargs='+',
                        help='Temperatures to process the NCMAT file', type=float, default=[293.6])
    parser.add_argument('-e', '--elastic_mode',
                        help='Approximation used for the elastic component', type=str,
                        choices=available_elastic_modes,
                        default='scaled')
    parser.add_argument('-v', '--verbosity',
                        help='Controls how verbose should be the output',
                        type=int, default=1)
    parser.add_argument('-l', '--luxury',
                        help='Set the NCrystal vdoslux parameter used to generate the library. 3 is normal, 4 is fine and 5 is very fine.',
                        type=int, default=3, choices=range(1, 6))
    parser.add_argument('-g', '--gif',
                        help='Include the generalized information file (MF=7/MT=451)',
                        action='store_true')
    parser.add_argument('-i', '--isotopic',
                        help='Expand each scatterer element into its isotopes',
                        action='store_true')
    parser.add_argument('-m', '--mats',
                        help='JSON dictionary containing material number assignement for each element, e.g. \'{"C":37, "H": 38}\'',
                        type=json.loads)
    parser.add_argument('--alab',
                        help='Set the ALAB parameter in MF1/MT451', type=str,
                        default=endf_defaults.alab)
    parser.add_argument('--auth',
                        help='Set the AUTH parameter in MF1/MT451', type=str,
                        default=endf_defaults.auth)
    parser.add_argument('--libname',
                        help='Set the LIBNAME parameter in MF1/MT451', type=str,
                        default=endf_defaults.libname)
    parser.add_argument('--nlib',
                        help='Set the NLIB parameter in MF1/MT451', type=str,
                        default=endf_defaults.nlib)
    parser.add_argument('--smin',
                        help='Set the minimum value of S(alpha, beta) stored in MF7/MT4',
                        type=float, default=endf_defaults.smin)
    parser.add_argument('--lasym', help='Write symmetric S(a,b) table', 
                       type=int, default=0, choices=range(0, 4))
    return parser.parse_args(args)

if __name__ == '__main__':

    options = parse_args()
    fn = options.input
    name = options.name
    temperatures = options.temperatures
    elastic_mode = options.elastic_mode
    verbosity = options.verbosity
    mat_numbers = options.mats
    vdoslux = options.luxury
    include_gif = options.gif
    isotopic_expansion = options.isotopic
    params = EndfParameters()
    params.alab = options.alab
    params.auth = options.auth
    params.smin = options.smin
    params.libname = options.libname
    params.nlib = options.nlib
    params.lasym = options.lasym
    
    file_names = ncmat2endf(fn, name, params, temperatures, mat_numbers, elastic_mode, include_gif, isotopic_expansion, vdoslux, verbosity)
    if verbosity > 0:
        print('Files created:')
        for fn, frac in file_names: print(f'  {fn}')
    sys.exit(0)
