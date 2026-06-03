"""
Copyright (C) 2025
Claudio Belotti (CNR-INO, Italy), claudio.belotti@cnr.it
Lorenzo Pattelli (INRiM, Italy), l.pattelli@inrim.it


    This file is part of PASTICHE:
    the PaRaMetriC Atmospheric Spectral Tool for Irradiance Calculation using Hourly ERA5 data.

    PASTICHE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PASTICHE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PASTICHE.  If not, see <http://www.gnu.org/licenses/>.


convert.py
"""

import numpy as np
from scipy.constants import R, N_A

def p2z(p, p0=1013.0):
    """Convert pressure in hPa to altitude in km, assuming sea-level pressure p0."""
    z = -16.0 * np.log10(p / p0)
    return z


def geopotential2height(phi):
    # https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.geopotential_to_height.html
    Re = 6371*1e3 #average Earth radius m
    g0 = 9.80665 # standard gravity ms**-2
    z = phi*Re/(g0*Re-phi)
    return z/1000 #km


def height2geopotential(z):
    #https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.height_to_geopotential.html
    Re = 6371*1e3 #average Earth radius m
    g0 = 9.80665 # standard gravity ms**-2
    z *= 1e3 #m
    phi = g0*Re*z/(Re+z)
    return phi


def q2mmr(q):
    '''
    Calcola il mass mixing ratio a partire da umidità specifica.
    q = mv/m = mv/(mv+md), mv massa vapore, md massa aria dry, m massa aria
    MMR = mv/md = q/(1-q)
    '''
    return q / (1 - q)


def mmr2vmr(MMR, species):
    '''
    Calcola il volume mixing ratio in ppm 
    a partire dal mass mixing ratio.
    
    '''
    species = species.lower()
    molar_mass = {
        'o3': 47.9982,
        'h2o': 18.0153 }
    md = 28.9644
    vmr = MMR * md / molar_mass[species]
    return vmr * 1e+06

def n_density2vmr(C, T, p):
    '''
    Convert number density or number concentration [molecules cm**-3] to VMR.

    C [molecules cm**-3]
    T [K]
    p [hPa]
    '''
    C *= 1e+06 # [molecules m**-3]
    p *= 100 # [Pa]
    
    vmr = C*R*T/(N_A * p)
    return vmr * 1e+06
