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


utils.py
"""

import numpy as np
import planck


def relative_humidity_from_dewpoint(temperature, dewpoint):
    """
    relative_humidity_from_dewpoint(temperature, dewpoint)
    see: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_dewpoint.html
    """
    e = saturation_vapor_pressure(dewpoint)
    e_s = saturation_vapor_pressure(temperature)
    return e / e_s


def saturation_vapor_pressure(temperature):
    """
    saturation_vapor_pressure(temperature [K]) [Pa]
    see: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.saturation_vapor_pressure.html#metpy.calc.saturation_vapor_pressure
    """
    return 611.2 * np.exp(17.67 * (temperature - 273.15) / (temperature - 29.65))


def __normalise_lon(longitude):
    """Return longitude(s) normalized from [0, 360) to [-180, 180)."""
    longitude_array = np.asarray(longitude)
    normalised = np.where(longitude_array >= 180.0, longitude_array - 360.0, longitude_array)
    if np.isscalar(longitude):
        return normalised.item()
    return normalised


def __calc_emissivity(Fsd, Fsn, Ts):
    epsilon = Fsn / (Fsd - planck.stefan_boltzmann(Ts))
    return epsilon


def __calc_eq_T(F):
    return planck.inv_stefan_boltzmann(F)


def nanrmse(X, Y):
    return np.sqrt(np.nanmean((X - Y) ** 2))
