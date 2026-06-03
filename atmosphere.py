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


atmosphere.py
"""

import numpy as np
import convert
import era5
from ig2 import IG2
import logger
import rrtm.input


def clwc2cwp(clwc, q, p, T, DZ):
    ## FIXME usare __intermediate dove appropriato
    """
    clwc: mass mixing ratio [kg/kg]
    q: specific humidity [kg/kg]
    p: pressure [hPa]
    T: temperature [K]
    DZ: cloud thickness [m]
    """

    Mw = 18.01528e-3  # molar mass [kg/mol]
    Md = 28.9652e-3  # molecular mass of dry air  [kg/mol]
    eps_air = Mw / Md

    Rgas = 8.314472e-2  # [m3 hPa K-1 mol-1]
    Rdry = Rgas / Md

    Rmoist = Rdry * (1 - ((1 - eps_air) / eps_air) * q)

    # P pressione in hPa, T temp in K, DZ spessore strato in m.
    cwp = clwc * p / (Rmoist * T) * DZ * 1000  # [g/m2]
    return cwp


def clwc2cwp_simple(clwc, p):
    g = 9.80665  # m/s**2
    dp = 100 * np.abs(np.diff(p))  # Pa
    clwc_i = __intermediate(clwc) * 1000  # g/kg
    cwp = clwc_i * dp / g
    return cwp


def __intermediate(x):
    """
    utility function to get the intermediate values in a vector
    """
    return 0.5 * (x[0:-1] + x[1:])


def get_co2(lt, date):
    p = IG2.by_date(date, "pre", "in/ig2").by_lat(lt).data
    co2 = IG2.by_date(date, "co2", "in/ig2").by_lat(lt).data
    y = date.year
    if y > 2012:
        co2 = co2 + 2.73 * (y - 2012)  # see email M.Ridolfi 20240325T1749
    elif y < 2002:
        co2 = co2 - 2.73 * (2002 - y)  # solo per test FIXME

    return np.column_stack((p, co2))


def prepare_atmosphere(
    description,
    ps,
    p,
    zs,
    z,
    Ts,
    T,
    q,
    o3,
    model,
    ig2_co2,
    cf=None,
    clwc=None,
    ciwc=None,
    T2m=None,
    add_one_at_surf=False,
):
    if p[0] < p[1]:
        p = np.flipud(p)
        z = np.flipud(z)
        T = np.flipud(T)
        q = np.flipud(q)
        o3 = np.flipud(o3)
        if cf is not None:
            cf = np.flipud(cf)
        if clwc is not None:
            clwc = np.flipud(clwc)
        if ciwc is not None:
            ciwc = np.flipud(ciwc)

    # converto a vmr per via di bug di RTTMG
    # vedi https://github.com/AER-RC/RRTMG_LW/issues/14
    h2o = convert.mmr2vmr(convert.q2mmr(q), species="h2o")
    o3 = convert.mmr2vmr(o3, species="o3")

    ps /= 100  # -> hPa
    ps = np.round(ps)
    h = era5.z2h(z)
    h /= 1000  # km
    hs = era5.z2h(zs)
    hs /= 1000  # km
    logger.debug(f"ps={ps}, h[:2]={np.round(h[:2] * 1000)}, hs={round(hs * 1000)}")

    if add_one_at_surf:  # add one layer at surface
        # indices of levels above psurf and hsurf
        # hmsk = np.nonzero(
        #     np.logical_and( levels<ps, np.round(h,decimals=3)>np.round(hs,decimals=3)))
        hmsk = np.nonzero(p < ps)
        p = p[hmsk]
        h = h[hmsk]
        T = T[hmsk]
        h2o = h2o[hmsk]
        o3 = o3[hmsk]

        if T2m is None:
            logger.debug(f"no T2m, setting T0 to Ts={Ts:.1f}.")
            T0 = Ts
        else:
            T0 = T2m
            logger.debug("setting T0=%.1f to T2m=%.1f (Ts=%.1f).", T0, T2m, Ts)

        logger.debug(
            "adding one level at surface hs=%.2f m, ps=%.2f, T0=%.2f K, h2o=%.4e, o3=%.4e",
            1000 * hs,
            ps,
            T0,
            h2o[0],
            o3[0],
        )
        p = np.insert(p, 0, ps)
        h = np.insert(h, 0, hs)
        T = np.insert(T, 0, T0)
        h2o = np.insert(h2o, 0, h2o[0])
        o3 = np.insert(o3, 0, o3[0])
    else:
        logger.debug("no level added at surface")

    co2 = np.interp(p, ig2_co2[:, 0], ig2_co2[:, 1])
    gases = [h2o, co2, o3]
    atm = rrtm.input.Atmosphere(
        description,
        h,
        p,
        "A",
        T,
        "A",
        Ts,
        gases,
        "AAA" + str(model) * 4,
        model,
        semiss=0.95,
    )

    cld = None
    if clwc is not None:
        if add_one_at_surf:
            clwc = clwc[hmsk]
            ciwc = ciwc[hmsk]
            cf = cf[hmsk]
            clwc = np.insert(clwc, 0, 0.0)
            ciwc = np.insert(ciwc, 0, 0.0)
            cf = np.insert(cf, 0, 0.0)

        # clwc, p have N values
        # cwp has N-1 values
        clwp = clwc2cwp_simple(clwc, p)
        ciwp = clwc2cwp_simple(ciwc, p)
        cwp = clwp + ciwp
        th = 1e-3  # threshold cloud/clear-sky
        iscloud = cwp > th
        if any(iscloud):
            logger.debug("There's a cloud")
            layers = np.argwhere(iscloud)[:, 0]
            cfi = __intermediate(cf)[iscloud]
            fi = ciwp / cwp
            cld = rrtm.input.Cloud(
                2, layers, cfi, cwp[iscloud], fi[iscloud], radice=10.0
            )
    return atm, cld
