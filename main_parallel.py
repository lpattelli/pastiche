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


main_parallel.py:
Cycles through ERA5 datasets and calculates fluxes for each pixel and timesteps
saves results on netCDF4 files.
"""

import json
import datetime
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from netCDF4 import Dataset
from pathlib import Path
import era5
import export
import logger
from atmosphere import prepare_atmosphere, get_co2
from utils import __normalise_lon  # , get_description

import rrtm.input
import rrtm.run

SKIP_SEA_PIXELS = False


def process_one_timestep(atmosphere, cloud):
    """
    Calls the method to run RRTM
    and returns the results
    """
    try:
        rrtm_results = rrtm.run.run_and_read_results(atmosphere, cloud)
    except Exception as e:
        rrtm_results = None
        logger.error("RRTM error for %s", atmosphere.description)
        logger.error("%s", e)
    return rrtm_results


with open(sys.argv[1]) as config_file:
    config = json.load(config_file)
myplace = config["label"]
mydates = config["dates"]

dg = Dataset("in/era5/era5_geopotential_surface.nc", "r")  # time-invariant
dl = Dataset("in/era5/era5_land_sea_mask.nc", "r")  # time-invariant

mydates.sort()
started_outer = datetime.datetime.now()
for mydate in mydates:
    logger.info(mydate)
    if "TMY" in config:
        mm = mydate[4:6]
        pl_filename = f"in/era5/{myplace}/TMY/{mm}/pl_{mydate}.nc"
        sl_filename = f"in/era5/{myplace}/TMY/{mm}/surf_{mydate}.nc"
        outdir = f"out/{myplace}/TMY/{mm}"
    else:
        pl_filename = f"in/era5/{myplace}/pl_{mydate}.nc"
        sl_filename = f"in/era5/{myplace}/surf_{mydate}.nc"
        outdir = f"out/{myplace}"

    Path(outdir).mkdir(parents=True, exist_ok=True)
    ofname = Path(outdir).joinpath(f"{mydate}.nc")

    try:
        dp = Dataset(pl_filename, "r")
        ds = Dataset(sl_filename, "r")
    except FileNotFoundError as error:
        logger.warning("input file(s) not found %s", error)
        continue

    assert all(dp["latitude"][:] == ds["latitude"][:])
    assert all(dp["longitude"][:] == ds["longitude"][:])

    if "valid_time" in ds.variables:
        # the new era5 cds exports time as 'valid_time'
        # the format has also changed
        (Ntm, Nlt, Nln) = (
            len(ds["valid_time"]),
            len(ds["latitude"]),
            len(ds["longitude"]),
        )
        t = [datetime.datetime.utcfromtimestamp(x) for x in ds["valid_time"][:]]
    else:
        (Ntm, Nlt, Nln) = (len(ds["time"]), len(ds["latitude"]), len(ds["longitude"]))
        t = era5.time2datetime(ds["time"][:])

    if "pressure_level" in dp.variables:
        # the new era5 cds exports level as 'pressure_level'
        # they are also upside down
        levels = dp["pressure_level"][:]  # hPa
    else:
        levels = dp["level"][:]  # hPa

    N_BANDS = 17
    t0 = t[0]

    started = datetime.datetime.now()

    futures = []
    fut_map = {}
    with ProcessPoolExecutor(max_workers=16) as pool:
        for ilt in range(Nlt):
            rrtm_climatology = rrtm.input.atmospheric_model(ds["latitude"][ilt], t0)
            ig2_co2 = get_co2(ds["latitude"][ilt], t0)
            lt = ds["latitude"][ilt]
            for iln in range(Nln):
                ln = ds["longitude"][iln]
                (glt, gln), c2, d2 = era5.closest_gridpoint((lt, ln))
                assert ds["latitude"][ilt] == dg["latitude"][glt]
                assert ds["longitude"][iln] == __normalise_lon(dg["longitude"][gln])

                zs = dg["z"][0, glt, gln]
                lsm = dl["lsm"][0, glt, gln]
                if SKIP_SEA_PIXELS and (lsm < 0.5):
                    logger.info("lsm=%.2f skipping, (%d,%d)", lsm, ilt, iln)
                    continue

                for itm in range(Ntm):
                    description = f"ERA5 [{itm},{ilt},{iln}] t={t[itm]}, lat={lt}, lon={ln}, mdl={rrtm_climatology}"
                    # get_description(
                    #     ds, itm, ilt, iln,
                    #     rrtm_climatology
                    # )

                    ps = ds["sp"][itm, ilt, iln]
                    z = dp["z"][itm, :, ilt, iln]  # m**2 s**-2
                    T = dp["t"][itm, :, ilt, iln]
                    q = dp["q"][itm, :, ilt, iln]
                    o3 = dp["o3"][itm, :, ilt, iln]
                    cf = dp["cc"][itm, :, ilt, iln]
                    clwc = dp["clwc"][itm, :, ilt, iln]
                    ciwc = dp["ciwc"][itm, :, ilt, iln]
                    Ts = ds["skt"][itm, ilt, iln]
                    T2m = ds["t2m"][itm, ilt, iln]

                    atm, cld = prepare_atmosphere(
                        description,
                        ps,
                        levels,
                        zs,
                        z,
                        Ts,
                        T,
                        q,
                        o3,
                        rrtm_climatology,
                        ig2_co2,
                        cf=cf,
                        clwc=clwc,
                        ciwc=ciwc,
                        T2m=T2m,
                        add_one_at_surf=True,
                    )

                    fut = pool.submit(
                        process_one_timestep,
                        atm,
                        cld,
                    )
                    fut_map[fut] = (ilt, iln, itm)
                    futures.append(fut)

    Xsd = np.ma.MaskedArray(np.empty((Ntm, Nlt, Nln, N_BANDS)), mask=True)
    Xsu = np.ma.MaskedArray(np.empty((Ntm, Nlt, Nln, N_BANDS)), mask=True)
    Xsn = np.ma.MaskedArray(np.empty((Ntm, Nlt, Nln, N_BANDS)), mask=True)
    Xtu = np.ma.MaskedArray(np.empty((Ntm, Nlt, Nln, N_BANDS)), mask=True)
    for f in as_completed(futures):
        X = f.result()
        (ilt, iln, itm) = fut_map[f]
        if X is None:
            logger.warning("Empty result for (%d,%d,%d)", itm, ilt, iln)
            continue

        for bb in range(N_BANDS):
            Y = X[bb]
            Xsd[itm, ilt, iln, bb] = Y.downward_flux[-1]
            Xsu[itm, ilt, iln, bb] = Y.upward_flux[-1]
            Xsn[itm, ilt, iln, bb] = Y.net_flux[-1]
            Xtu[itm, ilt, iln, bb] = Y.upward_flux[0]

    if np.all(Xsd.mask):
        logger.warning("all values of Xsd are masked")
    # logger.debug(Xsd[:,10,10,0])

    export.nc(sl_filename, ofname, Xsd, Xsu, Xsn, Xtu)
    export.my_plots(ofname, ilt, iln, myplace, mydate)
    dp.close()
    ds.close()

    elapsed = datetime.datetime.now() - started
    logger.info("elapsed: %s", elapsed)

elapsed_outer = datetime.datetime.now() - started_outer
logger.info("elapsed main loop: %s", elapsed_outer)
