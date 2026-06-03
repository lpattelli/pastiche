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


export.py
"""

import os
import numpy as np
import convert
import era5
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import logger
from datetime import datetime
from myconstants import VERSION
from utils import relative_humidity_from_dewpoint


def __compute_scale_and_offset(min, max, n):
    # stretch/compress data to the available packed range
    scale_factor = (max - min) / (2**n - 1)
    # translate the range to be symmetric about zero
    add_offset = min + 2 ** (n - 1) * scale_factor
    return (scale_factor, add_offset)


def __copy_nc(src, dst):
    """
    Copies structure and data from another file
    """
    # to_include = ['longitude', 'latitude',
    #               'time', 'valid_time',
    #               't2m', 'skt',
    #               'tcc', 'tcwv',
    #               'msnlwrf', 'msnlwrfcs', # srf net lw
    #               'msnswrf', 'msnswrfcs', # srf net sw
    #               'mtnlwrf', 'mtnlwrfcs', # toa net lw
    #               'msdwlwrf', 'msdwlwrfcs', # srf down lw
    #               'msdwswrf', 'msdwswrfcs', # srf down sw
    #               ]
    to_include = [
        "longitude",
        "latitude",
        "time",
        "valid_time",
        "t2m",
        "skt",
        "tcc",
        "cbh",
        "tcwv",
        "u10",
        "v10",  # 10m U and V components
        "stl3",
        "stl4",  # soil temperature level 3 and 4
        "avg_snlwrf",
        "avg_snlwrfcs",  # srf net lw
        "avg_snswrf",
        "avg_snswrfcs",  # srf net sw
        "avg_tnlwrf",
        "avg_tnlwrfcs",  # toa net lw
        "avg_sdlwrf",
        "avg_sdlwrfcs",  # srf down lw
        "avg_sdswrf",
        "avg_sdswrfcs",  # srf down sw
    ]
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None)
        )

    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name in to_include:
            dst.createVariable(name, variable.datatype, variable.dimensions)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            dst[name][:] = src[name][:]


def __append_history(src, dst):
    # append history information
    now = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    history = f"{now} generated with PASTICHE {VERSION}: export.nc\n"
    history += src.getncattr("history")
    dst.setncattr("history", history)


def __add_flux_var(dst, flux, name, long_name):
    """Add calculated flux field"""
    if "valid_time" in dst.dimensions:
        name_of_time_dimension = "valid_time"
    else:
        name_of_time_dimension = "time"
    f = dst.createVariable(
        name,
        "i2",
        (name_of_time_dimension, "latitude", "longitude", "lw_bands"),
        fill_value=-32767,
    )
    f.long_name = long_name
    f.units = "W m**-2"
    scale_factor, add_offset = __compute_scale_and_offset(
        np.nanmin(flux), np.nanmax(flux), 16
    )
    if np.isfinite(scale_factor) and np.isfinite(add_offset):
        f.scale_factor = scale_factor
        f.add_offset = add_offset
        logger.debug("adding scale and offset")
    else:
        logger.debug("scale or offset are not finite")
    logger.debug("scale=%f offset=%f", scale_factor, add_offset)
    # logger.debug("flux[0,10,10,0]=%f", flux[0,10,10,0])
    f[:] = flux
    # logger.debug("f[0,10,10,0]=%f", f[0,10,10,0])


def __add_relative_humidity(src, dst):
    """Add calculated relative humidity field"""
    t2m = src["t2m"][:]
    td = src["d2m"][:]
    rh = 100.0 * relative_humidity_from_dewpoint(t2m, td)
    if "valid_time" in dst.dimensions:
        name_of_time_dimension = "valid_time"
    else:
        name_of_time_dimension = "time"
    r = dst.createVariable(
        "r", "i2", (name_of_time_dimension, "latitude", "longitude"), fill_value=-32767
    )
    r.long_name = (
        "Relative humidity calculated from 2 m temperature and 2 m dewpoint temperature"
    )
    r.standard_name = "relative humidity"
    r.units = "%"
    scale_factor, add_offset = __compute_scale_and_offset(
        np.nanmin(rh), np.nanmax(rh), 16
    )
    if np.isfinite(scale_factor) and np.isfinite(add_offset):
        r.scale_factor = scale_factor
        r.add_offset = add_offset
        logger.debug("adding scale and offset")
    else:
        logger.debug("scale or offset are not finite")
    logger.debug("scale=%f offset=%f", scale_factor, add_offset)
    r[:] = rh


def nc(
    src_file, dst_file, flux_dwn_at_srf, flux_up_at_srf, flux_net_at_srf, flux_up_at_toa
):
    """
    It first copies structure and data from another file;
    then adds spectral dimensions and data
    """
    with netCDF4.Dataset(src_file, "r") as src, netCDF4.Dataset(dst_file, "w") as dst:
        __copy_nc(src, dst)

        dst.createDimension("lw_bands", 17)  # 16 + total
        dst.createDimension("band_limits", 2)

        lw_bands = dst.createVariable("lw_bands", int, ("lw_bands",))
        lw_bands[:] = [x for x in range(17)]
        lw_band_limits = dst.createVariable(
            "lw_band_limits", "i2", ("lw_bands", "band_limits")
        )
        rrtm_bands = np.array(
            [
                [10, 3250],
                [10, 350],
                [350, 500],
                [500, 630],
                [630, 700],
                [700, 820],
                [820, 980],
                [980, 1080],
                [1080, 1180],
                [1180, 1390],
                [1390, 1480],
                [1480, 1800],
                [1800, 2080],
                [2080, 2250],
                [2250, 2380],
                [2380, 2600],
                [2600, 3250],
            ]
        )  # RRTM band limits
        lw_band_limits.units = "cm**-1"
        lw_band_limits[:] = rrtm_bands

        __add_flux_var(
            dst,
            flux_dwn_at_srf,
            "sd",
            "RRTM-calculated surface downward long-wave radiation flux",
        )
        __add_flux_var(
            dst,
            flux_up_at_srf,
            "su",
            "RRTM-calculated surface upward long-wave radiation flux",
        )
        __add_flux_var(
            dst,
            flux_net_at_srf,
            "sn",
            "RRTM-calculated surface net long-wave radiation flux",
        )
        __add_flux_var(
            dst,
            flux_up_at_toa,
            "tu",
            "RRTM-calculated TOA upward long-wave radiation flux",
        )
        __append_history(src, dst)

        if "d2m" in src.variables.keys():
            # if the 2m_dewpoint_temperature is present
            __add_relative_humidity(src, dst)

        # if 'valid_time' in dst.dimensions:
        #     backup_time = dst['valid_time'][:]
        #     dst.renameDimension('valid_time','time')
        #     dst.renameVariable('valid_time','time')
        #     dst['time'][:] = backup_time


def ridolfi_atm(p, T, q, o3, cc, ciwc, clwc, zs, ps, Ts, T2m, t, description, label):
    h2o = convert.mmr2vmr(convert.q2mmr(q), species="h2o")  # vmr [ppmv]
    o3 = convert.mmr2vmr(o3, species="o3")  # vmr [ppmv]
    X = np.column_stack((p, T, h2o, o3, cc, 1e3 * ciwc, 1e3 * clwc))

    hs = era5.z2h(zs)  # [m]
    ps /= 100  # [hPa]

    if p[0] > p[1]:
        # data are ordered from the ground up
        # we need to flip them
        X = np.flipud(X)

    ofname = os.path.join(
        "out/ridolfi", f"{label}_atm_{t.strftime('%Y%m%dT%H%M%S')}.dat"
    )
    hdr = description
    hdr += "\n"
    hdr += f"hs = {hs} [m]\n"
    hdr += f"ps = {ps} [hPa]\n"
    hdr += f"Ts = {Ts} [K]\n"
    hdr += f"T2m = {T2m} [K]\n"
    hdr += "\n"
    hdr += f"dims = {X.shape}\n"
    hdr += "p [hPa],T [K], h2o [ppmv], o3 [ppmv], cc [0-1], ciwc [g kg**-1], clwc [g kg**-1]"
    np.savetxt(ofname, X, header=hdr)
    return ofname


def ridolfi_fluxes(X, era_fluxes, era_fluxes_cs, t, description, label):
    Y = X[0]
    Z = np.column_stack(
        (Y.level, Y.pressure, Y.downward_flux, Y.upward_flux, Y.net_flux)
    )
    for Y in X[1:]:
        Z = np.column_stack((Z, Y.downward_flux, Y.upward_flux, Y.net_flux))

    hdr = description + "\n"
    hdr += "ERA5 fluxes all sky [Wm-2]:\n"
    hdr += f"Fsd = {era_fluxes[0]}\n"
    hdr += f"Fsu = {era_fluxes[1]}\n"
    hdr += f"Fsn = {era_fluxes[2]}\n"
    hdr += f"Ftu = {era_fluxes[3]}\n"
    hdr += "ERA5 fluxes clear-sky [Wm-2]:\n"
    hdr += f"Fsd_cs = {era_fluxes_cs[0]}\n"
    hdr += f"Fsu_cs = {era_fluxes_cs[1]}\n"
    hdr += f"Fsn_cs = {era_fluxes_cs[2]}\n"
    hdr += f"Ftu_cs = {era_fluxes_cs[3]}\n"
    hdr += f"dims = {Z.shape}, "
    hdr += "columns are:\n"
    hdr += "level, p [hPa], "
    hdr += "Band tot.: Downw. F [Wm-2], Upw. F [Wm-2], Net F [Wm-2], \n"
    hdr += "Band[1-16] Downw. F [Wm-2], Upw. F [Wm-2], Net F [Wm-2]\n"

    ofname = os.path.join(
        "out/ridolfi", f"{label}_fluxes_{t.strftime('%Y%m%dT%H%M%S')}.dat"
    )
    np.savetxt(ofname, Z, header=hdr)
    return ofname


def inrim_fluxes(
    hh,
    Xsd,
    Xsu,
    Xsn,
    Fsd,
    Fsu,
    Fsn,
    Fsdc,
    Fsuc,
    Fsnc,
    Ts,
    T2m,
    date,
    description,
    label,
):
    Z = np.column_stack((hh, Xsd, Xsu, Xsn, Fsd, Fsu, Fsn, Fsdc, Fsuc, Fsnc, Ts, T2m))
    hdr = description + "\n"
    hdr += f"dims = {Z.shape}, "
    hdr += "columns are:\n"
    hdr += "Time [hours since 1900-01-01 00:00], RRTM Fsd x17 bands, RRTM Fsu x 17 bands, RRTM Fsn x 17 bands\n"
    hdr += "ERA Fsd, ERA Fsu, ERA Fsn, "
    hdr += "ERA Fsd_cs, ERA Fsu_cs, ERA Fsn_cs,\n"
    hdr += "Ts [K], T2m [K]\n"
    hdr += "all fluxes in [Wm-2]"
    ofname = os.path.join("out/inrim", f"{label}_fluxes_{date}.dat")
    np.savetxt(ofname, Z, header=hdr)


def my_plots(ncfile, ilt, iln, myplace, mydate):
    ts = 0.5
    ds = netCDF4.Dataset(ncfile, "r")
    t = ds["time"][:]
    t -= t[0]
    # lt = ds['latitude'][ilt]
    # ln = ds['longitude'][iln]

    # Fsd = ds['msdwlwrf'][:]
    # Fsdc = ds['msdwlwrfcs'][:]
    # Fsn = -ds['msnlwrf'][:]
    # Fsnc = -ds['msnlwrfcs'][:]
    Fsd = ds["avg_sdlwrf"][:]
    Fsdc = ds["avg_sdlwrfcs"][:]
    Fsn = -ds["avg_snlwrf"][:]
    Fsnc = -ds["avg_snlwrfcs"][:]
    Fsu = Fsd + Fsn
    Fsuc = Fsdc + Fsnc
    # Ftn = -ds['mtnlwrf'][:]
    # Ftnc = -ds['mtnlwrfcs'][:]
    # Ftn = -ds['avg_tnlwrf'][:]
    # Ftnc = -ds['avg_tnlwrfcs'][:]

    Xsd = ds["sd"][..., 0]
    Xsn = ds["sn"][..., 0]
    Xsu = ds["su"][..., 0]
    # Xtu = ds['tu'][..., 0]
    plt.rcParams.update({"font.size": 14})

    f, ax = plt.subplots(3, 1, figsize=(12.8, 9.6))
    # f.suptitle(f'{ncfile}, ({lt},{ln}), [{ilt},{iln}]')
    ax[0].plot(
        t - ts,
        Fsd[:, ilt, iln],
        label="ERA downwelling @surface",
        marker="s",
        linestyle="none",
    )
    ax[0].plot(
        t - ts,
        np.interp(t - ts, t, Xsd[:, ilt, iln], left=np.nan, right=np.nan),
        label="calculated downwelling @surface",
        marker="o",
        color="C1",
    )
    ax[0].plot(t - ts, Fsdc[:, ilt, iln], label="ERA clear sky", color="C2")
    ax[0].set_ylabel("Fdown / Wm**-2")

    ax[1].plot(
        t - ts,
        Fsu[:, ilt, iln],
        label="ERA upwelling @surface",
        marker="s",
        color="C0",
        linestyle="none",
    )
    ax[1].plot(
        t - ts,
        np.interp(t - ts, t, Xsu[:, ilt, iln], left=np.nan, right=np.nan),
        label="calculated upwelling @surface",
        marker="o",
        color="C1",
    )
    ax[1].plot(t - ts, Fsuc[:, ilt, iln], label="ERA clear sky", color="C2")
    # ax[1].plot(t-ts, B, label="Planck_i", linestyle='None', marker='.', markerfacecolor='None')
    ax[1].set_ylabel("Fup / Wm**-2")

    ax[2].plot(
        t - ts,
        Fsn[:, ilt, iln],
        label="-ERA net @surface",
        marker="s",
        color="C0",
        linestyle="none",
    )
    ax[2].plot(
        t - ts,
        np.interp(t - ts, t, Xsn[:, ilt, iln], left=np.nan, right=np.nan),
        label="calculated net @surface",
        marker="o",
        color="C1",
    )
    ax[2].plot(t - ts, Fsnc[:, ilt, iln], label="-ERA clear sky", color="C2")
    ax[2].set_ylabel("Fnet / Wm**-2")

    # ax[3].plot(t-ts, Ftn[:, ilt, iln], label="-ERA tn", marker='.', color='grey')
    # ax[3].plot(t-ts, np.interp(t-ts, t, Xtu[:, ilt, iln], left=np.nan, right=np.nan),
    #            label="Xtui", marker='o', linestyle='none', color='green')
    # # ax[3].plot(t-ts, Ftnc[:, ilt, iln], label="-ERA tn cs")
    # ax[3].set_ylabel('top up / Wm**-2')
    ax[2].set_xlabel("time / hour")

    for aa in ax:
        aa.grid()
        aa.legend()

    plot_dir = os.path.join("figs", myplace, "fluxes")
    if not os.path.exists(plot_dir):
        logger.debug("creating directory: %s", plot_dir)
        os.makedirs(plot_dir)
    f.savefig(os.path.join(plot_dir, f"{mydate}.png"))
    plt.close(f)

    # f,ax = plt.subplots(1, 4, figsize=(12.8,9.6))
    # ax[0].scatter(Fsd[~Xsd.mask].ravel(), Xsd[~Xsd.mask].ravel(), marker='.', color='blue', alpha=0.01)
    # ax[0].axline((1,1), slope=1)

    # ax[1].scatter(Fsu[~Xsu.mask].ravel(), Xsu[~Xsu.mask].ravel(), marker='.', color='red', alpha=0.01)
    # ax[1].axline((1,1), slope=1)

    # ax[2].scatter(Fsn[~Xsn.mask].ravel(), Xsn[~Xsn.mask].ravel(), marker='.', color='green', alpha=0.01)
    # ax[2].axline((1,1), slope=1)

    # ax[3].scatter(Ftn[~Xtu.mask].ravel(), Xtu[~Xtu.mask].ravel(), marker='.', color='grey', alpha=0.01)
    # ax[3].axline((1,1), slope=1)

    # for aa in ax:
    #     aa.grid()
    #     aa.legend()

    # plt.close(f)

    # f,ax = plt.subplots(1, 4, figsize=(12.8,9.6))
    # ax[0].hist2d(Fsd[~Xsd.mask].ravel(), Xsd[~Xsd.mask].ravel(), cmap='Blues', bins=100)
    # ax[0].axline((1,1), slope=1)

    # ax[1].hist2d(Fsu[~Xsu.mask].ravel(), Xsu[~Xsu.mask].ravel(), cmap='Reds', bins=100)
    # ax[1].axline((1,1), slope=1)

    # ax[2].hist2d(Fsn[~Xsn.mask].ravel(), Xsn[~Xsn.mask].ravel(), cmap='Greens', bins=100)
    # ax[2].axline((1,1), slope=1)

    # ax[3].hist2d(Ftn[~Xtu.mask].ravel(), Xtu[~Xtu.mask].ravel(), cmap='Greys', bins=100)
    # ax[3].axline((1,1), slope=1)

    # for aa in ax:
    #     aa.grid()
    #     aa.legend()

    # plt.close(f)

    # f,ax = plt.subplots(3,1, figsize=(12.8,9.6))
    # f.suptitle(f'{mydate}, {myplace}')
    # ax[0].plot(hh-hh[0], ds['skt'][:, ilt, iln], label="skt")
    # ax[0].plot(hh-hh[0], ds['t2m'][:, ilt, iln], label="t2m")
    # ax[0].plot(hh-hh[0], dp['t'][:, -1, ilt, iln], label="t[-1]")
    # ax[0].set_ylabel('T / K')
    # ax[0].legend()

    # ax[1].plot(hh-hh[0]-ts, planck.inv_stefan_boltzmann(Xsd[:,0]), label='Tsd')
    # ax[1].plot(hh-hh[0]-ts, planck.inv_stefan_boltzmann(Xsdi), label='Tsdi')
    # ax[1].plot(hh-hh[0]-ts, planck.inv_stefan_boltzmann(Fsd), label='T ERA sd')
    # ax[1].set_ylabel('Tsky / K')
    # ax[1].set_xlim(ax[0].get_xlim())
    # ax[1].legend()

    # rl = ax[2].plot(hh-hh[0], dp['r'][:, -1, ilt, iln], label="r")
    # ax[2].set_ylabel('rh / %')
    # axq = ax[2].twinx()
    # ql = axq.plot(hh-hh[0], dp['q'][:, -1, ilt, iln], label="q", color='C1')
    # axq.set_ylabel('q / kgkg**-1')
    # lns = rl+ql
    # labs = [l.get_label() for l in lns]
    # ax[2].legend(lns, labs, loc=0)
    # ax[2].set_xlabel('time / hour')

    # for aa in ax:
    #     aa.grid()

    # f.savefig(os.path.join('figs', myplace['label'], 'temps', f'{mydate}.png'))
    # plt.close(f)

    # f, ax = plt.subplots(2,1, figsize=(12.8,9.6))
    # f.suptitle(f'{mydate}, {myplace}')
    # # ax[0].plot(hh-hh[0], ds['tcwv'][:, ilt, iln], label="col. wv")
    # # ax[0].set_ylabel("tcwv / kg m**-2")
    # ax[0].plot(hh-hh[0], ds['tclw'][:, ilt, iln], label="cloud lw")
    # ax[0].plot(hh-hh[0], ds['tciw'][:, ilt, iln], label="cloud iw")
    # ax[0].set_ylabel("tc[il]w / kgm**-2")
    # ax[1].bar(hh-hh[0], ds['tcc'][:, ilt, iln], label="cloud fr.")
    # ax[1].set_ylim(0, 1)
    # ax[1].set_ylabel("cloud fraction")
    # acb = ax[1].twinx()
    # p2, = acb.plot(hh-hh[0], ds['cbh'][:, ilt, iln], color='C1', label='cloud bh')
    # acb.set_ylabel("cloud base / m")
    # ax[1].set_xlabel('time / hour')

    # acb.tick_params(axis='y', colors=p2.get_color())
    # for aa in ax:
    #     aa.grid()

    # f.savefig(os.path.join('figs', myplace['label'], 'clouds', f'{mydate}.png'))
    # plt.close(f)
