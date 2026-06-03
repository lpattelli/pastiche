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


pl_parallel_request_single_day_colab.py
"""

import sys
import os
import datetime
import cdsapi
import multiprocessing.dummy
import json


api_url = os.environ.get('CDSAPI_URL')
api_key = os.environ.get('CDSAPI_KEY')

USE_TMY = True
ALL_HOURS = [f"{x:02d}:00" for x in range(24)]


def pl(d):
    """
    pressure levels request
    """

    global config
    year = d.year
    month = d.month
    day = d.day
    c = cdsapi.Client(key=f"{api_key}", url=api_url)
    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": "reanalysis",
        "format": "netcdf",
        "variable": [
            "geopotential",
            "ozone_mass_mixing_ratio",
            "relative_humidity",
            "specific_humidity",
            "temperature",
            "fraction_of_cloud_cover",
            "specific_cloud_ice_water_content",
            "specific_cloud_liquid_water_content",
        ],
        "pressure_level": [
            "1",
            "2",
            "3",
            "5",
            "7",
            "10",
            "20",
            "30",
            "50",
            "70",
            "100",
            "125",
            "150",
            "175",
            "200",
            "225",
            "250",
            "300",
            "350",
            "400",
            "450",
            "500",
            "550",
            "600",
            "650",
            "700",
            "750",
            "775",
            "800",
            "825",
            "850",
            "875",
            "900",
            "925",
            "950",
            "975",
            "1000",
        ],
        "year": str(year),
        "month": f"{month:02}",
        "day": f"{day:02}",
        "time": config["hours"],
        "area": config["area"],
        "data_format": "netcdf_legacy",
    }
    target = (
        f"./in/era5/{config['label']}/TMY/{month:02}/pl_{year}{month:02}{day:02}.nc"
    )

    # do nothing if file was downloaded already
    if not os.path.isfile(target):
        c.retrieve(dataset, request, target)
    else:
        print(target + " exists: skipping")


if __name__ == "__main__":
    """
    source ./.venv/bin/activate
    python3 pl_parallel_request_single_day.py ./config/las_vegas.json
    """

    with open(sys.argv[1]) as config_file:
        config = json.load(config_file)
    for month in range(1, 13):
        os.makedirs(f"./in/era5/{config['label']}/TMY/{month:02}", exist_ok=True)

    dates = config["dates"]

    ds = []
    for date in dates:
        try:
            ds.append(datetime.datetime.strptime(date, "%Y%m%d"))
        except ValueError:
            print(f"Ignoring invalid: {date}")
    with multiprocessing.dummy.Pool(5) as pool:
        result_pl = pool.map(pl, ds)
