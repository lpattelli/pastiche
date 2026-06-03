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


ig2.py
"""

import pathlib
import re
import os
import datetime
import numpy as np


class IG2(object):
    altitudes = np.arange(120, 0 - 1, -1)
    n_altitudes = len(altitudes)
    latitude_bands = np.array(
        [
            [-90, -65],
            [-65, -20],
            [-20, 0],
            [0, 20],
            [20, 65],
            [65, 90],
        ]
    )
    n_bands = len(latitude_bands)

    def __init__(self):
        self.profiles = []

    def append(self, profile):
        self.profiles.append(profile)

    def __repr__(self):
        return f"IG2: {len(self.profiles)} profiles"

    def __str__(self):
        return f"IG2: {len(self.profiles)} profiles"

    def by_lat(self, lt, night=False):
        idx = np.where(
            (lt >= self.latitude_bands[:, 0]) & (lt <= self.latitude_bands[:, 1])
        )[0][0]
        if night:
            idx += 6
        return self.profiles[idx]

    @classmethod
    def from_file(cls, fname):
        yyyy = pathlib.Path(fname).parts[-3]

        profile_header = re.compile(
            r"# (\w+) profile \((\S+)\) for lat (-?\d+\.?0*) degrees, month (\d+\.?0*):"
        )

        with open(fname, "r") as fin:
            lines = fin.readlines()

        this_ig2 = IG2()
        profile = None
        npr = -1
        dn = None
        lines_iter = iter(lines)
        for line in lines_iter:
            line = line.strip()
            if line == "# DN:":
                nextline = next(lines_iter).split()
                dn = [bool(int(x)) for x in nextline]
            ph = profile_header.match(line)
            if ph:
                npr += 1
                night = (
                    (npr >= IG2.n_bands)
                    and (dn is not None)
                    and (dn[npr - IG2.n_bands])
                )
                profile = IG2.Profile(
                    species=ph.group(1),
                    units=ph.group(2),
                    latitude=float(ph.group(3)),
                    year=yyyy,
                    month=int(float(ph.group(4))),
                    night=night,
                )
                data = []
                continue
            if profile is not None:  # i.e. profile is not empty
                for x in line.split():
                    data.append(float(x))
                if len(data) == IG2.n_altitudes:
                    profile.data = np.array(data)
                    if (
                        npr < IG2.n_bands
                    ) or profile.night:  # add night profile only if meaningful
                        this_ig2.append(profile)
                    profile = None
        return this_ig2

    @staticmethod
    def by_date(d, species, archive="./ig2"):
        """
        d yyyymmdd
        TODO d could be datetime
        """
        if type(d) is datetime.datetime:
            y = d.year
            m = d.month
        else:
            d = str(d)
            y = int(d[:4])
            m = int(d[4:6])
        T = [
            "january",
            "january",
            "april",
            "april",
            "april",
            "july",
            "july",
            "july",
            "october",
            "october",
            "october",
            "january",
        ]
        if y > 2020:
            y = 2020
        elif y < 2002:
            y = 2002
        p = os.path.join(archive, str(y), T[m - 1], f"IG_{species.upper()}.DAT")
        return IG2.from_file(p)

    class Profile(object):
        def __init__(
            self, species, units, latitude, year, month, night=None, data=None
        ):
            self.species = species
            self.units = units
            self.latitude = latitude
            self.year = year
            self.month = month
            self.night = night
            self.data = data

        # def __str__(self):
        #     return f'{self.species} profile ({self.units}) for lat {self.latitude} degrees, year {self.year}, month {self.month}'

        def __repr__(self):
            s = f"{self.species} profile ({self.units}, {len(self.data)} values) for lat {self.latitude} degrees, year {self.year}, month {self.month}"
            if self.night:
                s += ", night"
            s += "."
            return s


if __name__ == "__main__":
    fname = "in/ig2/2012/april/IG_CO2.DAT"
    my_ig2 = IG2.from_file(fname)

    print(my_ig2)
    print(my_ig2.profiles)
