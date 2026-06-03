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


rrtm/output.py
"""

from collections import namedtuple
import logging

logger = logging.getLogger("rrtm:output")
logger.setLevel(logging.WARNING)

RRTM_data_block = namedtuple("RRTM_data_block", "w1 w2 level pressure upward_flux downward_flux net_flux heating_rate")

def __get_float(s):
    try:
        f = float(s)
    except ValueError as e:
        logging.warning(f'float conversion error ({s}), {e}')
        f = float('nan')
        # input("Press Enter to continue...")
    return f

def read(datafile_path="OUTPUT_RRTM"):
    blocks = []
    b = None
    with open(datafile_path, 'r') as f:
        for x in f:
            x = x.strip()
            p = x.split()
            if x.startswith("Wavenumbers"):
                # we are in a header block
                w1 = float(p[1])
                w2 = float(p[3])
                f.readline()
                f.readline()
                b = RRTM_data_block(w1, w2, [], [], [], [], [], [])
            elif len(x)==0:
                pass
            elif x.startswith("Modules"):
                break
            else:
                if b is None:
                    # if we are here we didn't go through a header block
                    raise Exception('Unrecognised file format', f'line="{x}", len(x)={len(x)}')
                try:
                    l = int(p[0])
                except (ValueError, IndexError):
                    logger.debug('Skipping unparsed output line: %r (len=%d)', x, len(x))
                    continue
                else:
                    b.level.append(l)
                    b.pressure.append(__get_float(p[1]))
                    b.upward_flux.append(__get_float(p[2]))
                    b.downward_flux.append(__get_float(p[3]))
                    b.net_flux.append(__get_float(p[4]))
                    b.heating_rate.append(__get_float(p[5]))
                    if l==0:
                        blocks.append(b)
                        b = None
                    
    return blocks
            

if __name__ == "__main__":
    data = read()
    print(f'loaded {len(data)} data blocks')
    for ii, b in enumerate(data):
        print(f'{ii:>2}: {b.w1:6.1f} - {b.w2:6.1f} cm-1')
        assert(min(b.level)==0)

    
