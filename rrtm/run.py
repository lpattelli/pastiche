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


rrtm/run.py
"""

import logger
import os
import subprocess
import contextlib
import tempfile
import rrtm.input
import rrtm.output


RRTM_EX = os.environ.get('PASTICHE_RRTM_EXE', '/content/RRTM_LW/rrtm_v3.3.1_linux_ifx')
logger.info('using %s as RRTM executable', RRTM_EX)

def run_and_read_results(atm, cld=None, keep_temp_files=False):
    with (
            contextlib.nullcontext(tempfile.mkdtemp(suffix='_rrtm', dir='.'))
            if keep_temp_files
            else tempfile.TemporaryDirectory(suffix='_rrtm')
    ) as tmpdirname:
        logger.debug(tmpdirname)
        try: # mv previous results
            os.rename(os.path.join(tmpdirname, 'OUTPUT_RRTM'),
                      os.path.join(tmpdirname, 'OUTPUT_RRTM.bak'))
        except FileNotFoundError:
            pass

        if not os.path.exists(RRTM_EX):
            raise FileNotFoundError(
                f'RRTM executable not found: {RRTM_EX}. '
                'Set PASTICHE_RRTM_EXE to the path of the RRTM_LW executable.'
            )
        os.symlink(RRTM_EX, os.path.join(tmpdirname, 'rrtm'))
        rrtm.input.write(atm, infile_path=os.path.join(tmpdirname,'INPUT_RRTM'),
                         cld=cld, use_pressure=True)
        result = subprocess.run(['./rrtm'],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmpdirname, check=True)
        if result.stderr != b'':
            logger.error(result)
            raise RuntimeError('RRTM executable wrote to stderr; see log for details')
        # print(os.listdir(tmpdirname))
        try:
            X = rrtm.output.read(datafile_path=os.path.join(tmpdirname,'OUTPUT_RRTM'))
        except Exception as e:
            logger.error(e)
            raise
        # input("Press Enter to continue...")
    return X
