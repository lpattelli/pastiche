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


rrtm/input.py
"""

from collections import namedtuple
import numpy as np
import os
import logging

logger = logging.getLogger("rrtm:input")
logger.setLevel(logging.WARNING)

Atmosphere = namedtuple(
    "Atmosphere",
    "description z p jcharp T jchart Ts vmolk jchark model semiss",
    defaults=(1, 1.0),
)
Cloud = namedtuple("Cloud", "inflag layers cf cwp fracice radice")

def atmospheric_model(latitude, date):
    # 1  tropical model
    # 2  midlatitude summer model
    # 3  midlatitude winter model
    # 4  subarctic summer model
    # 5  subarctic winter model
    issummer = date.month > 4 and date.month <= 10
    if latitude < 0:
        issummer = not issummer
    alt = abs(latitude)
    if alt <= 30:
        model = 1
    elif alt <= 60:
        model = 2 if issummer else 3
    else:
        model = 4 if issummer else 5
    return model

def write(atm, infile_path="./MY_INPUT_RRTM", cld=None, use_pressure=False):
    hdr = ['0        1         2         3         4         5         6         7         8         9\n', \
           '123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789\n']
    spc = ' '
    user_model = 0
    iatm = '1'
    numangs = '3'
    ixsect = '0'
    iscat = '0'
    nbands = 99 #0 # 1-16 or 0: 10-3250 cm-1, 99: 17 spectral intervals

    ireflect = '0' # Lambertian reflection
    if np.isscalar(atm.semiss):
        iemis = '1' # each band has the same emissivity
        semiss = [atm.semiss]
        logger.debug('iemis: %s', iemis)
        logger.debug('semiss: %s', semiss)
    elif len(atm.semiss)==1:
        iemis = '1' # each band has the same emissivity
        semiss = [atm.semiss[0]]
        logger.debug('iemis: %s', iemis)
        logger.debug('semiss: %s', semiss)
    elif len(atm.semiss)==16:
        iemis = '2' # band-dependent surface emissivity
        semiss = atm.semiss
        logger.debug('iemis: %s', iemis)
        logger.debug('semiss: %s', semiss)
    else:
        logger.error('semiss should be either a scalar or 16-value array: %s', atm.semiss)
        raise ValueError
    
    if cld is None:
        icld = '0'
    else:
        # icld = '1' # random overlap assumption
        icld = '2' # maximum/random overlap assumption
        wdir = os.path.dirname(infile_path)        
        __write_in_cld(cld, in_cld_path=os.path.join(wdir,'IN_CLD_RRTM'))
        
    if use_pressure:
        ibmax = -1*len(atm.p) # <0 boundaries to be specified in pressure [mbars]
        upper_layers = [] #[0.5, 0.2, 0.1, 0.05, 0.02]
        immax = -1*(len(atm.p) + len(upper_layers))
    else:
        ibmax = len(atm.z) # >0 boundaries in [km]
        upper_layers = range(round(atm.z[-1]+1), 100, 10)
        immax = len(atm.z) +len(upper_layers)
        
    with open(infile_path, 'w') as f:
        f.writelines(hdr)
        f.write(f'$ {atm.description}\n') # RECORD 1.1
        f.writelines((spc*49, iatm, spc*19, ixsect, spc*12, iscat, spc, numangs, f'{nbands:>5}', spc*4, icld, '\n')) # R1.2
        f.writelines((f'{atm.Ts:10.2F}', spc, iemis, spc*2, ireflect)) # R1.4 
        for ii in range(len(semiss)):
            # f.write(f'{semiss[ii]:.4f}'.lstrip('0')) # R1.4
            f.write(f'{semiss[ii]:7.4f}') # R1.4
        f.write('\n')    
        f.writelines((f'{user_model:>5}',spc*5, f'{ibmax:>5}', '    1    1    7    0', '\n')) # R3.1

        if ibmax<0:
            zbnd = atm.p
            pbound = atm.p[0]
            ptoa = 1 #.2 #.1
            r3_2 = f'{pbound:10.3E}{ptoa:10.3E}\n'
        else:
            zbnd = atm.z
            hbound = atm.z[0]
            htoa = 100.0
            r3_2 = f'{hbound:10.3E}{htoa:10.3E}\n'
        f.write(r3_2) # R3.2

        if ibmax==0:
            f.write(f'{0:10d}{0:10d}{0:10d}{0:10d}{0:10d}\n') # R3.3A
        else:
            # either pressure or altitude levels
            # depending on whether ibmax < 0 or > 0
            for ii in range(0, len(zbnd), 8):
                f.writelines([f'{z:10.3E}' for z in zbnd[ii:ii+8]]) # R3.3B
                f.write('\n')

        f.write(f'{immax:>5} {atm.description}\n') # R3.4
        for ii in range(len(atm.p)): # loop on levels
            f.writelines(
                (f'{atm.z[ii]:10.3E}', f'{atm.p[ii]:10.3E}', f'{atm.T[ii]:10.3E}', spc*5,
                  f'{atm.jcharp}{atm.jchart}', spc*3, f'{atm.jchark}')) # R3.5
            for jj in range(len(atm.vmolk)): # loop on molecules
                if jj%8==0: f.write('\n')
                f.write(f'{atm.vmolk[jj][ii]:10.3E}') # R3.6.1 ... R3.6.N
            f.write('\n')
            
        for xx in upper_layers:
            if immax<0:
                # xx in pressure
                f.writelines((spc*10, f'{xx:10.3f}', spc*15, f'{atm.model}'*2 , spc*3, f'{atm.model}A'+ f'{atm.model}'*5, '\n')) # R3.5
            else:
                # xx in altitude
                f.writelines((f'{xx:10.3f}', spc*25, f'{atm.model}'*2 , spc*3, f'{atm.model}A'+ f'{atm.model}'*5, '\n')) # R3.5
            f.write(f'{0.:10.3f}{423.3:10.3E}\n')
        
        f.write('%%%%%\n\n')

        
def __write_in_cld(cld, in_cld_path='./MY_IN_CLD_RRTM'):
    hdr = ['0        1         2         3         4         5         6\n', \
           '123456789-123456789-123456789-123456789-123456789-123456789-\n']
    spc = ' '
    inflag = cld.inflag
    iceflag = 0
    liqflag = 0
    radice = cld.radice
    with open(in_cld_path, 'w') as f:
        f.writelines((spc*3,f'{inflag:2d}',
                      spc*4, f'{iceflag:1d}',
                      spc*4, f'{liqflag:1d}', '\n')) # RC1.1
        for ii, ll in enumerate(cld.layers):
            lay = ll+1
            cldfrac = cld.cf[ii]
            my_cwp = cld.cwp[ii]
            fracice = cld.fracice[ii]
            f.writelines((spc*2,f'{lay:3d}{cldfrac:10.3E}{my_cwp:10.3E}{fracice:10.3E}{radice:10.3E}','\n')) # RC1.2
        f.write('%\n')
        f.writelines(hdr)
    pass

if __name__ == "__main__":
    import numpy as np
    zpt = np.loadtxt("my_iccrim_zpt")
    m_atm = np.loadtxt("my_iccrim_atm")
    gases = [m_atm[:,mm] for mm in range(m_atm.shape[1])]
    
    atm = Atmosphere('SPECTRE ICRCCM experiment profile sonde, ozone sonde, raman, RASS', zpt[:,0], zpt[:,1], 'A', zpt[:,2], 'B', 290.93, gases, 'CAAAAA6', 6)
    write(atm)
    
    
