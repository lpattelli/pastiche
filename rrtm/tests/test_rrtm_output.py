"""
Copyright (C) 2025
Claudio Belotti (CNR-INO, Italy), claudio.belotti@cnr.it
Lorenzo Pattelli (INRiM, Italy), l.pattelli@inrim.it


    This file is part of PASTICHE.

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


rrtm/tests/test_rrtm_output.py
"""

import unittest
from unittest.mock import patch, mock_open

import rrtm.output


class TestReadRRTMOutputFile(unittest.TestCase):
    def setUp(self):

        file_content_mock = """ Wavenumbers:   10.0 - 3250.0 cm-1, ATM      1
 LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NETFLUX       HEATING RATE
             mb          W/m2          W/m2           W/m2          degree/day
   3      918.2        405.2754      297.4865       107.7889836           -1.82311
   2      948.5        407.6543      308.3956        99.2587188           -2.37758
   1      966.8        405.7214      310.5301        95.1912704           -1.87381
   0      973.6        405.0444      310.9624        94.0819309           -1.36983

 Wavenumbers:   10.0 -  350.0 cm-1, ATM      1
 LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NETFLUX       HEATING RATE
             mb          W/m2          W/m2           W/m2          degree/day
   3      918.2         53.5015       53.4717         0.0297738           -0.00172
   2      948.5         53.9346       53.9352        -0.0005956           -0.00846
   1      966.8         53.7083       53.7311        -0.0227927           -0.01023
   0      973.6         53.5944       53.6062        -0.0118274            0.01354

  Modules and versions used in this calculation:

                rrtm.f:  $Revision: 3.6 $                     rrtatm.f:  $Revision: 9.4 $  
                 rtr.f:  NOT USED                             rtrcld.f:  NOT USED"""

        self.fake_file_path = 'file/path/mock'

        ## see here
        ## https://github.com/otrabalhador/python-testing-by-examples/blob/master/docs/en/mocking/examples/reading-writing-on-files.md
        ## for details
        with patch('rrtm.output.open'.format(__name__),
                   new=mock_open(read_data=file_content_mock)) as _file:
            self.actual = rrtm.output.read(self.fake_file_path)
            self.file_path = _file


    def test_file_path(self):
        self.file_path.assert_called_once_with(self.fake_file_path, 'r')

        
    def test_number_of_read_blocks(self):
        self.assertEqual(2, len(self.actual))

        
    def test_band_limits(self):
        self.assertEqual(10.0, self.actual[0].w1)
        self.assertEqual(3250.0, self.actual[0].w2)
        self.assertEqual(10.0, self.actual[1].w1)
        self.assertEqual(350.0, self.actual[1].w2)

        
    def test_read_level_and_p_values(self):
        expected_p_values = [918.2, 948.5, 966.8, 973.6]
        for aa in self.actual:
            self.actual_p_values = aa.pressure
            for ep, ap in zip(expected_p_values, self.actual_p_values):
                self.assertEqual(ep,ap)

            self.assertEqual(3, max(aa.level))
            self.assertEqual(0, min(aa.level))


    def test_fluxes_values(self):
        block_index = 1 # 2nd block
        level_index = 1 # 2nd level
        fields = ['upward_flux', 'downward_flux', 'net_flux', 'heating_rate']
        values = [53.9346, 53.9352, -0.0005956, -0.00846]
        for ii,fieldname in enumerate(fields):
            # iterate over fields and compare with expected
            self.assertEqual(values[ii], getattr(self.actual[block_index], fieldname)[level_index])
        
