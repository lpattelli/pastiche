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


rrtm/tests/test_rrtm_input.py
"""

import unittest
from unittest.mock import patch, mock_open

import rrtm.input


class TestWriteRRTMInputFile(unittest.TestCase):

    def setUp(self):
        z = [0.0, 10.338, 20.602]
        p = [973.60, 248.40, 47.61]
        T = [17.78, -44.61, -59.28]
        Ts = 290.93
        gases = [[6.491E+00, 3.036E-02, 2.790E-03],
                 [3.500E+02, 3.500E+02, 3.500E+02],
                 [2.668E-02, 7.557E-02, 2.336E+00],
                 [3.089E-01, 3.060E-01, 3.060E-01],
                 [1.892E-01, 1.430E-01, 1.430E-01],
                 [1.836E+00, 1.700E+00, 1.700E+00]]
        
        self.atm = rrtm.input.Atmosphere('SPECTRE ICRCCM experiment profile sonde, ozone sonde, raman, RASS', z, p, 'A', T, 'B', Ts, gases, 'CAAAAA6', 6, semiss=1.)
        fake_file_path = "fake/file/path"
        with patch('rrtm.input.open', mock_open()) as mocked_file:
            rrtm.input.write(self.atm, fake_file_path)
            
            # assert if opened file on write mode 'w'
            mocked_file.assert_called_once_with(fake_file_path, 'w')
            self.magic_mock = mocked_file()
            # assert if write(content) was called from the file opened
            # in other words, assert if the specific content was written in file
            #mocked_file().write.assert_called_once_with(atm)
    

    def test_file_writer(self):
        a_call_list = self.magic_mock.write.call_args_list
        b_call_list = self.magic_mock.writelines.call_args_list
        self.assertEqual(self.atm.Ts, float(b_call_list[2].args[0][0]))
        n_gases = len(self.atm.vmolk)
        n_levels = len(self.atm.vmolk[0])
        for ii in range(n_levels):
            for jj in range(n_gases):
                a_call_args = a_call_list[7 + ii * (n_gases + 2) + jj].args[0]
                self.assertEqual(self.atm.vmolk[jj][ii], float(a_call_args))
        
        for ii in range(3):
            a_call = self.magic_mock.writelines.call_args_list[ii+5]
            a_call_args = a_call.args[0]
            self.assertAlmostEqual(self.atm.z[ii], float(a_call_args[0]), places=2)
            self.assertEqual(self.atm.p[ii], float(a_call_args[1]))
            self.assertEqual(self.atm.T[ii], float(a_call_args[2]))
            self.assertEqual(' '*5, a_call_args[3])
            self.assertEqual('AB', a_call_args[4])
            self.assertEqual(' '*3, a_call_args[5])
            self.assertEqual('CAAAAA6', a_call_args[6])
    # def test_file_path(self):
    #     self.file_path.assert_called_once_with(self.fake_file_path, 'r')

        
