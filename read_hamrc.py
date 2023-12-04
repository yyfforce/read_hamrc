#!/usr/bin/env python

'''
Authors: Yuefeng Yin
This code is revised to be compatible with Python 3.*.

In this code, TB hamiltonian is stored in a dict called HR.hr.
All rvectors as key.

Delete most of the functions in the original code.
Only keeping (and revising) these:
    1. read HR to a dict
    2. copy a HR class instance
    3. write HR based on HR class

The name of the module (as a single file): read_hamrc.py
name for "Read Hamiltonian Compact version"

This code is designed to :
* Be easy to read
* Take no changes to raw data
* Be flexible for incorporating into other codes or further developing. 

NOTES ON ORIGIN:

This code is programmed by partly modifying the following code:
https://github.com/quanshengwu/wannier_tools/blob/master/utility/wannhr_symm/lib/read_hamr.py
(written in Python2 by Changming Yue, IOP (now at SUSTech))

For acknowledgement, please refer to the follwing link:
https://github.com/quanshengwu/wannier_tools



NOTES ON LICENCE APPLIED:

This code is distributed under GNU GPL v3.0 licence.
The full texts are available at:
https://www.gnu.org/licenses/gpl-3.0.en.html


This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


'''
#L1
import numpy as np

class HR():
    """ 
    The Wannier tight-binding (TB) Hamiltonian in real space.

    For the format of wannier90_hr.dat, please refer to the Wannier90 user guide.
    Can be accessed here:  https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf

    """

    def __init__(self, nwann, nrpt, irpt0, rpts, deg_rpt, hr):
        self.nwann=nwann
        self.nrpt = nrpt
        self.irpt0 = irpt0 
        self.deg_rpt = deg_rpt
        self.rpts = rpts
        self.hr = hr


    @staticmethod
    def from_file_dict(fname='wannier90_hr.dat'):
    
        """
        Generate a TB Hamiltonian from Wannier90 output file "case_hr.dat".

        Args:
            fname (string): the file that contains the Wannier90 output file 
            "case_hr.dat"

        This is a "initialization function" for the Class (see the return data type (HR)).
        """

        try: 
            with open(fname, 'r') as f:
                f.readline()
                nwann = int(f.readline().strip()) 
                nrpt = int(f.readline().strip())
                nline = np.int64(np.ceil(nrpt/15.0))
                tmp = []
                for i in range(nline):
                    tmp.extend(f.readline().strip().split())   
                tmp = [np.int64(item) for item in tmp]
                deg_rpt = np.array(tmp, dtype=np.int64)
                # read hr for each r-point
                rpts = np.zeros((nrpt,3), dtype=np.int64)
                hr={}
                nline_per_rvector = nwann*nwann
                for i in range(nrpt):
                    for j in range(nline_per_rvector):
                        rx, ry, rz, hr_i, hr_j, hr_real, hr_imag = f.readline().strip().split() 
                        if j==0:
                            hr[int(rx), int(ry), int(rz)] = np.zeros((nwann, nwann), dtype=np.complex64)
                        rpts[i,:] = int(rx), int(ry), int(rz) 
                        if int(rx) == 0 and int(ry) == 0 and int(rz) == 0:
                            irpt0 = i  
                        hr_i = int(hr_i) - 1
                        hr_j = int(hr_j) - 1
                        hr[int(rx),int(ry),int(rz)][hr_i,hr_j] = np.float64(hr_real) + np.float64(hr_imag) * 1j
                return HR(nwann, nrpt, irpt0, rpts, deg_rpt, hr)

        except IOError:
            print('read error')


    @staticmethod
    def copy_hr(other):

        """
        copy instance of HR
        """

        return HR(other.nwann, other.nrpt, other.irpt0, np.copy(other.rpts), np.copy(other.deg_rpt), np.copy(other.hr))


    def get_hr_output_dict(self,fname='wannier90_hr.dat'):
        """
        return a copy of wannier90_hr.dat
        """
        #print(self.rpts)
        with open(fname, 'w') as f:
            line="Deep copy of wannier90_hr.dat "+"\n"+"          "+ str(self.nwann) + "\n" + "        "+ str(self.nrpt) + "\n"
            f.write(line)
            nl = np.int32(np.ceil(self.nrpt/15.0))
            for l in range(nl):
                line="    "+'    '.join([str(np.int32(i)) for i in self.deg_rpt[l*15:(l+1)*15]])
                f.write(line)
                f.write('\n')
            for irpt in range(self.nrpt):
                rx = self.rpts[irpt,0];ry = self.rpts[irpt,1];rz = self.rpts[irpt,2]
                if irpt==0:
                    print(rx,ry,rz)
                for jatomorb in range(self.nwann):
                    for iatomorb in range(self.nwann):
                       rp =self.hr[rx,ry,rz][iatomorb,jatomorb].real
                       ip =self.hr[rx,ry,rz][iatomorb,jatomorb].imag
                       line="{:5d}{:5d}{:5d}{:5d}{:5d}{:20.6f}{:20.6f}\n".format(rx,ry,rz,iatomorb+1,jatomorb+1,rp,ip)    
                       f.write(line)

    