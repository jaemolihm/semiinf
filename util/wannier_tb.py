"""
This module reads the tight-binding Hamiltonian calculated by wannier90.x,
and calculates the band structure.
Note that hr[ir,i,j] = <w_i,R=0|H|w_j,R=ir>
"""
import os
import sys
import copy
import pickle
import numpy as np
import matplotlib.pyplot as plt

class TBdict:

    def __init__(self, nw=None, nrpts=0, hr=None, rvec=None, ndegen=None):
        if hr is not None:
            if hr.shape != (nrpts, nw, nw):
                raise ValueError('hr.shape must be (nrpts, nw, nw)')
        if rvec is not None:
            if rvec.shape != (3, nrpts):
                raise ValueError('rvec.shape must be (3, nrpts)')
        if ndegen is not None:
            if ndegen.shape != (nrpts,):
                raise ValueError('ndegen.shape must be (nrpts,)')
        self.nw = nw
        self.nrpts = nrpts
        self.hr = hr
        self.rvec = rvec
        self.ndegen = ndegen

    @staticmethod
    def read_array(filename, dims, dtype=np.float):
        with open(filename, 'rb') as f:
            array_out = np.fromfile(f, dtype=dtype).reshape(dims, order='F')
        return array_out

    @classmethod
    def read_hr_from_dat(cls, filename):
        """Read TBdict object from _hr.dat formatted output of Wannier90"""
        with open(filename, 'r') as f:
            f.readline()
            nw = int(f.readline().split()[0])
            nrpts = int(f.readline().split()[0])
            ndegen = np.zeros((nrpts), dtype=int)
            for ir in range(1+(nrpts-1)//15):
                ndegen_add = f.readline().split()
                ndegen_add = np.array([int(x) for x in ndegen_add])
                if ir == (nrpts-1)//15: # last line
                    ndegen[ir*15:] = ndegen_add.copy()
                else:
                    ndegen[ir*15:(ir+1)*15] = ndegen_add.copy()

            rvec = np.zeros((3, nrpts), dtype=int)
            hr = np.zeros((nrpts, nw, nw), dtype=complex)
            for ir in range(nrpts):
                for iw in range(nw):
                    for jw in range(nw):
                        data = f.readline().split()
                        assert jw == int(data[3])-1
                        assert iw == int(data[4])-1
                        if (iw == 0) & (jw == 0):
                            rvec[:,ir] = [int(x) for x in data[0:3]]
                        hr[ir,jw,iw] = float(data[5]) + 1j*float(data[6])
            return cls(nw=nw, nrpts=nrpts, ndegen=ndegen, rvec=rvec, hr=hr)

    @classmethod
    def read_hr_from_bin(cls, seedname):
        """Read TBdict object from binary output of Wannier90"""
        rvec = cls.read_array(seedname + '_irvec.bin', (3, -1), int)
        nrpts = rvec.shape[1]
        ndegen = cls.read_array(seedname + '_ndegen.bin', (nrpts,), int)
        hr = cls.read_array(seedname + '_hr.bin', (-1, nrpts), complex)
        nw = int(np.sqrt(hr.shape[0]))
        hr = hr.reshape((nw, nw, nrpts))
        hr = np.moveaxis(hr, 2, 0)
        return cls(nw=nw, nrpts=nrpts, ndegen=ndegen, rvec=rvec, hr=hr)

    @classmethod
    def read_hr_from_pkl(cls, filename):
        """Read TBdict object from file"""
        with open(filename+'.pkl', 'rb') as handle:
            return pickle.load(handle)

    @classmethod
    def interpolate_tb(cls, tb1, tb2, alpha):
        """
        Lineraly interpolate between tb1 and tb2 by factor alpha.
        tb = tb1 * (1.0 - alpha) + tb2 * alpha
        """
        import copy

        if tb1.nw != tb2.nw:
            raise ValueError("nw must be identical to interpolate")
        if tb1.nrpts != tb2.nrpts:
            raise ValueError("nrpts must be identical to interpolate")
        if np.any(tb1.rvec != tb2.rvec):
            # TODO: allow change order
            raise ValueError("rvec must be identical to interpolate")
        if np.any(tb1.ndegen != tb2.ndegen):
            # TODO: allow change order
            raise ValueError("ndegen must be identical to interpolate")

        hr_interpol = tb1.hr * (1.0 - alpha) + tb2.hr * alpha
        return cls(nw=tb1.nw, nrpts=tb1.nrpts, rvec=tb1.rvec,
                      ndegen=tb1.ndegen, hr=hr_interpol)

    def write_pkl(self, filename):
        """Write the TBdict object as file"""
        with open(filename+'.pkl', 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def get_hk(self, kvec):
        """Fourier transform TBdict from real space to k space to get H(k)"""
        hk = np.zeros((self.nw, self.nw), dtype=complex)
        for ir in range(self.nrpts):
            rdotk = np.inner(self.rvec[:,ir], kvec)
            hk += self.hr[ir,:,:] * np.exp(1j * 2 * np.pi * rdotk) / self.ndegen[ir]

        herm_check = np.linalg.norm(hk - hk.conjugate().T)
        if herm_check > 1E-6:
            print("WARNING: Hermiticity error: ", herm_check)

        return hk

    def get_ir_ind(self, rvec_target):
        """Find the ir index with rvec = rvec_target"""
        for ir in range(self.nrpts):
            if np.all(self.rvec[:, ir] == rvec_target):
                return ir
        return None

    def subtract_fermi(self, efermi):
        """
        Subtract fermi energy from the diagonal (onsite) tight-binding
        matrix elements
        """
        found = False
        for ir in range(self.nrpts):
            if np.all(self.rvec[:,ir] == [0, 0, 0]):
                found = True
                break
        if not found:
            raise ValueError("rvec == [0, 0, 0] is not found")
        self.hr[ir,:,:] -= np.eye(self.nw) * efermi

    def get_bands(self, kvecs, verbose=False):
        """
        For given array of k points, get H(k) and calculate the eigenvalues
        and return the band structure.

        Args:
            kvecs (3, num_k): k points to calculate band structure
            verbose (logical): If True, print out the progress of k-point loop
        """
        num_k = kvecs.shape[1]
        energy = np.zeros((self.nw, num_k))

        for ik in range(num_k):
            if verbose:
                if ik % 10 == 0: print(f"ik = {ik:3d}", flush=True)

            ham = self.get_hk(kvecs[:,ik])
            eigval, eigvec = np.linalg.eigh(ham)
            energy[:,ik] = eigval

        return energy

    def write_hr_dat(self, filename):
        """Write _hr.dat file in the format of wannier90.x output."""
        import datetime

        with open(filename, 'w') as f:
            f.write(f'Written by wannier_tb.py of Jae-Mo Lihm, {datetime.datetime.now()}\n')
            f.write(f'{self.nw:15d}\n')
            f.write(f'{self.nrpts:15d}\n')

            for ir in range(self.nrpts):
                f.write(f'{self.ndegen[ir]:5d}')
                if (ir+1) % 15 == 0:
                    f.write('\n')
            if self.nrpts % 15 != 0:
                f.write('\n')

            for ir in range(self.nrpts):
                for iw in range(self.nw):
                    for jw in range(self.nw):
                        f.write(f'{self.rvec[0,ir]:5d}{self.rvec[1,ir]:5d}'
                                f'{self.rvec[2,ir]:5d}{jw+1:5d}{iw+1:5d}'
                                f'{self.hr[ir,jw,iw].real:12.6f}{self.hr[ir,jw,iw].imag:12.6f}\n')

# Example usage
# def grep_fermi(filename):
#     with open(filename, 'r') as f:
#         for line in f:
#             if 'highest occupied' in line:
#                 return (float(line.split()[-2]) + float(line.split()[-1])) / 2.0

# nw = 18
# prefix = 'BiTeI'

# tbs = {}
# scalelist = [0.95, 1.00, 1.05]
# for scale in scalelist:
#     folder = f'scale_{scale:.2f}/'

#     pkl_filename = folder + prefix + "_hr"
#     try:
#         tb = TBdict.read_hr_from_pkl(pkl_filename)
#     except FileNotFoundError:
#         tb = TBdict.read_hr_from_bin(folder + prefix)
#         # tb.write_pkl(pkl_filename)

#     efermi = grep_fermi(folder + prefix + '.scf.out')
#     tb.subtract_fermi(efermi)
#     tbs[f'{scale:.2f}'] = tb


# nkdiv = 50
# factor = 0.2
# kvecs_x = np.concatenate((np.linspace(1/3*factor, 0.0, nkdiv, False), np.linspace(0.0, 0.5*factor, nkdiv, True)))
# kvecs_y = np.concatenate((np.linspace(1/3*factor, 0.0, nkdiv, False), np.linspace(0.0, 0.0*factor, nkdiv, True)))
# kvecs_z = np.concatenate((np.linspace(0.5, 0.5, nkdiv, False), np.linspace(0.5, 0.5, nkdiv, True)))
# kvecs = np.vstack((kvecs_x, kvecs_y, kvecs_z))

# scaleplot = np.linspace(0.98, 1.02, 11, True)

# for iscale, scale in enumerate(scaleplot):
#     c = plt.get_cmap('gist_rainbow')(iscale / len(scaleplot))
#     if scale < 1.00:
#         alpha = (1.0 - scale) / 0.05
#         tb = TBdict.interpolate_tb(tbs['1.00'], tbs['0.95'], alpha)
#     elif scale > 1.00:
#         alpha = (scale - 1.0) / 0.05
#         tb = TBdict.interpolate_tb(tbs['1.00'], tbs['1.05'], alpha)

#     energy = tb.get_bands(kvecs)

#     for ib in range(nw):
#         label = f"scale={scale:.4f}" if ib == 0 else None
#         plt.plot(energy[ib, :], c=c, label=label)
# plt.axhline(y=0, c='k', lw=1)
# plt.xlim([0, kvecs.shape[1]])
# plt.ylim([-0.5, 0.5])
# plt.legend()
# plt.show()



