import os
import pdb
import unittest

import scipy as sp
import scipy.linalg as LA

from limix.mtset import MTSet

dir_name = os.path.dirname(os.path.realpath(__file__))
base_folder = os.path.join(dir_name, 'data')


class TestMTSet(unittest.TestCase):
    """test class for optimization"""

    def setUp(self):
        sp.random.seed(0)
        self.Y = sp.loadtxt(os.path.join(base_folder, 'Y.txt'))
        self.XX = sp.loadtxt(os.path.join(base_folder, 'XX.txt'))
        self.Xr = sp.loadtxt(os.path.join(base_folder, 'Xr.txt'))
        self.N, self.P = self.Y.shape
        self.write = False

    def test_mtSetNull(self):
        fbasename = 'mtSetNull'
        setTest = MTSet(Y=self.Y, R=self.XX)
        nullMTInfo = setTest.fitNull(cache=False)
        ext = {'Cg': nullMTInfo['Cg'], 'Cn': nullMTInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSetNull_eigenCache(self):
        fbasename = 'mtSetNull'
        S, U = LA.eigh(self.XX)
        setTest = MTSet(Y=self.Y, S_R=S, U_R=U)
        nullMTInfo = setTest.fitNull(cache=False)
        ext = {'Cg': nullMTInfo['Cg'], 'Cn': nullMTInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSet_base(self):
        fbasename = 'mtSet_base'
        setTest = MTSet(self.Y, R=self.XX)
        optInfo = setTest.optimize(self.Xr)
        ext = {'Cr': optInfo['Cr'], 'Cg': optInfo['Cg'], 'Cn': optInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSet_eigenCache(self):
        fbasename = 'mtSet_base'
        S, U = LA.eigh(self.XX)
        setTest = MTSet(Y=self.Y, S_R=S, U_R=U)
        optInfo = setTest.optimize(self.Xr)
        ext = {'Cr': optInfo['Cr'], 'Cg': optInfo['Cg'], 'Cn': optInfo['Cn']}
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSetPCnull(self):
        fbasename = 'mtSetPCnull'
        setTest = MTSet(Y=self.Y)
        nullMTInfo = setTest.fitNull(cache=False)
        ext = {'Cn': nullMTInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSetPCnull_fixed(self):
        fbasename = 'mtSetPCnull_fixed'
        setTest = MTSet(Y=self.Y, F=self.Xr)
        nullMTInfo = setTest.fitNull(cache=False)
        ext = {'Cg': nullMTInfo['Cg'], 'Cn': nullMTInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSetPC_base(self):
        fbasename = 'mtSetPC_base'
        setTest = MTSet(Y=self.Y)
        optInfo = setTest.optimize(self.Xr)
        ext = {'Cr': optInfo['Cr'], 'Cn': optInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def test_mtSetPC_fixed(self):
        fbasename = 'mtSetPC_fixed'
        setTest = MTSet(Y=self.Y, F=self.Xr[:, :2])
        optInfo = setTest.optimize(self.Xr)
        ext = {'Cr': optInfo['Cr'], 'Cn': optInfo['Cn']}
        if self.write:
            self.saveStuff(fbasename, ext)
        RV = self.assess(fbasename, ext)
        self.assertTrue(RV)

    def saveStuff(self, fbasename, ext):
        """ util function """
        base = os.path.join(base_folder, 'res_' + fbasename + '_')
        for key in list(ext.keys()):
            sp.savetxt(base + key + '.txt', ext[key])

    def loadStuff(self, fbasename, keys):
        """ util function """
        RV = {}
        base = os.path.join(base_folder, 'res_' + fbasename + '_')
        for key in keys:
            RV[key] = sp.loadtxt(base + key + '.txt')
        return RV

    def assess(self, fbasename, ext):
        """ returns a bool vector """
        real = self.loadStuff(fbasename, list(ext.keys()))
        RV = sp.all([((ext[key] - real[key])**2).mean() < 1e-6
                     for key in list(ext.keys())])
        if not RV:
            pdb.set_trace()
        return RV


if __name__ == '__main__':
    unittest.main()
