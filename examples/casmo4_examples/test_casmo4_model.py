import unittest

from casmo4_model import CASMO4

class TestCASMO4(unittest.TestCase):
    
    def setUp(self):
        self.c = CASMO4(900.0, 580.0, 900.0, 4.0, degree=1, run=True)
    
    def test_rho(self):
        self.assertAlmostEqual(self.c.rho(None, 0, 900, 580, 900), 
                                0.19365177897136179)
        self.assertAlmostEqual(self.c.rho(None, 60, 900, 580, 900), 
                                -0.23938653389214884) 
        self.assertAlmostEqual(self.c.rho(None, 30, 1000, 580, 900), 
                                -0.025583634185665657)     
        self.assertAlmostEqual(self.c.rho(None, 30, 1000, 580, 1000), 
                                -0.033270696474043807)          
        self.assertAlmostEqual(self.c.rho(None, 30, 1000, 550, 1000), 
                                -0.023438146498884996) 
    
    def test_m2(self):
        self.assertAlmostEqual(self.c.m2(None, 0, 900, 580, 900), 
                               60.925593196101815)
        self.assertAlmostEqual(self.c.m2(None, 0, 900, 600, 900), 
                               66.815075565753588) 
        self.assertAlmostEqual(self.c.m2(None, 30, 1000, 600, 800), 
                               64.73392153778633)     
        
if __name__ == '__main__':
    
    unittest.main()