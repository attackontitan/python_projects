
import BS
import unittest

class mytest( unittest.TestCase):
    
    def setup(self):
        pass
    
    def tearDown(self):
        pass
    
    def testsum(self):
        if BS.callput == "Call":
            self.assertEqual(BS.bsformula("Call", 100.0, 120.0, 0.0, 1.0, 0.5, 0.), (13.108976212521462, 0.45436395362548687, 39.63292178326284))
        elif BS.callput == "Put":
            self.assertEqual(BS.bsformula("Put", 100.0, 120, 0.0, 1.0, 0.5, 0.), (33.10897621252146, -0.5456360463745131, 39.63292178326284))
        
        
if __name__=='__main__':    
    unittest.main()