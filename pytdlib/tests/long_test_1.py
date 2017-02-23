import unittest
import sys

if(sys.argv[1]!="long"):
	sys.exit(77)

sys.argv=sys.argv[:1]

class dummy(unittest.TestCase):
    def dummy(self):
        pass

if __name__ == '__main__':
    unittest.main()
