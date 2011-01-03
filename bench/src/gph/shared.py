# shared imports and defs

#############################################################################
# imports

from java.awt import *
from java.awt.image import *
from java.lang import *
from java.io import *
from java.nio import *
from java.util import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
#from edu.mines.jtk.la import *
from edu.mines.jtk.lapack import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from gph import *
  
#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
