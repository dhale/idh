import sys
from math import *

from java.awt import *
from java.io import *
from java.lang import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

from tp import *

dataDir = "/data/seis/tp/WellLogs/LAS_log_files/"

def main(args):
  #loadWellLogs(dataDir)
  loadWellLogs(dataDir+"Deeper_LAS_files")
  #loadWellLogs(dataDir+"Deeper_LAS_files/490252304700")
  #loadWellLog(dataDir+"Shallow_LAS_files/49025103080000_486907.las")

def loadWellLog(fileName):
  WellLog.load(File(fileName))

def loadWellLogs(dirName):
  logs = WellLogs()
  logs.load(dirName)

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
