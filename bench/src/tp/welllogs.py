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

# Directories and files for well logs, headers, directional surveys
#wlDir = "/data/seis/tp/WellLogs/LAS_log_files/"
wlDir = "/data/seis/tp/WellLogs/LAS_log_files/Deeper_LAS_files/"
#wlDir = "/data/seis/tp/WellLogs/LAS_log_files/test/"
whFile = "/data/seis/tp/WellLogs/WellHeaders.txt"
dsFile = "/data/seis/tp/WellLogs/DirectionalSurveys.txt"
wlFile = "/data/seis/tp/resamp/tp3logs.dat"

def main(args):
  #loadWellLogs()
  #loadDirectionalSurveys()
  #loadWellHeaders()
  viewCoordinates()

def loadWellLogs():
  wldata = WellLog.Data(wlDir,whFile,dsFile)
  wldata.printInfo()
  wldata.writeBinary(wlFile)

def loadWellHeaders():
  whdata = WellHeader.Data(whFile)
  whdata.printInfo()

def loadDirectionalSurveys():
  dsdata = DirectionalSurvey.Data(dsFile)
  dsdata.printInfo()

# A few test wells.
names = {
  4902510993:"33-X-23",
  4902510925:"56-TpX-3",
  4902510902:"62-TpX-11"}

def viewCoordinates():
  wldata = WellLog.Data()
  wldata.readBinary(wlFile)
  wldata.printInfo()
  spz1 = SimplePlot()
  sp23 = SimplePlot()
  zlist = []
  x1list = []
  x2list = []
  x3list = []
  for log in wldata.getAll():
    zlist.append(log.z)
    x1list.append(log.dx1dz())
    x2list.append(log.x2)
    x3list.append(log.x3)
  tv = PointsView(zlist,x1list)
  spz1.add(tv)
  tv = PointsView(x2list,x3list)
  sp23.add(tv)

def dx1(log):
  return Array.max(log.x1)-Array.min(log.x1)

def dx2(log):
  return Array.max(log.x2)-Array.min(log.x2)
  #return log.x2[0]-log.x2[log.n-1]

def dx3(log):
  return Array.max(log.x3)-Array.min(log.x3)
  #return log.x3[0]-log.x3[log.n-1]

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
