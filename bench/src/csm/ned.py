"""
Reads and displays NED data in grid float format.
Author: Dave Hale, Colorado School of Mines
Version: 2012.05.10
"""
from imports import *

#############################################################################
def main(args):
  goPlot()

def goTest():
  e = elevation(37.2694,-107.0092) # Pagosa Springs, Colorado
  print "elevation =",e

def goPlot():
  fileName,s2,s1 = samplingA()
  n1,n2 = s1.count,s2.count
  e = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readFloats(e)
  j2 = n2-1
  for i2 in range(n2/2):
    ei2 = e[i2]
    e[i2] = e[j2]
    e[j2] = ei2
    j2 -= 1
  ais.close()
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  pv = sp.addPixels(s1,s2,e)
  pv.setColorModel(ColorMap.JET)
  sp.setHLabel("Longitude (degrees)")
  sp.setVLabel("Latitude (degrees)")
  sp.addColorBar("Elevation (m)")

def elevation(lat,lon):
  for sampling in [samplingA,samplingB]:
    fileName,slat,slon = sampling()
    e = readElevation(fileName,slat,slon,lat,lon)
    if e: break
  return e

def samplingA():
  fileName = "/data/seis/csm/fc2012/ned/n38w108/floatn38w108_13.flt"
  slat = Sampling(10812,0.0000925925926, 36.99944444444)
  slon = Sampling(10812,0.0000925925926,-108.0005555556)
  return fileName,slat,slon

def samplingB():
  fileName = "/data/seis/csm/fc2012/ned/n38w107/floatn38w107_13.flt"
  slat = Sampling(10812,0.0000925925926, 36.99944444444)
  slon = Sampling(10812,0.0000925925926,-107.0005555556)
  return fileName,slat,slon

def readElevation(fileName,slat,slon,lat,lon):
  ilat = slat.indexOfNearest(lat)
  ilon = slon.indexOfNearest(lon)
  lati = slat.getValue(ilat)
  loni = slon.getValue(ilon)
  nlat = slat.count
  if lati<slat.first or lati>slat.last:
    return None
  if loni<slon.first or loni>slon.last:
    return None
  ais = ArrayInputStream(fileName)
  ais.skipBytes(4*(ilon+(nlat-1-ilat)*slon.count))
  e = ais.readFloat()
  ais.close()
  return e

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
