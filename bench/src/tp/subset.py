"""
Makes directories with subsets of seismic time or depth images.
For example, the directory with name "subz_401_4_600" contains a
subset of a seismic depth image with 400 samples, 4 m sampling 
interval, and the depth of the first sample is 600 m.
Author: Dave Hale, Colorado School of Mines
Version: 2009.07.25
"""
from imports import *

#############################################################################
def main(args):
  #makeSubset("st",st,Sampling(251,0.004,0.500))
  #makeSubset("sz",sz,Sampling(51,0.004,1.400))
  #makeSubset("sz",sz,Sampling(401,0.004,0.400))
  makeSubset("sz",sz,Sampling(401,0.004,0.600))

csmDir = "/data/seis/tp/csm/"
seismictDir = csmDir+"seismict/"
seismiczDir = csmDir+"seismicz/"
st = Sampling(1501,0.002,0.000) # time sampling of input
sz = Sampling(2762,0.002,0.000) # depth sampling of input
s2 = Sampling(357,0.025,0.000) # x2 sampling (both input and output)
s3 = Sampling(161,0.025,0.000) # x3 sampling (both input and output)

def makeSubset(what,s1i,s1o):
  n1i,d1i,f1i = s1i.count,s1i.delta,s1i.first
  n1o,d1o,f1o = s1o.count,s1o.delta,s1o.first
  n1s,d1s,f1s = str(n1o),str(round(d1o*1000)),str(round(f1o*1000))
  k1i = round(d1o/d1i)
  j1i = round((f1o-f1i)/d1i)
  tz = what[1]
  fileName = "tp"+what+".dat"
  inDir = csmDir+"seismic"+tz+"/"
  outDir = inDir+"sub"+tz+"_"+n1s+"_"+d1s+"_"+f1s+"/"
  inFile = inDir+fileName
  outFile = outDir+fileName
  File(outDir).mkdir()
  xi = zerofloat(n1i)
  xo = zerofloat(n1o)
  ais = ArrayInputStream(inFile)
  aos = ArrayOutputStream(outFile)
  for i3 in range(s3.count):
    for i2 in range(s2.count):
      ais.readFloats(xi)
      copy(n1o,j1i,k1i,xi,0,1,xo)
      aos.writeFloats(xo)
  ais.close()
  aos.close()
  makeImageFrame(outFile,s1o,s2,s3)

def makeImageFrame(file,s1,s2,s3):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(file)
  ais.readFloats(x)
  ais.close()
  frame = SimpleFrame.asImagePanels(s1,s2,s3,x)
  frame.setSize(1200,900)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,(n2*d2)/(n1*d1))
  view.setScale(2.0)
  view.setAzimuth(-50.0)
  return frame

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=="__main__":
  SwingUtilities.invokeLater(RunMain()) 
