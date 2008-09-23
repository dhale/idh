#############################################################################
# Test Sibson interpolation

import sys
from org.python.util import PythonObjectInputStream
from math import *
from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mesh import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

from fmm import *

#############################################################################
# global parameters

#pngDir = "./png"
pngDir = None

n1,n2 = 251,357 # image dimensions (as for Teapot dome image)
ms = 50 # nominal spacing between nodes

#############################################################################
# tests

def main(args):
  #doShowSibson("u")
  #doShowSibson("p")
  doInterpolateTest()

def doInterpolateTest():
  s1,s2 = Sampling(n1),Sampling(n2)
  #x1,x2 = makeUniformSubsampling(s1,s2,ms)
  #x1,x2 = makePerturbedSubsampling(s1,s2,ms)
  x1,x2 = makeRandomSubsampling(s1,s2,ms)

  # make mesh
  mesh = makeTriMesh(x1,x2,sinusoid)
  fs = sampleFunction(sinusoid,s1,s2)
  plotfab(fs,mesh,None,False,False,"fs")

  # interpolate with natural neighbors and plot
  fi0 = interpolateNatural(mesh,s1,s2)
  plotfab(fi0,mesh,None,False,False,"fi0")

  # interpolate with natural neighbors and plot
  fi1 = interpolateNatural1(mesh,s1,s2)
  plotfab(fi1,mesh,None,False,False,"fi1")
  return

  # interpolate with approximate natural neighbors and plot
  tsmall = 0.01
  fj = interpolateNaturalApproximate(mesh,s1,s2,tsmall)
  plotfab(fj,mesh,None,False,False,"fj")
  fji = Array.sub(fj,fi)
  print " fji max =",Array.max(Array.abs(fji))
  plotfab(fji,mesh,None,False,False,"fji")

  # interpolate with approximate natural neighbors and plot
  tsmall = 10.0
  fk = interpolateNaturalApproximate(mesh,s1,s2,tsmall)
  plotfab(fk,mesh,None,False,False,"fk")
  fki = Array.sub(fk,fi)
  print " fki max =",Array.max(Array.abs(fki))
  plotfab(fki,mesh,None,False,False,"fki")

def doShowSibson(sampling):

  # mesh A with either uniform or perturbed sampling
  s1,s2 = Sampling(n1),Sampling(n2)
  if sampling=="u":
    x1,x2 = makeUniformSubsampling(s1,s2,ms)
  else:
    x1,x2 = makePerturbedSubsampling(s1,s2,ms)
  mesha = makeTriMesh(x1,x2,sinusoid)

  # mesh B is mesh A with one extra sample
  x1min,x1max = Array.min(x1),Array.max(x1);
  x2min,x2max = Array.min(x2),Array.max(x2);
  #y1,y2 = 0.5*(x1max+x1min),0.5*(x2max+x2min)
  y1,y2 = 0.38*(x1max+x1min),0.52*(x2max+x2min)
  meshb = makeTriMesh(x1,x2,sinusoid)
  node = TriMesh.Node(y1,y2)
  meshb.addNode(node)
  fmap = meshb.getNodePropertyMap("f")
  fmap.put(node,Float(sinusoid(y1,y2)))

  # plot sampled function with nodes and triangles from mesh A
  fs = sampleFunction(sinusoid,s1,s2)
  plotfab(fs,None,None,False,False,sampling+"fs")
  plotfab(fs,mesha,None,False,False,sampling+"fsa")
  plotfab(fs,mesha,None,True,False,sampling+"fsat")

  # interpolate with nearest neighbors using mesh A and plot
  fi = interpolateNearest(mesha,s1,s2)
  plotfab(fi,mesha,None,False,False,sampling+"fca")
  plotfab(fi,mesha,None,False,True,sampling+"fcap")
  plotfab(fi,mesha,meshb,False,True,sampling+"fcabp")

  # interpolate with natural neighbors using mesh A and plot
  fi = interpolateNatural(mesha,s1,s2)
  plotfab(fi,mesha,None,False,False,sampling+"fia")
  plotfab(fi,mesha,None,True,False,sampling+"fiat")
  plotfab(fi,mesha,None,False,True,sampling+"fiap")
  plotfab(fi,meshb,None,False,True,sampling+"fibp")
  plotfab(fi,mesha,meshb,False,True,sampling+"fiabp")

def roundToNearestSample(s,x):
  n = len(x)
  for i in range(n):
    x[i] = s.valueOfNearest(x)

def quadratic(x1,x2):
  return 1.0+x1+x2+x1*x1+x1*x2+x2*x2

def sinusoid(x1,x2):
  return sin(0.02*x1)*sin(0.02*x2)

def makePerturbedSubsampling(s1,s2,ms):
  x1,x2 = makeUniformSubsampling(s1,s2,ms)
  n1,n2 = len(x1[0]),len(x1)
  d1,d2 = x1[0][1]-x1[0][0],x2[1][0]-x2[0][0]
  r = Random(314159)
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      dx1 = d1*(r.nextFloat()-0.5)
      dx2 = d2*(r.nextFloat()-0.5)
      x1[i2][i1] += dx1
      x2[i2][i1] += dx2
      x1[i2][i1] = s1.valueOfNearest(x1[i2][i1])
      x2[i2][i1] = s2.valueOfNearest(x2[i2][i1])
  return x1,x2

def makeRandomSubsampling(s1,s2,ms):
  x1,x2 = makeUniformSubsampling(s1,s2,ms)
  n1,n2 = len(x1[0]),len(x1)
  d1,d2 = x1[0][1]-x1[0][0],x2[1][0]-x2[0][0]
  f1,f2 = s1.first,s2.first
  l1,l2 = f1+d1*(n1-1),f2+d2*(n2-1)
  r = Random(314159)
  for i2 in range(1,n2-1):
    for i1 in range(1,n1-1):
      x1[i2][i1] = f1+(l1-f1)*r.nextFloat()
      x2[i2][i1] = f2+(l2-f2)*r.nextFloat()
      x1[i2][i1] = s1.valueOfNearest(x1[i2][i1])
      x2[i2][i1] = s2.valueOfNearest(x2[i2][i1])
  return x1,x2

def makeUniformSubsampling(s1,s2,ms):
  n1,n2 = s1.count,s2.count 
  d1,d2 = s1.delta,s2.delta 
  f1,f2 = s1.first,s2.first 
  m1 = (n1-1)/ms+1
  m2 = (n2-1)/ms+1
  d1 = d1*float(n1-1)/float(m1-1)
  d2 = d2*float(n2-1)/float(m2-1)
  n1,n2 = m1,m2
  x1 = Array.zerofloat(n1,n2)
  x2 = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      x1[i2][i1] = f1+i1*d1
      x2[i2][i1] = f2+i2*d2
      x1[i2][i1] = s1.valueOfNearest(x1[i2][i1])
      x2[i2][i1] = s2.valueOfNearest(x2[i2][i1])
  # adjust outermost samples so that we do not interpolate outside the
  # mesh and to reduce the likelihood of skinny triangles on convex hull
  a1 = 0.0001*d1
  a2 = 0.0001*d2
  b1 = 1.0/(d1*(n1-1))
  b2 = 1.0/(d2*(n2-1))
  for i1 in range(n1):
    c2 = a2*sin(pi*b1*(x1[0][i1]-f1));
    x2[   0][i1] -= c2
    x2[n2-1][i1] += c2
  for i2 in range(n2):
    c1 = a1*sin(pi*b2*(x2[i2][0]-f2));
    x1[i2][   0] -= c1
    x1[i2][n1-1] += c1
  return x1,x2

def makeTriMesh(x1,x2,function):
  mesh = TriMesh()
  fmap = mesh.getNodePropertyMap("f")
  x1min,x1max = Array.min(x1),Array.max(x1)
  x2min,x2max = Array.min(x2),Array.max(x2)
  x1min -= 0.15*(x1max-x1min)
  x1max += 0.15*(x1max-x1min)
  x2min -= 0.15*(x2max-x2min)
  x2max += 0.15*(x2max-x2min)
  mesh.setOuterBox(x1min,x2min,x1max,x2max)
  n1,n2 = len(x1[0]),len(x1)
  for i2 in range(n2):
    for i1 in range(n1):
      x1i = x1[i2][i1]
      x2i = x2[i2][i1]
      node = TriMesh.Node(x1i,x2i)
      mesh.addNode(node)
      fki = function(x1i,x2i)
      fmap.put(node,Float(fki))
  return mesh

def sampleFunction(function,s1,s2):
  n1,n2 = s1.count,s2.count 
  fs = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      fs[i2][i1] = function(x1,x2)
  return fs

def interpolateNatural(mesh,s1,s2):
  fmap = mesh.getNodePropertyMap("f")
  n1,n2 = s1.count,s2.count 
  fi = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      fi[i2][i1] = mesh.interpolateSibson(x1,x2,fmap,0.0)
  return fi

def interpolateNatural1(mesh,s1,s2):
  fmap = mesh.getNodePropertyMap("f")
  fgmap = mesh.getNodePropertyMap("fg")
  mesh.estimateGradients(fmap,fgmap);
  n1,n2 = s1.count,s2.count 
  fi = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      fi[i2][i1] = mesh.interpolateSibson1(x1,x2,fgmap,0.0)
  return fi

def interpolateNearest(mesh,s1,s2):
  fmap = mesh.getNodePropertyMap("f")
  n1,n2 = s1.count,s2.count 
  fi = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      node = mesh.findNodeNearest(x1,x2)
      fi[i2][i1] = fmap.get(node)
  return fi

def interpolateNaturalApproximate(mesh,s1,s2,tsmall):
  meshc = copyMesh(mesh)
  fmap = meshc.getNodePropertyMap("f")
  n1,n2 = s1.count,s2.count 
  fi = Array.zerofloat(n1,n2)
  tmap = distanceMap(meshc,s1,s2)
  heap = TimeHeap2(TimeHeap2.Type.MAX,n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      heap.insert(i1,i2,tmap[i2][i1])
  plotfab(tmap,mesh,None,False,False,"tmapk")
  while not heap.isEmpty():
    e = heap.remove()
    i1,i2,ti = e.i1,e.i2,e.t
    x1 = s1.getValue(i1)
    x2 = s2.getValue(i2)
    fi[i2][i1] = meshc.interpolateSibson(x1,x2,fmap,0.0)
    if ti<tsmall:
      continue
    tmap[i2][i1] = 0.0
    node = TriMesh.Node(x1,x2)
    meshc.addNode(node)
    fmap.put(node,Float(fi[i2][i1]))
    it = 1+int(0.5*ti) # assume sampling intervals both equal one
    j2lo = max(0,i2-it)
    j2hi = min(n2,i2+it+1)
    j1lo = max(0,i1-it)
    j1hi = min(n1,i1+it+1)
    tiny = 0.0001
    for j2 in range(j2lo,j2hi):
      x2 = s2.getValue(j2)
      for j1 in range(j1lo,j1hi):
        x1 = s1.getValue(j1)
        node = meshc.findNodeNearest(x1,x2)
        tj = distance(x1,x2,node.x(),node.y())
        if tj<tmap[j2][j1]-tiny:
          tmap[j2][j1] = tj
          heap.reduce(j1,j2,tj)
  plotfab(tmap,mesh,None,False,False,"tmaps")
  return fi

def distanceMap(mesh,s1,s2):
  n1,n2 = s1.count,s2.count 
  d = Array.zerofloat(n1,n2)
  for i2 in range(n2):
    x2 = s2.getValue(i2)
    for i1 in range(n1):
      x1 = s1.getValue(i1)
      node = mesh.findNodeNearest(x1,x2)
      d[i2][i1] = distance(x1,x2,node.x(),node.y())
  return d

def distance(x1,y1,x2,y2):
  dx = x2-x1
  dy = y2-y1
  return sqrt(dx*dx+dy*dy)

def copyMesh(mesha):
  meshc = TriMesh()
  fmapa = mesha.getNodePropertyMap("f")
  fmapc = meshc.getNodePropertyMap("f")
  nodes = mesha.getNodes()
  while nodes.hasNext():
    nodea = nodes.next()
    nodec = TriMesh.Node(nodea.x(),nodea.y())
    f = fmapa.get(nodea)
    meshc.addNode(nodec)
    fmapc.put(nodec,Float(f))
  return meshc

#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

def plotf(f):
  plotfab(f,None,None)

def plotfa(f,mesh,tris=False,polys=True):
  plotfab(f,mesh,None,tris,polys)

def plotfab(f,mesha=None,meshb=None,tris=False,polys=True,png=None):
  ppo = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  ppa = PlotPanel.AxesPlacement.NONE
  panel = PlotPanel(1,1,ppo,ppa);
  #if png:
  #  panel.setTitle(png)
  mosaic = panel.getMosaic();
  tile = mosaic.getTile(0,0)
  pv = panel.addPixels(f);
  pv.setColorModel(jet)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  tmvo = TriMeshView.Orientation.XDOWN_YRIGHT;
  if meshb:
    tmvb = TriMeshView(meshb);
    tmvb.setOrientation(tmvo);
    tmvb.setLineWidth(1)
    tmvb.setMarkWidth(6)
    tmvb.setTrisVisible(tris);
    tmvb.setPolysVisible(polys);
    tmvb.setMarkColor(Color.WHITE);
    tmvb.setTriColor(Color.WHITE);
    tmvb.setPolyColor(Color.WHITE);
    tile.addTiledView(tmvb);
  if mesha:
    tmva = TriMeshView(mesha);
    tmva.setOrientation(tmvo);
    tmva.setTrisVisible(tris);
    tmva.setPolysVisible(polys);
    tmva.setLineWidth(1)
    tmva.setMarkWidth(6)
    tmva.setMarkColor(Color.BLACK);
    tmva.setTriColor(Color.BLACK);
    tmva.setPolyColor(Color.BLACK);
    tile.addTiledView(tmva);
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  frame.setFontSize(24)
  frame.setSize(1000,720);
  frame.setVisible(True);
  if png and pngDir:
    frame.paintToPng(200,6,pngDir+"/"+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
