#############################################################################
# Matern covariance and approximation based on smoothing filters

from shared import *

#pngDir = "./png/"
pngDir = None

def main(args):
  goFigures()
  #goCompare()

def goCompare():
  n1,n2 = 201,201
  d1,d2 = 2.0,2.0
  f1,f2 = 0.0,0.0
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  sigma,shape,scale = 3.0,1.0,20.0*d1
  mc = MaternCovariance(sigma,shape,scale)
  sc = SmoothCovariance(sigma,shape,scale,2)
  j2,j1 = n2/2,n1/2
  p = zerofloat(n1,n2)
  p[j2][j1] = 1.0
  et = makeIdentityTensors(s1,s2)
  sc.testSpd(n1,n2,et)
  sc.apply(s1,s2,et,p)
  sc.applyInverse(s1,s2,et,p)
  pa = p[j2]
  pm = zerofloat(n1)
  ps = zerofloat(n1)
  xj = s1.getValue(j1)
  for i1 in range(n1):
    xi = s1.getValue(i1)
    r = abs(xi-xj)
    pm[i1] = mc.evaluate(r)
    ps[i1] = sc.evaluate(r)
  SimplePlot.asPoints(s1,pm)
  SimplePlot.asPoints(s1,ps)
  SimplePlot.asPoints(s1,pa)

def goFigures():
  ndim = 2
  sigma = 1.0
  scale = 1.0
  rmin,rmax = 0.0,4.0*scale
  nr = 1001
  dr = (rmax-rmin)/(nr-1)
  fr = rmin
  sr = Sampling(nr,dr,fr)
  shapes = [0.5,0.75,1.0,1.5]
  nshape = len(shapes)
  cm = zerofloat(nr,nshape)
  cs = zerofloat(nr,nshape)
  for scaled in [False,True]:
    for ishape in range(nshape):
      shape = shapes[ishape]
      if scaled:
        mc = MaternCovariance(sigma,shape,scale)
        sc = SmoothCovariance(sigma,shape,scale,ndim)
      else:
        mc = MaternCovariance(sigma,shape,scale*2*sqrt(shape))
        sc = SmoothCovariance(sigma,shape,scale*2*sqrt(shape),ndim)
      for ir in range(nr):
        r = sr.getValue(ir)
        cm[ishape][ir] = mc.evaluate(r)
        cs[ishape][ir] = sc.evaluate(r)
    if scaled:
      suffix = "s"
    else:
      suffix = "r"
    plotCovariances(sr,cm,"cm"+suffix)
    plotCovariances(sr,cs,"cs"+suffix)
  
def plotCovariances(sx,ys,png=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  styles = [PointsView.Line.DASH_DOT,
            PointsView.Line.DASH,
            PointsView.Line.SOLID,
            PointsView.Line.DOT]
  colors = [Color.LIGHT_GRAY,Color.GRAY,Color.DARK_GRAY,Color.BLACK]
  for iy in range(len(ys)): 
    pv = sp.addPoints(sx,ys[iy])
    pv.setLineWidth(2)
    #pv.setLineColor(colors[iy%len(colors)])
    pv.setLineStyle(styles[iy%len(styles)])
    sp.setHLabel("distance")
    sp.setVLabel("covariance")
    sp.setHLimits(0.0,4.0);
    sp.setVLimits(0.0,1.03);
    sp.setFontSizeForPrint(8,240.0)
    sp.setSize(790,500);
  if pngDir and png:
    sp.paintToPng(720,240.0/72.0,pngDir+png+".png")

def makeIdentityTensors(s1,s2):
  n1,n2 = s1.count,s2.count;
  au = fillfloat(1.0,n1,n2)
  av = fillfloat(1.0,n1,n2)
  u1 = fillfloat(1.0,n1,n2)
  u2 = zerofloat(n1,n2)
  return EigenTensors2(u1,u2,au,av)

#############################################################################
if __name__ == "__main__":
  run(main)
