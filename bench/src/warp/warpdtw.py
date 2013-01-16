#############################################################################
# Dynamic time warping

from imports import *
from warp import *

#############################################################################

#pngDir = "./png/
pngDir = None
seed = 31
nrms = 0.0

def main(args):
  #goSequences()
  #goShifts()
  #goSinoImages()
  goSinoShifts()
  #goCubic()

def goCubic():
  xk = rampfloat(0.0,0.2,26)
  uk = [0.000,  0.140,  0.496,  0.628,  0.752, 
        0.860,  0.968,  1.076,  1.172,  1.272, 
        1.352,  1.452,  1.512,  1.580,  1.716, 
        1.844,  1.844,  2.052,  2.112,  2.256, 
        2.436,  2.436,  2.664,  2.736,  2.832, 
        3.000]
  ci = CubicInterpolator(CubicInterpolator.Method.MONOTONIC,xk,uk)
  n = 1251
  x = rampfloat(0.0,0.004,n)
  u = zerofloat(n)
  ci.interpolate(x,u)
  SimplePlot.asPoints(x,u)

def goSinoImages():
  s1,s2,f,g = getSinoImages()
  SimplePlot.asPixels(s1,s2,f)
  SimplePlot.asPixels(s1,s2,g)

def goSinoShifts():
  s1,s2,f,g = getSinoImages()
  print "d1 =",s1.delta," d2 =",s2.delta
  sf = Sampling(1+int((5.0-s1.first)/s1.delta),s1.delta,0.0)
  sg = s1
  s1 = sf
  n1,n2 = s1.count,s2.count
  f = copy(n1,n2,f)
  smin,smax = 0.0,sg.last-sf.last
  dtw = DynamicWarpingR(smin,smax,s1,s2)
  ss = dtw.getSamplingS()
  ns = ss.count
  e = zerofloat(ns,n1)
  firstTime = False
  if firstTime:
    for i2 in range(n2): #range(n2/2-11,n2/2+12):
      ei = dtw.computeErrors(sf,f[i2],sg,g[i2])
      add(ei,e,e)
    dtw.normalizeErrors(e)
    aos = ArrayOutputStream("/data/seis/sino/esino.dat")
    aos.writeFloats(e)
    aos.close()
  else:
    ais = ArrayInputStream("/data/seis/sino/esino.dat")
    ais.readFloats(e)
    ais.close()
  plotc(s1,ss,e,       cmax=0.5)
  dtw.setStrainLimits(0.1,1.6,-0.05,0.05) # 25 ms/km
  dtw.setSmoothness(0.60,2.0)
  u = dtw.findShifts(e)
  plotc(s1,ss,e,None,u,cmax=0.5)
  vpvs = gamma1(s1,u)
  SimplePlot.asPoints(s1,vpvs)
  return
  u = dtw.findShifts(sf,f,sg,g)
  plotc(s1,ss,e,None,u[n2/2],cmax=0.5)
  h = zerofloat(n1,n2)
  for i2 in range(n2):
    h[i2] = dtw.applyShifts(sg,g[i2],u[i2])
  vpvs = gamma2(s1,u)
  def plotf(s1,s2,f,cmin=0.0,cmax=0.0,cmap=gray):
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    pv = sp.addPixels(s1,s2,f)
    if cmin<cmax:
      pv.setClips(cmin,cmax)
    pv.setColorModel(cmap)
    if (s1.last<6.0):
      sp.setVLimits(0.5,2.5)
    else:
      sp.setVLimits(0.8,4.0)
    sp.addColorBar();
    sp.setSize(700,700)
  plotf(s1,s2,u,cmin=0.5,cmax=1.5,cmap=prism)
  plotf(s1,s2,vpvs,cmin=1.4,cmax=2.2,cmap=jet)
  plotf(sf,s2,f)
  plotf(sf,s2,h)
  plotf(sg,s2,g)

def gamma1(st,u):
  n1 = len(u)
  dudt = zerofloat(n1)
  dt = st.delta
  h = zerofloat(3)
  h[0] =  0.5/dt
  h[1] =  0.0
  h[2] = -0.5/dt
  Conv.conv(3,-1,h,n1,0,u,n1,0,dudt)
  dudt[0] = (u[1]-u[0])/dt
  dudt[n1-1] = (u[n1-1]-u[n1-2])/dt
  vpvs = add(1.0,mul(2.0,dudt))
  return vpvs

def gamma2(st,u):
  n1,n2 = len(u[0]),len(u)
  dudt = zerofloat(n1,n2)
  dt = st.delta
  h = zerofloat(3)
  h[0] =  0.5/dt
  h[1] =  0.0
  h[2] = -0.5/dt
  for i2 in range(n2):
    Conv.conv(3,-1,h,n1,0,u[i2],n1,0,dudt[i2])
    dudt[i2][0] = (u[i2][1]-u[i2][0])/dt
    dudt[i2][n1-1] = (u[i2][n1-1]-u[i2][n1-2])/dt
  vpvs = add(1.0,mul(2.0,dudt))
  return vpvs

def goShifts():
  sf,f,sg,g,st,s = makeSequences()
  dt = st.delta
  s0min,s0max = minmax0(st,s)
  s1min,s1max = minmax1(st,s)
  s2min,s2max = minmax2(st,s)
  u0,rmin,rmax = 0.0,0.0,2.0*s1max
  u2max = 2.0*max(abs(s2min),abs(s2max))
  nr = 1+int((rmax-rmin)/(u2max*dt))
  nr = max(nr,2)
  print "nr =",nr
  umin,umax = s0min,s0max+20*dt
  dtw = DynamicTimeWarping(st,umin,umax)
  ss = dtw.getShiftSampling()
  e = dtw.computeErrors(sf,f,sg,g)
  u = dtw.findShifts(0.0,nr,rmin,rmax,e)
  printShiftStats("s",st,s)
  printShiftStats("u",st,u)
  plotc(st,ss,e,s,u,cmin=0,cmax=1.6)
  plotc(st,ss,e,s,  cmin=0,cmax=1.6)
  plotc(st,ss,e,    cmin=0,cmax=1.6)
  plotfg(sf,f,sg,g)

def makeSequences():
  nt,dt,ft = 501,0.004,0.000
  st = Sampling(nt,dt,ft)
  a,b = 0.8,1.5*2.0*PI/(st.last-st.first)
  u = makeShiftsSine(st,a,b)
  #u = makeShiftsRandom(st)
  umin,umax = min(u),max(u)
  ng,dg,fg = int(nt+umax/dt),dt,ft
  sg = Sampling(ng,dg,fg)
  fpeak = 0.5/8 # 1/8 of Nyquist
  dtw = DynamicTimeWarping(st,umin,umax)
  g = makeRandomEvents(ng,seed=seed); 
  f = dtw.applyShifts(sg,g,st,u)
  f = addRickerWavelet(fpeak,f)
  g = addRickerWavelet(fpeak,g)
  f = addNoise(nrms,fpeak,f,seed=seed+1)
  g = addNoise(nrms,fpeak,g,seed=seed+2)
  return st,f,sg,g,st,u

def minmax0(st,u):
  return min(u),max(u)
def minmax1(st,u):
  nt = st.count
  dt1 = st.delta
  u1 = zerofloat(nt-2)
  Conv.conv(3,-1,(0.5/dt1, 0.0/dt1,-0.5/dt1),nt,0,u,nt-2,1,u1);
  return min(u1),max(u1)
def minmax2(st,u):
  nt = st.count
  dt1 = st.delta
  dt2 = dt1*dt1
  u2 = zerofloat(nt-2)
  Conv.conv(3,-1,(1.0/dt2,-2.0/dt2, 1.0/dt2),nt,0,u,nt-2,1,u2);
  return min(u2),max(u2)
def printShiftStats(pre,st,u):
  u0min,u0max = minmax0(st,u)
  u1min,u1max = minmax1(st,u)
  u2min,u2max = minmax2(st,u)
  print pre+"0: min =",u0min," max =",u0max
  print pre+"1: min =",u1min," max =",u1max
  print pre+"2: min =",u2min," max =",u2max

def goSequences():
  sf,f,sg,g,st,u = makeSequences()
  plotfg(sf,f,sg,g)

def makeShiftsSine(st,a,b):
  nt,dt,ft = st.count,st.delta,st.first
  t = rampfloat(ft,dt,nt)
  u = mul(0.5,add(t,mul(a/b,sin(mul(b,t)))))
  return u

def makeShiftsRandom(st):
  nt,dt,ft = st.count,st.delta,st.first
  r = randfloat(Random(seed),nt)
  u = zerofloat(nt)
  for i in range(1,nt):
    u[i] = u[i-1]+r[i]*dt
  return u

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),11.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

def plotShifts(e,u=None):
  nl,n1 = len(e[0]),len(e)
  sp = SimplePlot()
  sp.setSize(1400,500)
  sp.setHLimits(0,n1-1)
  pv = sp.addPixels(s1,sl,transpose(e))
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setColorModel(ColorMap.JET)
  pv.setPercentiles(2,98)
  if u:
    pv = sp.addPoints(s1,u)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(3)

def getSinoImages():
  dataDir = "/data/seis/sino/"
  n1,d1,f1 = 2001,0.004,0.0
  n2,d2,f2 =  721,0.015,0.0
  s1 = Sampling(n1,d1,f1)
  s2 = Sampling(n2,d2,f2)
  f = readImage(dataDir+"z260.dat",n1,n2)
  g = readImage(dataDir+"x260.dat",n1,n2)
  gain(100,f)
  gain(100,g)
  return s1,s2,f,g

#############################################################################
# utilities

def lowpass(f3db,f):
  """ low-pass filter with specified 3dB frequency, in cycles per sample """
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  """ stretch (supersample) by specified factor c time sampling of image f """
  n1,n2 = len(f[0]),len(f)
  si = SincInterp()
  g = zerofloat(n1)
  for i2 in range(n2):
    si.interpolate(n1,1.0,0.0,f[i2],n1,1.0/c,0.0,g)
    copy(g,f[i2])

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def noiseImage(n1,n2):
  r = Random(3)
  x = sub(randfloat(r,n1,n2),0.5)
  rgf = RecursiveGaussianFilter(2.0)
  for x2 in x:
    rgf.apply1(x2,x2)
  return x
 
#############################################################################
# plotting

gray = ColorMap.GRAY
jet = ColorMap.JET
prism = ColorMap.PRISM

def plotfg(sf,f,sg,g,png=None):
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setHeightElastic(1,25)
  fv = panel.addPoints(0,0,sf,f)
  gv = panel.addPoints(1,0,sg,g)
  fv.setLineWidth(2)
  gv.setLineWidth(2)
  panel.setHLimits(0,sg.first,sg.last)
  panel.setVLimits(0,-1.8,1.8)
  panel.setVLimits(1,-1.8,1.8)
  panel.setHLabel("Time (s)")
  panel.setVLabel(0,"PP")
  panel.setVLabel(1,"PS")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1000,470)
  frame.setVisible(True)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotc(st,ss,c,s=None,u=None,cmin=0.0,cmax=0.0,perc=None,png=None):
  panel = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.setHLabel("Time (s)")
  panel.setVLabel("Time shift (s)")
  cv = panel.addPixels(0,0,st,ss,transpose(c))
  #cv.setColorModel(ColorMap.JET)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc:
    cv.setPercentiles(100-perc,perc)
  elif cmin<cmax:
    cv.setClips(cmin,cmax)
  if s:
    sv = panel.addPoints(0,0,st,s)
    sv.setLineColor(Color.RED)
    sv.setLineStyle(PointsView.Line.DOT)
    sv.setLineWidth(3)
  if u:
    uv = panel.addPoints(0,0,st,u)
    uv.setLineColor(Color.RED)
    uv.setLineWidth(3)
  panel.addColorBar("Amplitude")
  frame = PlotFrame(panel)
  frame.setVisible(True)
  panel.setHLimits(0,st.first,st.last)
  panel.setVLimits(0,ss.first,ss.last)
  #panel.setHLimits(0,0.0,3.5)
  #panel.setVLimits(0,0.0,2.0)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setBackground(backgroundColor)
  frame.setFontSizeForPrint(8,240)
  frame.setSize(1430,760)
  if png and pngDir:
    png += "n"+str(int(10*nrms))
    png += "s"+str(int(10*strainMax))
    frame.paintToPng(720,3.33333,pngDir+"/"+png+".png")

def plotp(f,s1,clips=None,limits=None,title=None,cmap=None,cbar=None,png=None):
  n1,n2 = len(f[0]),len(f)
  #width,height,cbwm = 610,815,145
  width,height,cbwm = 550,325,70
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #sp.addGrid("H-").setColor(Color.BLACK)
  if clips:
    pv.setClips(clips[0],clips[1])
  if limits:
    sp.setVLimits(limits[0],limits[1])
  #if cmap:
  #  pv.setColorModel(cmap)
  if cbar:
    cone = cbar=="Amplitude"
    cbar = sp.addColorBar(cbar)
    if cone:
      cbar.setInterval(1)
  sp.setVInterval(0.2)
  if s1==s1f:
    sp.setVLabel("PP time (s)")
  else:
    sp.setVLabel("PS time (s)")
  sp.setHLabel("Distance (km)")
  sp.setFontSizeForPrint(8,240)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,3.33333,pngDir+png+".png")

def plots(f,s1,clips=None,limits=None,title=None,cmap=None,cbar=None,png=None):
  n1,n2 = len(f[0]),len(f)
  #width,height,cbwm = 610,815,145
  width,height,cbwm = 900,900,180
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.plotPanel.setColorBarWidthMinimum(cbwm)
  pv = sp.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #sp.addGrid("H-").setColor(Color.YELLOW)
  if clips:
    pv.setClips(clips[0],clips[1])
  if limits:
    sp.setVLimits(limits[0],limits[1])
  if title:
    sp.setTitle(title)
  if cmap:
    pv.setColorModel(cmap)
  if cbar:
    cone = cbar=="Amplitude"
    cbar = sp.addColorBar(cbar)
    if cone:
      cbar.setInterval(1)
  sp.setVInterval(1.0)
  if s1==s1f:
    sp.setVLabel("Z time (s)")
  else:
    sp.setVLabel("X time (s)")
  sp.setHLabel("Distance (km)")
  sp.setFontSizeForPrint(8,150)
  sp.setSize(width,height)
  sp.setVisible(True)
  if png and pngDir:
    sp.paintToPng(720,2.0,pngDir+png+".png")

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
