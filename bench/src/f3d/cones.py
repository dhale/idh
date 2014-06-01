"""
Displays images through conical faults in F3.
Author: Dave Hale, Colorado School of Mines
Version: 2013.03.02
"""
from f3utils import *

#############################################################################

slides = False
pngDir = "png/"
#pngDir = None

# P-wave velocities for sediments with cones ~ 2.2 km/s => dz ~ 1.1*dt 
# For aspect ratio = w:h = 1:1, make pixel_ratio = value_ratio/1.1
# For aspect ratio = w:h = 1:5, make pixel_ratio = value_ratio/5.5

coneSurfBounds = {
  # cone:(x1,x2,x3,h,r)
  0:(35,37,99,80,25),
  1:(14,99,101,108,29),
  2:(7,100,100,78,18),
  3:(10,100,43,69,16),
}
coneViewParams = {
  # cone:(k1,k2,k3,az,el,sc,tx,ty,tz)
  #0:(94,32,121,-34.4,20.4,5.3,0.159,-0.058,0.196),
  0:(94,32,121,-40.5,77.6,5.3,0.111,-0.096,0.255),
  1:(91,84, 85, 41.4,19.0,3.5,0.000,-0.014,0.007),
  2:(81,88,112,-25.7,27.8,3.3,0.000, 0.000,0.000),
  3:(72,89, 54,-54.1,18.7,3.1,0.136, 0.022,0.174),
}

#############################################################################
def main(args):
  #goConeSurfing()
  #goConeSubsets()
  goFigures()
  #makeColorBar3dForSlide(0.0,8.00*1.1*4.0,"Fault throw (m)")

def goConeSurfing():
  makeColorBar3d(0.0,8.00*1.1*4.0,"Fault throw (m)",h=800)
  makeColorBar3d(0.0,8.00*1.1*4.0,"Fault throw (m)",h=400)
  #for cone in [0,1,2,3]:
  for cone in [0]:
    surfCone(cone)
def surfCone(cone):
  scone = str(cone)
  setup("setc"+scone);
  g = readImage("g")
  gs = readImage("gs8")
  p2 = readImage("p2")
  p3 = readImage("p3")
  fl = readImage("fl")
  fp = readImage("fp")
  ft = readImage("ft")
  for csb in [coneSurfBounds[cone]]:
    fs = FaultSurfer3([fl,fp,ft])
    if csb:
      c1,c2,c3,hc,rc = csb
      qf = Util.QuadInsideCone(c1,c2,c3,hc,rc)
      fs.addQuadFilter(qf)
    fs.setThreshold(0.7)
    quads = fs.findQuads()
    for i in range(3): # tend to get one unorientable surface, so
      fs.shuffleQuads(quads) # shuffle to put seams away from cone
    quads = fs.linkQuads(quads)
    surfs = fs.findSurfs(quads)
    surfs = fs.getSurfsWithSize(surfs,2000)
    s = fs.findShifts(20.0,surfs,gs,p2,p3)
    print "s: min =",min(s)," max =",max(s)
    png = "csurf"+str(cone)
    if not csb:
      png += "n"
    plot3(gs,surfs=surfs,smax=8.0,cone=cone,png=png)

def goConeSubsets():
  setup("seta")
  def makeConeSubset(cone):
    x2c,x3c = getConeLocations()
    x2c,x3c = x2c[cone],x3c[cone]
    i2c = s2.indexOfNearest(x2c)
    i3c = s3.indexOfNearest(x3c)
    m1,m2,m3 = 126,201,201
    j1 = s1.indexOfNearest(1.04)
    j2 = max(0,min(n2-m2,i2c-(m2-1)/2))
    j3 = max(0,min(n3-m3,i3c-(m3-1)/2))
    print "cone: "+str(cone)
    print "s1.first = ",s1.getValue(j1)
    print "s2.first = ",s2.getValue(j2)
    print "s3.first = ",s3.getValue(j3)
    if True:
      for name in ["g","gs8","p2","p3","fl","fp","ft","flt","fpt","ftt"]:
        g = readImage(name)
        g = copy(m1,m2,m3,j1,j2,j3,g)
        dataSetDir = f3dDataDir+"setc"+str(cone)+"/"
        aos = ArrayOutputStream(dataSetDir+name+".dat")
        aos.writeFloats(g)
        aos.close()
  for cone in [0,1,2,3]:
    makeConeSubset(cone)
  
def goFigures():
  setup("seta")
  for name in ["g","gs8"]:
    #makeSliceThruPoints(name,"c1")
    #displaySlicesThruCones(name,"c1")
    if not slides:
      displayConesForPrint(name)
    else:
      displayConesForSlides()
 
def displayConesForPrint(name):
  #regional = True # for Rick's interpretation
  regional = False # narrower for composition
  g = readImage(name)
  nt = 126 # sampling of times with cones present
  dt = 0.004
  ft = 1.040
  st = Sampling(nt,dt,ft)
  nd = 73  # sampling of distance d from center of cone
  #nd = 119  # sampling of distance d from center of cone
  dd = s2.delta
  fd = -dd*(nd-1)/2
  sd = Sampling(nd,dd,fd)
  if regional:
    na = 1 # sampling of azimuth a; north = 0 degrees
  else:
    na = 6 # sampling of azimuth a; north = 0 degrees
  da = 30.0
  fa = 0.0
  sa = Sampling(na,da,fa)
  c2s,c3s = getConeLocations()
  if regional:
    ics = [0]
  else:
    ics = [0]
  for ic in ics:
    c2,c3 = c2s[ic],c3s[ic]
    h = makeConeSlices(c2,c3,st,sd,sa,g)
    def plota(ia,aspect11=False):
      sia = str(int(sa.getValue(ia)))
      sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
      #sp.setTitle("Azimuth = "+sia+" degrees")
      sp.setHLabel("Radial distance (km)")
      sp.setVLabel("Time (s)")
      #sp.setHLimits(sd.first,sd.last)
      #sp.setVLimits(1.05,1.55)
      if aspect11:
        if regional:
          sp.setSize(1400,344) # 1:1 scale, width = 3 km
        else:
          sp.setSize(1400,498) # 1:1 scale
      else:
        if regional:
          sp.setSize(498,472) # 1:5 scale, width = 3 km
        else:
          sp.setSize(498,712) # 1:5 scale
      if aspect11:
        wpt = 504.0
      else:
        wpt = 200.0
      sp.setFontSizeForPrint(8,wpt)
      pv = sp.addPixels(st,sd,h[ia])
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv.setClips(-1.0,1.0)
      if not regional:
        pv = sp.addPoints([st.first,st.last],[0,0])
        pv.setLineColor(Color.WHITE)
      pngFile="cone"+str(ic)
      if aspect11:
        pngFile += "r11a"
      else:
        pngFile += "r15a"
      if len(sia)==1:
        sia = "00"+sia
      elif len(sia)==2:
        sia = "0"+sia
      pngFile += sia+name+".png"
      if pngDir:
        sp.paintToPng(360,wpt/72,pngDir+pngFile)
    for ia in range(na):
      #plota(ia,True)
      plota(ia,False)
 
def displayConesForSlides():
  #regional = True # for Rick's interpretation
  regional = False # narrower for composition
  g1 = readImage("g")
  g2 = readImage("gs8")
  nt = 126 # sampling of times with cones present
  dt = 0.004
  ft = 1.040
  st = Sampling(nt,dt,ft)
  nd = 73  # sampling of distance d from center of cone
  #nd = 119  # sampling of distance d from center of cone
  dd = s2.delta
  fd = -dd*(nd-1)/2
  sd = Sampling(nd,dd,fd)
  na = 12 # sampling of azimuth a; north = 0 degrees
  da = 30.0
  fa = 0.0
  sa = Sampling(na,da,fa)
  c2s,c3s = getConeLocations()
  ics = [0,1,2,3]
  for ic in ics:
    c2,c3 = c2s[ic],c3s[ic]
    h1 = makeConeSlices(c2,c3,st,sd,sa,g1)
    h2 = makeConeSlices(c2,c3,st,sd,sa,g2)
    def plota(ia):
      sia = str(int(sa.getValue(ia)))
      pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
      #sp.setTitle("Azimuth = "+sia+" degrees")
      pp.setHLabel(0,"Radial distance (km)")
      pp.setHLabel(1,"Radial distance (km)")
      pp.setVLabel("Time (s)")
      #sp.setHLimits(sd.first,sd.last)
      #sp.setVLimits(1.05,1.55)
      pv1 = pp.addPixels(0,0,st,sd,h1[ia])
      pv1.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv1.setClips(-1.0,1.0)
      pv1 = pp.addPoints(0,0,[st.first,st.last],[0,0])
      pv1.setLineColor(Color.WHITE)
      pv2 = pp.addPixels(0,1,st,sd,h2[ia])
      pv2.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv2.setClips(-1.0,1.0)
      pv2 = pp.addPoints(0,1,[st.first,st.last],[0,0])
      pv2.setLineColor(Color.WHITE)
      pf = PlotFrame(pp)
      pf.setSize(900,718) # 1:5 scale
      pf.setFontSizeForSlide(1.0,0.9)
      pf.setVisible(True)
      pngFile="cone"+str(ic)
      pngFile += "r15a"
      if len(sia)==1:
        sia = "00"+sia
      elif len(sia)==2:
        sia = "0"+sia
      pngFile += sia+"ggs8.png"
      if pngDir:
        pf.paintToPng(256,8,pngDir+pngFile)
    for ia in range(na):
      plota(ia)

def displaySlicesThruCones(name,points):
  ss = getSamplingS(name,points)
  ns,ds,fs = ss.count,ss.delta,ss.first
  fcp = readImage2("slices/"+name+points,n1,ns)
  fk1 = readImage2(getF3dSlice1Name(name,309),n2,n3)
  x2c,x3c = getConeLocations()
  x2p,x3p = getPointSet(points)
  x2s,x3s = getCurveThruPoints(x2p,x3p)
  # Horizontal slice 1.24 seconds (k1 = 309)
  sp = SimplePlot()
  if name[0:2]=="fl": fk1 = neg(fk1)
  pv = sp.addPixels(s2,s3,fk1)
  if name[0:2]=="fl": 
    pv.setClips(-1,0)
  else:
    pv.setClips(-clip,clip)
  pv = sp.addPoints(x2s,x3s)
  if name[0:2]=="fl": 
    pv.setLineColor(Color.BLACK)
  else:
    pv.setLineColor(Color.WHITE)
  pv.setLineWidth(1.0)
  pv = sp.addPoints(x2c,x3c)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.HOLLOW_CIRCLE)
  pv.setMarkSize(48.0)
  if name[0:2]=="fl": 
    pv.setMarkColor(Color.BLACK)
  else:
    pv.setMarkColor(Color.WHITE)
  sp.setHLimits(s2.first,s2.last)
  sp.setVLimits(s3.first,s3.last)
  sp.setHInterval(5.0)
  sp.setVInterval(5.0)
  sp.setHLabel("Inline (km)")
  sp.setVLabel("Crossline (km)")
  if slides:
    sp.setFontSizeForSlide(1.0,0.9)
    sp.setSize(900,680)
    if pngDir:
      sp.paintToPng(256,8.0,pngDir+name+points+"23.png")
  else:
    wpt = 504.0
    sp.setFontSizeForPrint(8,wpt)
    sp.setSize(717,525)
    if pngDir:
      sp.paintToPng(360,wpt/72,pngDir+name+points+"23.png")
  # Vertical slice through cone points
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,ss,fcp)
  pv.setClips(-clip,clip)
  x1,xs = linesThruCones(x2c,x3c,x2s,x3s,ss)
  pv = sp.plotPanel.addPoints(x1,xs)
  pv.setLineColor(Color.WHITE)
  sp.setVLimits(1.0,s1.last)
  sp.setHLimits(ss.first,ss.last)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  if slides:
    sp.setFontSizeForSlide(1.0,0.9)
    sp.setSize(670,465)
    if pngDir:
      sp.paintToPng(256,8.0,pngDir+name+points+"1s.png")
  else:
    wpt = 504.0
    sp.setFontSizeForPrint(8,wpt)
    sp.setSize(670,255)
    if pngDir:
      sp.paintToPng(360,wpt/72,pngDir+name+points+"1s.png")

# Gets cone locations in CSM (x2,x3) coordinates
def getConeLocations():
  x2c = [ 0.910, 2.387, 4.533, 7.002]
  x3c = [11.896, 6.229, 6.140, 1.013]
  return x2c,x3c

def getPointSet(points):
  x2c,x3c = getConeLocations()
  if (points=="c1"):
    x2p = [x2c[0],x2c[1],x2c[3]]
    x3p = [x3c[0],x3c[1],x3c[3]]
  return x2p,x3p

# Gets finely sampled curve through points.
def getCurveThruPoints(x2s,x3s):
  ns = len(x2s)
  for i in range(2): # 2 iterations should be sufficient
    ds = zerofloat(ns)
    ds[0] = 0.0
    for js in range(1,ns):
      ds[js] = ds[js-1]+hypot(x2s[js]-x2s[js-1],x3s[js]-x3s[js-1])
    method = CubicInterpolator.Method.SPLINE
    ci2 = CubicInterpolator(method,ds,x2s)
    ci3 = CubicInterpolator(method,ds,x3s)
    smin,smax = ds[0],ds[-1]
    smin -= 0.500
    smax += 0.500
    ns = 1+int((smax-smin)/s2.delta)
    ds = (smax-smin)/(ns-1)
    sj = rampfloat(smin,ds,ns)
    x2s = zerofloat(ns)
    x3s = zerofloat(ns)
    ci2.interpolate(sj,x2s)
    ci3.interpolate(sj,x3s)
  return x2s,x3s

# Gets 2D seismic image along specified curve.
def getImageAlongCurve(name,x2s,x3s):
  f = readImage(name)
  ns = len(x2s)
  g = zerofloat(n1,ns)
  si = SincInterpolator()
  for js in range(ns):
    x2j = x2s[js]
    x3j = x3s[js]
    for j1 in range(n1):
      x1j = s1.getValue(j1)
      g[js][j1] = si.interpolate(s1,s2,s3,f,x1j,x2j,x3j)
  return f,g

# Writes a 2D image slice that passes through point set
def makeSliceThruPoints(name,points):
  x2p,x3p = getPointSet(points)
  x2s,x3s = getCurveThruPoints(x2p,x3p)
  f,fs = getImageAlongCurve(name,x2s,x3s)
  print "slice through points has",len(fs),"traces"
  writeImage2("slices/"+name+points,fs)

# Sampling of 2nd dimension for slices
def getSamplingS(name,points):
  f = File(getF3dDataSetDir()+"slices/"+name+points+".dat")
  n2 = f.length()/n1/4
  return Sampling(n2,d2,0.0)

def linesThruCones(x2c,x3c,x2s,x3s,ss):
  nc = len(x2c)
  x1,xs = [],[]
  for ic in range(nc):
    i = indexOfNearestPoint(x2c[ic],x3c[ic],x2s,x3s)
    if i>=0:
      xi = ss.getValue(i)
      x1.append([1.04,1.54])
      xs.append([xi,xi])
  return x1,xs

def indexOfNearestPoint(x2,x3,x2s,x3s):
  def distanceSquared(xa,ya,xb,yb):
    dx = xb-xa
    dy = yb-ya
    return dx*dx+dy*dy
  ns = len(x2s)
  jsmin = ns
  dsmin = Double.MAX_VALUE
  for js in range(ns):
    ds = distanceSquared(x2,x3,x2s[js],x3s[js])
    if ds<dsmin:
      dsmin = ds
      jsmin = js
  if dsmin>d2:
    jsmin = -1
  return jsmin

def makeConeSlices(c2,c3,st,sd,sa,g):
  nt,nd,na = st.count,sd.count,sa.count
  dt,dd,da = st.delta,sd.delta,sa.delta
  ft,fd,fa = st.first,sd.first,sa.first
  # output array is a 3D image[na][nd][n1]
  h = zerofloat(nt,nd,na)
  # survey azimuth is 88.4 degrees
  ps = toRadians(88.4)
  # interpolate slices
  si = SincInterpolator()
  for ja in range(na):
    aj = sa.getValue(ja)
    pj = toRadians(aj)
    cj = cos(pj-ps)
    sj = sin(pj-ps)
    for jd in range(nd):
      dj = sd.getValue(jd)
      x2 = c2+cj*dj;
      x3 = c3-sj*dj;
      for jt in range(nt):
        x1 = st.getValue(jt)
        h[ja][jd][jt] = si.interpolate(s1,s2,s3,g,x1,x2,x3)
  return h

def setup(dataset):
  global s1,s2,s3,n1,n2,n3,d1,d2,d3,f1,f2,f3,f3dDataDir,clip
  setupForDataSet(dataset)
  s1,s2,s3 = getSamplings()
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  f3dDataDir = getF3dDataDir()
  clip = 1.0

#############################################################################
# plotting

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.RED_WHITE_BLUE,alpha)
def jetRamp():
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,1.0/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.RED_WHITE_BLUE,a)

def plot3(f,g=None,gmin=None,gmax=None,gmap=None,
          xyz=None,surfs=None,smax=None,cone=None,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  sf = SimpleFrame()
  sf.setBackground(Color(255,255,255))
  if g==None:
    ipg = sf.addImagePanels(f)
    ipg.setClips(-1.0,1.0)
  else:
    ipg = ImagePanelGroup2(f,g)
    ipg.setClips1(-1.0,1.0)
    if gmap==None:
      gmap = jetFill(0.8)
    ipg.setColorModel2(gmap)
    if gmin and gmax:
      ipg.setClips2(gmin,gmax)
    sf.world.addChild(ipg)
  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if surfs:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    isurf = 0
    for surf in surfs:
      #surf.blocky()
      if smax:
        xyz,uvw,rgb = surf.getXyzUvwRgbShifts(smax)
      else:
        xyz,uvw,rgb = surf.getXyzUvwRgb()
      #qg = QuadGroup(False,xyz,rgb)
      qg = QuadGroup(True,xyz,rgb) #qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      #if isurf==11: sg.addChild(qg) # cone0
      #if isurf==6: sg.addChild(qg) # cone0
      #if isurf==5: sg.addChild(qg) # cone1
      #if isurf==6: sg.addChild(qg) # cone2
      #if isurf==3: sg.addChild(qg) # cone3
      sg.addChild(qg)
      print "isurf =",isurf," size =",surf.size()," qg =",qg
      isurf += 1
    sf.world.addChild(sg)
  #ipg.setSlices(209,12,18)
  #ipg.setSlices(200,0,0)
  #ipg.setSlices(80,9,13)
  #ipg.setSlices(80,9,209)
  #sf.setSize(1300,1100)
  sf.setSize(1040,1124)
  sf.setWorldSphere(n3/2,n2/2,n1/2,0.5*sqrt(n1*n1+n2*n2+n3*n3))
  if cone!=None:
    k1,k2,k3,az,el,sc,tx,ty,tz = coneViewParams[cone]
    ipg.setSlices(k1,k2,k3)
    sf.orbitView.setAzimuthAndElevation(az,el)
    sf.orbitView.setScale(sc)
    sf.orbitView.setTranslate(Vector3(tx,ty,tz))
    sf.orbitView.setAxesScale(1.0,1.0,5.0*1.1*0.004/0.025) # 1:5 aspect ratio
  #sf.viewCanvas.setBackground(sf.getBackground())
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")

def makeColorBar3d(cmin,cmax,clab,cint=None,h=800):
  cbar = ColorBar(clab)
  cbar.setFont(Font("Arial",Font.PLAIN,24))
  cbar.setBackground(Color.WHITE)
  if cint:
    cbar.setInterval(cint)
  cbar.setWidthMinimum(80)
  cmap = ColorMap(cmin,cmax,ColorMap.JET)
  cmap.addListener(cbar)
  frame = JFrame()
  frame.setSize(140,h)
  frame.add(cbar,BorderLayout.EAST)
  frame.setVisible(True)
  if pngDir:
    cbar.paintToPng(400,0.5,pngDir+"cbar"+str(h)+".png")

def makeColorBar3dForSlide(cmin,cmax,clab,cint=None):
  cbar = ColorBar(clab)
  cbar.setFont(Font("Arial",Font.PLAIN,40))
  cbar.setBackground(Color.WHITE)
  if cint:
    cbar.setInterval(cint)
  cbar.setWidthMinimum(140)
  cmap = ColorMap(cmin,cmax,ColorMap.JET)
  cmap.addListener(cbar)
  frame = JFrame()
  frame.setSize(170,700)
  frame.add(cbar,BorderLayout.EAST)
  frame.setVisible(True)
  if pngDir:
    cbar.paintToPng(400,0.5,pngDir+"cbar.png")

#############################################################################
run(main)
