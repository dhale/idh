from imports import *

pngDir = None
#pngDir = '/Users/sluo/Desktop/'
dataDir = '/data/sluo/pnz/'
n1,n2,n3 = 751,1001,501

#############################################################################

def main(args):
  #crop()
  show()
  #goFlattenS()
  #goFlattenR()

def show():
  f = read('f'); display(f,perc=99.5,name='f')
  gs = read('gs'); display(gs,perc=99.5,name='gs')
  ggs = read('ggs'); display(ggs,perc=99.5,name='ggs')
  #gr = read('gr'); display(gr,perc=99.5,name='gr')

def goFlattenR():
  computeNormals = True
  f = read('f')
  u1 = zeros(n1,n2,n3)
  u2 = zeros(n1,n2,n3)
  u3 = zeros(n1,n2,n3)
  ep = zeros(n1,n2,n3)
  if computeNormals:
    print "estimating normal vectors..."
    lof = LocalOrientFilter(8.0,2.0)
    lof.applyForNormalPlanar(f,u1,u2,u3,ep)
    #display(u1,cmap=jet,name='u1')
    #display(u2,cmap=jet,name='u2')
    #display(u3,cmap=jet,name='u3')
    #display(ep,cmap=jet,name='ep')
    write('u1',u1)
    write('u2',u2)
    write('u3',u3)
    write('ep',ep)
  read('u1',u1)
  read('u2',u2)
  read('u3',u3)
  read('ep',ep); pow(ep,6.0,ep)
  p = array(u1,u2,u3,ep)
  fl = FlattenerRT(8.0,8.0)
  s = fl.findShifts(p)
  g = FlattenerUtil.applyShiftsR(f,s)
  #display(s,cmap=jet,name='s')
  display(g,perc=99.5,name='g')
  display(f,perc=99.5,name='f')
  write('gr',g)
  #write('s',s)

def goFlattenS():
  computeSlopes = True
  #ffile = 'f'
  ffile = 'gs'
  f = read(ffile)
  p2,p3,ep = zeros(n1,n2,n3),zeros(n1,n2,n3),zeros(n1,n2,n3)
  """
  if computeSlopes:
    print "estimating slopes..."
    lsf = LocalSlopeFinder(8.0,2.0,5.0)
    lsf.findSlopes(f,p2,p3,ep)
    #display(p2,cmap=jet,name='p2')
    #display(p3,cmap=jet,name='p3')
    #display(ep,cmap=jet,name='ep')
    write('p2',p2)
    write('p3',p3)
    write('ep',ep)
  read('p2',p2)
  read('p3',p3)
  read('ep',ep); pow(ep,6.0,ep)
  """
  lsf = LocalSlopeFinder(8.0,2.0,5.0)
  lsf.findSlopes(f,p2,p3,ep)
  pow(ep,6.0,ep)
  fl = FlattenerS(8.0,8.0)
  s = fl.findShifts(p2,p3,ep)
  g = fl.applyShifts(f,s)
  display(s,cmap=jet,name='s')
  display(g,perc=99.5,name='g')
  display(f,perc=99.5,name='f')
  if ffile=='f':
    gfile = 'gs'
  elif ffile=='gs':
    gfile = 'ggs'
  write(gfile,g)
  #write('s',s)

def crop():
  f = zeros(751,1001,1001)
  read('pnz10',f,'/data/seis/pnz/dat/')
  g = copy(751,1001,501,0,0,500,f)
  #display(f,perc=99.5)
  display(g,perc=99.5)
  #write('f',g)

def read(name,image=None,dir=None):
  if not image:
    image = zeros(n1,n2,n3)
  if not dir:
    dir = dataDir
  fileName = dir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def write(name,image,dir=None):
  if not dir:
    dir = dataDir
  fileName = dir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

def zeros(n1,n2,n3):
  return zerofloat(n1,n2,n3)

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))

#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=gray,cmin=0,cmax=0,perc=100,cbar=None,name=None):
  pan = panel(cbar)
  pix = pan.addPixels(x)
  pix.setColorModel(cmap)
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  pix.setInterpolation(PixelsView.Interpolation.LINEAR)
  frame(pan,name)

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  #p.setHLabel("index i2")
  #p.setVLabel("index i1")
  cb = p.addColorBar()
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name=None):
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(800,600)
  if name:
    frame.setTitle(name)
  frame.setVisible(True)
  if name and pngDir:
    frame.paintToPng(360,3.0,pngDir+name+'.png')

#############################################################################
# graphics

def display(image,cmap=gray,cmin=0,cmax=0,perc=100,name=None):
  world = World()
  addImageToWorld(world,image,cmap,cmin,cmax,perc)
  makeFrame(world,name)

def display2(image1,image2,cmap1=gray,cmap2=gray,name=None):
  world = World()
  addImageToWorld(world,image1,cmap1)
  addImageToWorld(world,image2,cmap2)
  makeFrame(world,name)

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0,perc=100):
  ipg = ImagePanelGroup(image)
  ipg.setColorModel(cmap)
  #ipg.setSlices(k1,k2,k3)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def addColorBar(frame,label):
  cbar = ColorBar(label)
  cbar.setFont(cbar.getFont().deriveFont(64.0))
  frame.add(cbar,BorderLayout.EAST)
  #frame.viewCanvas.setBackground(frame.getBackground())
  return cbar

def makeFrame(world,name=None):
  frame = SimpleFrame(world)
  frame.setBackground(Color.WHITE)
  if name:
    frame.setTitle(name)
  view = frame.getOrbitView()
  #zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  #view.setAxesScale(1.0,1.0,zscale)
  view.setScale(2.0)
  #view.setAzimuth(azimuth)
  #view.setElevation(elevation)
  #view.setWorldSphere(BoundingSphere(BoundingBox(f3-1.0,f2,f1,l3,l2,l1)))
  frame.viewCanvas.setBackground(frame.getBackground())
  #frame.setSize(1020,750)
  frame.setSize(1000,800)
  frame.setVisible(True)
  return frame

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

#############################################################################
# Run the function main on the Swing thread
import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())
