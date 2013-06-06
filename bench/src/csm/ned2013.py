"""
Reads and displays NED data in grid float format.
Author: Dave Hale, Colorado School of Mines
Version: 2013.05.17
"""
from imports import *

# Seismic line 1 from Dave's Garmin hand-held GPS
s1StnLatLon = (
  (3201, 37.22639, -107.05465),
  (3210, 37.22683, -107.05379),
  (3220, 37.22729, -107.05278),
  (3230, 37.22775, -107.05177),
  (3240, 37.22822, -107.05089),
  (3250, 37.22871, -107.04986),
  (3260, 37.22919, -107.04890),
  (3270, 37.22966, -107.04790),
  (3280, 37.23009, -107.04695),
  (3290, 37.23054, -107.04601),
  (3300, 37.23101, -107.04504),
  (3310, 37.23151, -107.04407),
  (3320, 37.23199, -107.04310),
  (3330, 37.23247, -107.04211),
  (3340, 37.23298, -107.04115),
  (3350, 37.23346, -107.04015),
  (3360, 37.23396, -107.03918),
  (3365, 37.23416, -107.03872),
  (3370, 37.23437, -107.03820),
  (3380, 37.23487, -107.03718),
  (3390, 37.23526, -107.03624),
  (3400, 37.23566, -107.03529),
  (3410, 37.23606, -107.03417),
  (3415, 37.23628, -107.03366),
  (3420, 37.23647, -107.03318),
  (3430, 37.23687, -107.03217),
  (3440, 37.23733, -107.03117),
  (3450, 37.23778, -107.03021),
  (3460, 37.23812, -107.02924),
  (3470, 37.23850, -107.02811),
  (3480, 37.23874, -107.02711),
  (3490, 37.23935, -107.02611),
  (3495, 37.23956, -107.02564),
  (3500, 37.23986, -107.02527),
  (3510, 37.24027, -107.02438),
  (3520, 37.24078, -107.02334),
  (3530, 37.24112, -107.02264),
  (3540, 37.24138, -107.02149),
  (3550, 37.24138, -107.02005),
  (3560, 37.24102, -107.01891),
  (3570, 37.24107, -107.01786),
  (3580, 37.24113, -107.01682),
  #(3587, 37.24091, -107.01563),
  (3590, 37.24087, -107.01551),
  (3594, 37.24090, -107.01523),
  (3600, 37.24122, -107.01475),
  (3610, 37.24179, -107.01403),
  (3620, 37.24236, -107.01313),
  (3624, 37.24265, -107.01274),
  (3632, 37.24270, -107.01178),
  (3633, 37.24257, -107.01173),
  (3640, 37.24248, -107.01088),
  (3650, 37.24263, -107.00980),
  (3660, 37.24273, -107.00875),
  (3670, 37.24336, -107.00792),
  (3680, 37.24391, -107.00699),
  (3690, 37.24449, -107.00612),
  (3700, 37.24488, -107.00519),
  (3705, 37.24503, -107.00470),
  (3710, 37.24521, -107.00400),
  (3720, 37.24567, -107.00310),
  (3730, 37.24601, -107.00206),
  (3740, 37.24642, -107.00099),
  (3750, 37.24680, -106.99989),
  (3756, 37.24705, -106.99926),
  (3770, 37.24789, -106.99803),
  (3780, 37.24839, -106.99718),
  (3790, 37.24903, -106.99627),
  (3800, 37.2496 , -106.99540),
  (3810, 37.25017, -106.99458),
  (3820, 37.25074, -106.99355),
  (3826, 37.25100, -106.99313),
  (3827, 37.25111, -106.99304),
  (3840, 37.25183, -106.99194),
  (3850, 37.25228, -106.99092),
  (3860, 37.25272, -106.98994),
  (3874, 37.25325, -106.98872),
)

# Seismic line 1 from Paul's iPhone.
s1StnLatLonElePaulsiPhone = (
  (3201, 37.226466, -107.054654, 7410),
  (3210, 37.226679, -107.053858, 7381),
  (3220, 37.227362, -107.052861, 7387),
  (3230, 37.227805, -107.051918, 7362),
  (3240, 37.228217, -107.051046, 7344),
  (3250, 37.228747, -107.049923, 7334),
  (3260, 37.229188, -107.049027, 7321),
  (3270, 37.229699, -107.048050, 7309),
  (3280, 37.230157, -107.047010, 7306),
  (3290, 37.230609, -107.045987, 7294),
  (3300, 37.231129, -107.045005, 7270),
  (3310, 37.231524, -107.044136, 7263),
  (3320, 37.231990, -107.043087, 7267),
  (3330, 37.232466, -107.042080, 7249),
  (3340, 37.232922, -107.041163, 7262),
  (3350, 37.233371, -107.040225, 7244),
  (3360, 37.233880, -107.039116, 7199),
  (3365, 37.234126, -107.038662, 7197),
  (3370, 37.234266, -107.038297, 7178),
  (3380, 37.234761, -107.037229, 7146),
  (3390, 37.235154, -107.036126, 7099),
  (3400, 37.235621, -107.035131, 7118),
  (3410, 37.235993, -107.034204, 7118),
  (3420, 37.236384, -107.033187, 7097),
  (3430, 37.236840, -107.032033, 7121),
  (3440, 37.237291, -107.031097, 7114),
  (3450, 37.237772, -107.030176, 7169),
  (3460, 37.238130, -107.029174, 7148),
  (3470, 37.238528, -107.028123, 7140),
  (3480, 37.238914, -107.027051, 7139),
  (3490, 37.239351, -107.026112, 7248),
  (3495, 37.239557, -107.025715, 7318),
  (3500, 37.239811, -107.025303, 7356),
  (3510, 37.240829, -107.024328, 7356),
  (3520, 37.240829, -107.023342, 7356),
)


#############################################################################
def main(args):
  goS3000Srf()
  #goPlot()

def goS3000Srf():
  n = len(s1StnLatLon)
  snms,lats,lons,xms,yms = [],[],[],[],[]
  for (sn,lat,lon) in s1StnLatLon:
    xm,ym = utmFromLatLon(lat,lon)
    snms.append(sn)
    lats.append(lat)
    lons.append(lon)
    xms.append(xm)
    yms.append(ym)
  ss = Sampling(674,1.0,3201.0)
  ciLat = CubicInterpolator(snms,lats)
  ciLon = CubicInterpolator(snms,lons)
  ns = ss.count
  sns,xs,ys,zs = [],[],[],[]
  for js in range(ns):
    sn = ss.getValue(js)
    lat = ciLat.interpolate(sn)
    lon = ciLon.interpolate(sn)
    x,y = utmFromLatLon(lat,lon)
    z = elevation(lat,lon)
    z -= 15.0 # subtract 15 m to match survey and Paul's iPhone
    sns.append(sn)
    xs.append(x)
    ys.append(y)
    zs.append(z)
    pformat = "%12i%12.2f%12.2f%12.2f"
    print pformat % (sn,x,y,z)
  sp = SimplePlot()
  pv = sp.addPoints(xs,ys)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.POINT)
  pv = sp.addPoints(xms,yms)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.PLUS)
  pv.setMarkColor(Color.RED)
  """
  sp = SimplePlot()
  pv = sp.addPoints(snms,lats)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  sp = SimplePlot()
  pv = sp.addPoints(snms,lons)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  """
  sp = SimplePlot()
  pv = sp.addPoints(sns,zs)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkStyle(PointsView.Mark.POINT)

def goS1Geometry():
  n = len(s1StnLatLon)
  sns,xs,ys,zs = [],[],[],[]
  for (sn,lat,lon) in s1StnLatLon:
    x,y = utmFromLatLon(lat,lon)
    z = elevation(lat,lon)
    z -= 15.0 # subtract 15 m to match survey and Paul's iPhone
    sns.append(sn)
    xs.append(x)
    ys.append(y)
    zs.append(z)
    #print "%12i%12.2f%12.2f%12.2f" % (sn,x,y,ele)
  ss = Sampling(674,1.0,3201.0)
  cix = CubicInterpolator(sns,xs)
  ciy = CubicInterpolator(sns,ys)
  ciz = CubicInterpolator(sns,zs)
  ns = ss.count
  ffid = 39
  pattern = 1
  snlive1 = 3201
  cnlive1 = 1
  for js in range(ns):
    sn = ss.getValue(js)
    #if sn<3220: continue
    x = cix.interpolate(sn)
    y = ciy.interpolate(sn)
    z = ciz.interpolate(sn)
    ffid += 1
    #pformat = "%12i%12.2f%12.2f%12.2f%12d%12d%12d%12d"
    #print pformat % (sn,x,y,z,ffid,pattern,snlive1,cnlive1)
    pformat = "%12i%12.2f%12.2f%12.2f"
    print pformat % (sn,x,y,z)

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
  pv.setClips(2000,2300)
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

"""
Some examples of (lat,lon) => (easting,northing) conversion
(37.22639,-107.05465) => (317718.18,4121964.75)
(37.22639,-107.02330) => (320499.76,4121904.85)
(37.24073,-107.05465) => (317752.72,4123555.88)
(37.24073,-107.02330) => (320533.78,4123495.96)
"""
  
def utmFromLatLon(lat,lon):
  def atanh(x):
    return 0.5*log((1.0+x)/(1.0-x))
  a = 6378.137
  f = 1.0/298.257223563
  n = f/(2.0-f)
  e0 = 500.0
  n0 = 0.0
  k0 = 0.9996
  aa = a/(1.0+n)*(1.0+n*n*(1.0/4.0+n*n/64.0))
  a1 = n*(0.5-n*(2.0/3.0-n*5.0/16.0))
  a2 = n*n*(13.0/48.0-n*3.0/5.0)
  a3 = n*n*n*61.0/240.0
  st = 2.0*sqrt(n)/(1.0+n)
  lon0 = -105.0 # reference longitude for UTM Zone 13
  lon -= lon0
  lat *= PI/180.0
  lon *= PI/180.0
  t = sinh(atanh(sin(lat))-st*atanh(st*sin(lat)))
  ep = atan(t/cos(lon))
  np = atanh(sin(lon)/sqrt(1.0+t*t))
  sx  = a1*cos(2.0*ep)*sinh(2.0*np)
  sx += a2*cos(4.0*ep)*sinh(4.0*np)
  sx += a3*cos(6.0*ep)*sinh(6.0*np)
  sy  = a1*sin(2.0*ep)*cosh(2.0*np)
  sy += a2*sin(4.0*ep)*cosh(4.0*np)
  sy += a3*sin(6.0*ep)*cosh(6.0*np)
  x = e0+k0*aa*(np+sx)
  y = n0+k0*aa*(ep+sy)
  x *= 1000.0
  y *= 1000.0
  return x,y

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
