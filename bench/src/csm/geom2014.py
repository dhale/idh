"""
Script to make SPS (geometry) files for field camp 2014
The three files (sps, rps, and xps) are needed because
we have two receiver lines. One for the main line of
geophones and one for stations 1242-1341 with 3C (DSU) phones.
Author: Dave Hale, Colorado School of Mines
Version: 2014.05.18
"""

import sys

#############################################################################

surveyDir = "/data/seis/csm/fc2014/survey/"

def readGeom():
  f = open(surveyDir+"s2014GeomFinal.txt")
  lines = f.readlines()
  f.close()
  slist,xlist,ylist,zlist = [],[],[],[]
  for line in lines[1:]:
    s,x,y,z = line.split()
    slist.append(s)
    xlist.append(x)
    ylist.append(y)
    zlist.append(z)
  return slist,xlist,ylist,zlist

"""
R2014                ssss G1                  xxxxxxx.x yyyyyyyy.y zzzz.z
"""
def writeR(slist,xlist,ylist,zlist,fileName):
  pformat = "R%4s                %4s G1                  %9.1f%10.1f%6.1f\n"
  f = open(surveyDir+fileName,"w")
  n = len(slist)
  for i in range(n):
    s = slist[i]
    x = xlist[i]
    y = ylist[i]
    z = zlist[i]
    f.write(pformat % ("2013",s,float(x),float(y),float(z)))
    f.write(pformat % ("2014",s,float(x),float(y),float(z)))
    f.write(pformat % ("2015",s,float(x),float(y),float(z)))
  f.close()

"""
S2014                ssss                     xxxxxxx.x yyyyyyyy.y zzzz.z
"""
def writeS(slist,xlist,ylist,zlist,fileName):
  pformat = "S%4s                %4s                     %9.1f%10.1f%6.1f\n"
  f = open(surveyDir+fileName,"w")
  n = len(slist)
  for i in range(n):
    s = slist[i]
    x = xlist[i]
    y = ylist[i]
    z = zlist[i]
    f.write(pformat % ("2014",s,float(x),float(y),float(z)))
  f.close()

"""
x file format
             slin                sstn chnFchnL1rlin                recF    recL1
             2014                1010    1 24012014                1001    12401
             2014                1010    1 29832013                1242    13411
             2014                1010    2 29932013                1242    13411
             2014                1010    3 30032013                1242    13411
             2014                1011    1 24012014                1001    12401
             2014                1011    1 29832013                1242    13411
             2014                1011    2 29932013                1242    13411
             2014                1011    3 30032013                1242    13411
             ...
             2014                1120    1 24012014                1001    12401
             2014                1120    1 29832013                1242    13411
             2014                1120    2 29932013                1242    13411
             2014                1120    3 30032013                1242    13411
             2014                1121    1 24012014                1002    12411
             2014                1121    1 29832013                1242    13411
             2014                1121    2 29932013                1242    13411
             2014                1121    3 30032013                1242    13411
             ...
"""
def writeX(fileName):
  f2014 = "X            2014                %4s %4i%4i12014                %4i    %4i1\n"
  f2013 = "X            2014                %4s %4i%4i32013                %4i    %4i1\n"
  f = open(surveyDir+fileName,"w")
  n = len(slist)
  for sstn in range(1010,1901):
    if sstn<=1120:
      f.write(f2014 % (sstn,1,240,1001,1240)) # rolling on
    elif sstn<=1780:
      f.write(f2014 % (sstn,1,240,sstn-119,sstn+120)) # rolling
    else:
      f.write(f2014 % (sstn,1,240,1661,1900)) # rolling off
    f.write(f2013 % (sstn,1,298,1242,1341))
    f.write(f2013 % (sstn,2,299,1242,1341))
    f.write(f2013 % (sstn,3,300,1242,1341))
  f.close()

slist,xlist,ylist,zlist = readGeom()
writeR(slist,xlist,ylist,zlist,"geom.rps")
writeS(slist,xlist,ylist,zlist,"geom.sps")
writeX("geom.xps")
