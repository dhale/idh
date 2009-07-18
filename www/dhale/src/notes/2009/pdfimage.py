#!/usr/bin/env python
# Author: Dave Hale, Colorado School of Mines, 2009.07.08
try:
  import os, sys
  from LaunchServices import *
  from Quartz import *
  from CoreGraphics import * # Must import CoreGraphics after Quartz!?
except ImportError:
  print "Sorry, this program is for Macs only."
  sys.exit(-1)

if len(sys.argv)<2:
  print """
  Converts PDF files to JPEG, PNG, or TIFF (CMYK) images for figures in
  printed manuscripts. Specifically, each page of the input PDF file is 
  rendered into an image with resolution 720 dpi. Width and height (in
  pixels) of each output image equals twice the width and height (in 
  points) of the corresponding input PDF page. Output PNG and JPEG 
  images have RGB color; TIFF images have CMYK color. To save space,
  JPEG quality is low; for higher quality, use PNG or TIFF.

  usage: pdfimage foo.pdf [jpg] [png] [tif] (default is all three)
  """
  sys.exit(1)

# Conversion is performed in two steps. First, we use python bindings
# to Quartz to convert a PDF file to an image of the appropriate type.
# However, the function writeToFile for bitmap contexts does not let
# us specify image parameters, such as DPI or compression. Therefore,
# in the second step, we use the MAC OS X command-line program sips to 
# set those parameters in the image files.

def jpgWrite(bmc,ipage):
  jpgFile = pdfBase+suffix(ipage)+".jpg"
  bmc.writeToFile(jpgFile,kCGImageFormatJPEG)
  os.system( "sips" +
    " -s formatOptions low" +
    " -s dpiWidth 720 " +
    " -s dpiHeight 720 " +
    " "+jpgFile)

def pngWrite(bmc,ipage):
  pngFile = pdfBase+suffix(ipage)+".png"
  bmc.writeToFile(pngFile,kCGImageFormatPNG)
  os.system( "sips" +
    " -s dpiWidth 720 " +
    " -s dpiHeight 720 " +
    " "+pngFile)

def tifWrite(bmc,ipage):
  tifFile = pdfBase+suffix(ipage)+".tif"
  bmc.writeToFile(tifFile,kCGImageFormatTIFF)
  cmykProfile ="\"/System/Library/ColorSync/Profiles/Generic CMYK Profile.icc\""
  os.system( "sips" +
    " -s formatOptions lzw" +
    " -s dpiWidth 720 " +
    " -s dpiHeight 720 " +
    " -m "+cmykProfile+" " +
    " "+tifFile)

def suffix(ipage):
  if npage<=1:
    return ""
  else:
    return "_"+str(ipage)

pdfFile = sys.argv[1]
pdfBase,ext = os.path.splitext(pdfFile)
jpg,png,tif = False,False,False
for arg in sys.argv[1:]:
  if arg=="jpg": jpg = True
  if arg=="png": png = True
  if arg=="tif": tif = True
if not (jpg or png or tif):
  jpg = png = tif = True

cs = CGColorSpaceCreateWithName(kCGColorSpaceGenericRGB)
pdf = CGPDFDocumentCreateWithProvider(
        CGDataProviderCreateWithFilename(pdfFile))
npage = pdf.getNumberOfPages()
for ipage in range(1,npage+1):
  page = pdf.getPage(ipage)
  rect = page.getBoxRect(kCGPDFBleedBox)
  width = int(rect.getWidth())
  height = int(rect.getHeight())
  scale = 144.0/72.0 # render PDF with twice its 72 points/inch resolution
  width = int(width*scale+0.5)
  height = int(height*scale+0.5)
  bmc = CGBitmapContextCreateWithColor(width,height,cs,(1,1,1,1))
  bmc.scaleCTM(scale,scale)
  bmc.drawPDFDocument(rect,pdf,ipage)
  if jpg: jpgWrite(bmc,ipage)
  if png: pngWrite(bmc,ipage)
  if tif: tifWrite(bmc,ipage)
