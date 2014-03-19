package csm;

/**
 * For reading SEG-D files.
 * @author Dave Hale, Colorado School of Mines
 */
public class Segd {

  /** 
   * Returns integer value encoded in three bytes b[k:k+2]. 
   */
  public static int decode3(byte[] b, int k) {
    int b0 = b[k  ];
    int b1 = b[k+1];
    int b2 = b[k+2];
    if (b0<0) b0 += 256;
    if (b1<0) b1 += 256;
    if (b2<0) b2 += 256;
    return b0*65536+b1*256+b2;
  }

  /** 
   * Returns fixed-point value (in [0,99999.99]) encoded in bytes b[k:k+4].
   */
  public static double decode5(byte[] b, int k) {
    int b0 = b[k  ];
    int b1 = b[k+1];
    int b2 = b[k+2]; 
    int b3 = b[k+3];
    int b4 = b[k+4];
    if (b0<0) b0 += 256;
    if (b1<0) b1 += 256;
    if (b2<0) b2 += 256;
    if (b3<0) b3 += 256;
    if (b4<0) b4 += 256;
    return b0*65536.0+b1*256.0+b2+b3/256.0+b4/65536.0;
  }
}

/*
def bcd2(b,k):
  """ Returns binary-coded decimal integer from bytes k,k+1 in b."""
  return (1000*((b[k  ]>>4)&0xf)+100*(b[k  ]&0xf)+
            10*((b[k+1]>>4)&0xf)+  1*(b[k+1]&0xf))

def bin3(b,k):
  """ Returns binary integer from bytes k,k+1,k+2 in b."""
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2] 
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  return (b0<<16)|(b1<<8)|(b2)

def bin5(b,k):
  """ Returns binary integer from bytes k,k+1,...,k+4 in b."""
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2] 
  b3 = b[k+3] 
  b4 = b[k+4] 
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  if b3<0: b3 += 256
  if b4<0: b4 += 256
  return b0*65536.0+b1*256.0+b2+b3/256.0+b4/65536.0

def readSegd(segdFile):
  n1,n2 = 3001,309 # number of samples, number of traces
  gh = zerobyte(32) # general header
  th = zerobyte(20) # trace header
  the = zerobyte(32) # trace header extension
  csh = zerobyte(32) # channel set header
  ais = ArrayInputStream(segdFile,ByteOrder.BIG_ENDIAN)
  ais.readBytes(gh) # general header 1
  fn = bcd2(gh,0) # file number
  ais.readBytes(gh) # general header 2
  ais.readBytes(gh) # general header 3
  sln = bin5(gh,3) # source line number
  spn = bin5(gh,8) # source point number
  print "file =",segdFile
  #print "fn = ",fn," sln =",sln," spn =",spn
  cns = 0 # channel set number for seismic traces
  nct = 0 # total number of channels, including aux channels
  for ics in range(16): # for each channel set header, ...
    ais.readBytes(csh) # read channel set header
    cn = csh[1] # channel set number
    ct = (csh[10]>>4)&0xf # channel type (in high 4 bits)
    nc = bcd2(csh,8) # number of channels
    if nc>0: # if we have channels of this type, ...
      #print "cn =",cn," nc =",nc," ct =",ct
      if ct==1: # if seismic, ...
        cns = cn # remember channel set number for seismic
        ncs = nc # remember number of seismic channels
      nct += nc # count total number of channels
  #print "nct =",nct,"cns =",cns
  ais.skipBytes(1024) # skip extended header
  ais.skipBytes(1024) # skip external header
  rpf = 1
  rpl = 1
  f = None
  for ict in range(nct): # for all channels (including aux channels)
    ais.readBytes(th) # trace header
    cn = th[3] # channel set number
    ic = bcd2(th,4) # channel (trace) number
    ais.readBytes(the) # trace header extension 1
    rln = bin3(the,0) # receiver line number
    rpn = bin3(the,3) # receiver point number
    n1 = bin3(the,7) # number of samples
    #print "ic =",ic," rln =",rln," rpn =",rpn," n1 =",n1
    if ic==1:
      rpf = rpn
    elif ic==n2:
      rpl = rpn
    ais.skipBytes(6*len(the)) # skip trace header extensions 2-7
    if cn==cns: # if seismic channel, ...
      #print "ic =",ic," rln =",rln," rpn =",rpn
      if not f:
        f = zerofloat(n1,n2) # the traces
      ais.readFloats(f[ic-1]) # get the seismic trace
    else:
      ais.skipBytes(4*n1) # skip the auxiliary trace
  ais.close()
  f = mul(1.0e-14,f) # scale values to approximate range [-10,10]
  return sln,spn,rpf,rpl,f
*/
