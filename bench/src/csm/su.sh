#!/bin/sh

DATADIR=/data/seis/csm/fc2012/line141s10
#suaddhead <$DATADIR/shotsp.dat ntrpr=342 tsort=1 ns=4001 |
#sushw key=dt a=2000 |
#sushw key=tracf a=1 b=1 j=342 |
#sushw key=fldr a=1 c=1 j=342 |
#sushw key=sx a=735.0 c=15.0 j=342 |
#sushw key=gx a=0.000 b=15.0 j=342 |
#suchw key1=cdp key2=fldr key3=tracf b=1 c=1 d=1 |
#suchw key1=offset key2=gx key3=sx b=1 c=-1 \
#>$DATADIR/shotsp.su

#suxmovie n1=4001 n2=342 perc=95 <$DATADIR/shotsp.su

#susort cdp offset <$DATADIR/shotsp.su >$DATADIR/cmps.su

#suwind <$DATADIR/shotsp.su key=offset abs=1 min=0 max=10000 |
#suwind <$DATADIR/shotsp.su key=offset abs=1 min=100 max=10000 |
suwind <$DATADIR/shotsp.su key=fldr min=19 max=14 |
#sumute mode=2 ntaper=50 tmute=0.3 xmute=2.655 absolute=1 linvel=332 |
sufilter f=0,5,35,45 |
sugain tpow=2 |
#sugain agc=1 |
suximage perc=95

#suwind <$DATADIR/cmps.su tmin=0.0 tmax=4.0 |
#suwind key=offset abs=1 min=100 max=10000 |
#sufilter f=0,5,35,45 |
#sugain tpow=2 |
#sunmo tnmo=0,4 vnmo=3500,3500 |
#sugain agc=1 |
#sustack |
#sudipfilt dt=2.0 dx=7.5 \
#  slopes=-1.0,-0.85,0.85,1.0 amps=0,1,1,0 |
#sugain agc=1 |
#suximage

#suwind <$DATADIR/cmps.su tmin=0.0 tmax=3.0 |
#suwind key=offset abs=1 min=100 max=10000 |
#suwind key=cdp min=1001 max=1002 |
#suvelan nv=101 dv=25 fv=3000 |
#suxcontour
