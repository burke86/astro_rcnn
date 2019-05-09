pro validation_1C_measure,qqq,nn


  np=16384.0
  platescale=180000.0
  angle0=296250./platescale/7.0/9.0
  pixelsizemicron=10.0*4096.0/np
  buffer=round(np/48.0)
  pixelsize=pixelsizemicron/platescale*3600.0

  if nn lt 10 then ff='measure_'+string(qqq,format='(I1)')+'_'+string(nn,format='(I1)')+'.txt'
  if nn ge 10 then ff='measure_'+string(qqq,format='(I1)')+'_'+string(nn,format='(I2)')+'.txt'
  openw,1,ff
  for rx=0,4 do begin
     for ry=0,4 do begin
        if (rx eq 0 and ry eq 0) or $
           (rx eq 4 and ry eq 0) or $
           (rx eq 0 and ry eq 4) or $
           (rx eq 4 and ry eq 4) then goto,skipraft
        for sx=0,2 do begin
           for sy=0,2 do begin

              if qqq eq 0 then filter=2
              if qqq gt 0 and (nn eq 0 or nn eq 1 or nn eq 2 or nn eq 3) then filter=0
              if qqq gt 0 and (nn eq 4 or nn eq 5 or nn eq 6 or nn eq 7) then filter=2
              if qqq gt 0 and (nn eq 8 or nn eq 9 or nn eq 10 or nn eq 11) then filter=5

                 for a=-4,4 do begin
                    for b=-4,4 do begin

                 if qqq eq 0 and nn lt 10 then filename='lsst_e_120'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 1 and nn lt 10 then filename='lsst_e_140'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 2 and nn lt 10 then filename='lsst_e_150'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 3 and nn lt 10 then filename='lsst_e_180'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                                                        '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'

                 if qqq eq 0 and nn ge 10 then filename='lsst_e_12'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 1 and nn ge 10 then filename='lsst_e_14'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 2 and nn ge 10 then filename='lsst_e_15'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'
                 if qqq eq 3 and nn ge 10 then filename='lsst_e_18'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits.gz'


                 subdata=mrdfits(filename,0,status=status,/silent)
                 if status ne 0 then goto,skip
                 sss='rm '+filename
                 spawn,sss
                 if total(subdata) eq 0 then goto,skip

                       anglex=((rx-2)*27+(sx-1)*9+a)*angle0*!Pi/180.
                       angley=((ry-2)*27+(sy-1)*9+b)*angle0*!Pi/180.
                       measurepsf,subdata,rms,e1,e2,medx,medy,flux
                       if finite(rms) ne 1 or finite(e1) ne 1 or finite(e2) ne 1 or $
                          finite(medx) ne 1 or finite(medy) ne 1 or finite(flux) ne 1 then begin
                          print,filename,rms,e1,e2,medx,medy,flux
                          goto,finished
                       endif
                       print,qqq,nn,rx,ry,sx,sy,a,b,rms,e1,e2,medx,medy,flux
                       printf,1,rx,ry,sx,sy,a,b,rms,e1,e2,medx,medy,flux,format='(I3,I3,I3,I3,I3,I3,F,F,F,F,F,F)'
                    endfor
                 endfor
                 skip:
              endfor
        endfor
skipraft:
     endfor
  endfor
  
  close,1
  
        finished:
        

end
