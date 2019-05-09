pro validation_1C_cutout,qqq,nn

  np=16384.0
  platescale=180000.0
  angle0=296250./platescale/7.0/9.0
  pixelsizemicron=10.0*4096.0/np
  buffer=round(np/48.0)
  pixelsize=pixelsizemicron/platescale*3600.0

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

                 if qqq eq 0 and nn lt 10 then filename='lsst_e_120'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 1 and nn lt 10 then filename='lsst_e_140'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 2 and nn lt 10 then filename='lsst_e_150'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 3 and nn lt 10 then filename='lsst_e_180'+string(nn,format='(I1)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                                           '_E000'
                 if qqq eq 0 and nn ge 10 then filename='lsst_e_12'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 1 and nn ge 10 then filename='lsst_e_14'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 2 and nn ge 10 then filename='lsst_e_15'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                          '_E000'
                 if qqq eq 3 and nn ge 10 then filename='lsst_e_18'+string(nn,format='(I2)')+'_f'+string(filter,format='(I1)')+'_R'+$
                          string(rx,format='(I1)')+string(ry,format='(I1)')+'_S'+$
                          string(sx,format='(I1)')+string(sy,format='(I1)')+$
                                           '_E000'

                 filename2=filename+'.fits.gz'
                 data=mrdfits(filename2,0,status=status,/silent)
                 if status ne 0 then begin
                    print,filename2
                    goto,skip
                 endif
                 
                 sq=size(data)
                 if sq(1) ne np or sq(2) ne np then begin
                    print,filename2,sq(1),sq(2)
                    goto,skip
                 endif

                 for a=-4,4 do begin
                    for b=-4,4 do begin
                       anglex=((rx-2)*27+(sx-1)*9+a)*angle0*!Pi/180.
                       angley=((ry-2)*27+(sy-1)*9+b)*angle0*!Pi/180.
                       vcorr=sqrt(1.0+tan(anglex)*tan(anglex)+tan(angley)*tan(angley)/cos(anglex)/cos(anglex))
                       xx=floor((tan(anglex)*platescale/(!Pi/180.) - (((rx-2)*3+(sx-1))*42250.0+((rx-2.0)*250.0)))/pixelsizemicron+np/2)
                       yy=floor((tan(angley)/cos(anglex)*platescale/(!Pi/180.) - (((ry-2)*3+(sy-1))*42250.0+((ry-2.0)*250.0)))/pixelsizemicron+np/2)
                       subdata=data((xx-buffer):(xx+buffer),(yy-buffer):(yy+buffer))
                       if total(subdata) eq 0 then begin
                          print,'Error: ',xx,yy
                       endif
                       
                       filename3=filename+'_'+string(a+4,format='(I1)')+string(b+4,format='(I1)')+'.fits'
                       mwrfits,subdata,filename3
                       sss='gzip -f '+filename3
                       spawn,sss
                       print,filename3,qqq,nn,rx,ry,sx,sy,a,b
                    endfor
                 endfor
                 skip:
              endfor
        endfor


skipraft:

     endfor
  endfor




end

