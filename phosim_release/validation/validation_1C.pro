;;
;; @package phosim
;; @file validation_1C.pro
;; @brief validation task 1C
;;
;; @brief Created by:
;; @author John R. Peterson (Purdue)
;;
;; @brief Modified by:
;;
;; @warning This code is not fully validated
;; and not ready for full release.  Please
;; treat results with caution.
;;

pro validation_1C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  np=16384.0
  platescale=180000.0
  angle0=296250./platescale/7.0/9.0
  pixelsizemicron=10.0*4096.0/np
  buffer=round(np/48.0)
  pixelsize=pixelsizemicron/platescale*3600.0

  maxshift=0.0
  ccolor=[60,60,60,60,150,150,150,150,250,250,250,250]
  lline=[0,1,2,3,0,1,2,3,0,1,2,3]

  for qqq=0,3 do begin

  if qqq eq 0 then ss2='Validation Task 1C; '+vers
  if qqq eq 1 then ss2='Validation Task 1E; '+vers
  if qqq eq 2 then ss2='Validation Task 1F; '+vers
  if qqq eq 3 then ss2='Validation Task 1I; '+vers


  ssigma=dblarr(12,135,135)
  se1=dblarr(12,135,135)
  se2=dblarr(12,135,135)
  smag=dblarr(12,135,135)
  sdx=dblarr(12,135,135)
  sdy=dblarr(12,135,135)

  ussigma=dblarr(12,135,135)
  use1=dblarr(12,135,135)
  use2=dblarr(12,135,135)
  usmag=dblarr(12,135,135)
  usdx=dblarr(12,135,135)
  usdy=dblarr(12,135,135)

  tussigma=dblarr(12,135,135)
  tuse1=dblarr(12,135,135)
  tuse2=dblarr(12,135,135)
  tusmag=dblarr(12,135,135)
  tusdx=dblarr(12,135,135)
  tusdy=dblarr(12,135,135)

  mussigma=dblarr(12)
  muse1=dblarr(12)
  muse2=dblarr(12)
  musmag=dblarr(12)
  musdx=dblarr(12)
  musdy=dblarr(12)

  px=dblarr(135,135)
  py=dblarr(135,135)
  xxx=(findgen(135)-67.)*angle0
  yyy=(findgen(135)-67.)*angle0
  xx1=(findgen(15)-7.)*angle0*9.0
  yy1=(findgen(15)-7.)*angle0*9.0

  for nn=0,11 do begin
     if nn lt 10 then ff='measure_'+string(qqq,format='(I1)')+'_'+string(nn,format='(I1)')+'.txt'
     if nn ge 10 then ff='measure_'+string(qqq,format='(I1)')+'_'+string(nn,format='(I2)')+'.txt'
     sss='/bin/ls -l '+ff+' > temp'
     spawn,sss
     sss='wc temp > temp2'
     spawn,sss
     readcol,'temp2',aaa,/silent
     if aaa(0) eq 0 then goto,skiprun
     readcol,'temp',aa,ab,ac,ad,ae,format='(A,I,A,A,F)',/silent
     if ae(0) eq 0 then goto,skiprun
     readcol,ff,vrx,vry,vsx,vsy,va,vb,vrms,ve1,ve2,vmedx,vmedy,vflux,/silent
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

                       anglex=((rx-2)*15+(sx-1)*5+a)*angle0*!Pi/180.
                       angley=((ry-2)*15+(sy-1)*5+b)*angle0*!Pi/180.
                       g=where(rx eq vrx and ry eq vry and sx eq vsx and sy eq vsy and $
                               a eq va and b eq vb)
                       if N_elements(g) gt 1 then begin
                          print,'multiple image'
                          goto,skip
                       endif
                       if g(0) eq -1 then begin
                          print,'missing image:  ',qqq,nn,rx,ry,sx,sy,a,b
                          goto,skip
                       endif
                       rms=vrms(g(0))
                       e1=ve1(g(0))
                       e2=ve2(g(0))
                       medx=vmedx(g(0))
                       medy=vmedy(g(0))
                       flux=vflux(g(0))
                       if finite(rms) ne 1 or finite(e1) ne 1 or finite(e2) ne 1 or $
                          finite(medx) ne 1 or finite(medy) ne 1 or finite(flux) ne 1 then begin
                          print,filename,rms,e1,e2,medx,medy,flux
                          goto,finished
                       endif

                       xx=((tan(anglex)*platescale/(!Pi/180.) - (((rx-2)*3+(sx-1))*42250.0+((rx-2.0)*250.0)))/pixelsizemicron)+np/2
                       yy=((tan(angley)/cos(anglex)*platescale/(!Pi/180.) - (((ry-2)*3+(sy-1))*42250.0+((ry-2.0)*250.0)))/pixelsizemicron)+np/2
                       dxx=(xx-floor(xx))-0.5
                       dyy=(yy-floor(yy))-0.5
                       medx=medx-dxx
                       medy=medy-dyy
                       ix=rx*27+sx*9+(a+4)
                       iy=ry*27+sy*9+(b+4)
                       px(ix,iy)=tan(anglex)/pixelsizemicron*pixelsize/(!Pi/180.)
                       py(ix,iy)=tan(angley)/cos(anglex)/pixelsizemicron*pixelsize/(!Pi/180.)
                       ssigma(nn,ix,iy)=rms*2.35*pixelsize    ; FWHM arcsecond
                       se1(nn,ix,iy)=e1
                       se2(nn,ix,iy)=e2
                       sdx(nn,ix,iy)=medx*pixelsize            ; arcsecond
                       sdy(nn,ix,iy)=medy*pixelsize            ; arcsecond
                       smag(nn,ix,iy)=-2.5*alog10(flux>1)        ; magnitudes
                       if qqq eq 0 and (nn le 7) then magshift=16.0
                       if qqq eq 0 and (nn gt 7) then magshift=16.0
                       if qqq eq 1 or qqq eq 2 then magshift=16.0
                       if qqq eq 3 and (nn le 7) then magshift=15.0
                       if qqq eq 3 and (nn gt 7) then magshift=15.0
                       magshift=16.0-magshift
                       smag(nn,ix,iy)=smag(nn,ix,iy)+magshift
;                       print,nn,rx,ry,sx,sy,a,b,rms,e1,e2,medx,medy,flux,max(abs(sdx)),max(abs(sdy))
                    endfor
                 endfor
                 skip:
              endfor
        endfor


skipraft:

     endfor
  endfor
  skiprun:
endfor

print,max(abs(sdx)),max(abs(sdy))
  sigorig=ssigma

;determine global values
  chipx=dblarr(12,15,15)
  chipy=dblarr(12,15,15)
  zeropoint=dblarr(12,15,15)
  e1common=dblarr(12,15,15)
  e2common=dblarr(12,15,15)
  seeingcommon=dblarr(12,15,15)
  count=dblarr(12,15,15)
 for nn=0,11 do begin
     for ix=0,134 do begin
        for iy=0,134 do begin
           if sigorig(nn,ix,iy) ne 0 then begin
              qx=floor(ix/9.0) & qy=floor(iy/9.0)
              chipx(nn,qx,qy)=chipx(nn,qx,qy)+sdx(nn,ix,iy)
              chipy(nn,qx,qy)=chipy(nn,qx,qy)+sdy(nn,ix,iy)
              zeropoint(nn,qx,qy)=zeropoint(nn,qx,qy)+smag(nn,ix,iy)
              e1common(nn,qx,qy)=e1common(nn,qx,qy)+se1(nn,ix,iy)
              e2common(nn,qx,qy)=e2common(nn,qx,qy)+se2(nn,ix,iy)
              seeingcommon(nn,qx,qy)=seeingcommon(nn,qx,qy)+ssigma(nn,ix,iy)
              count(nn,qx,qy)=count(nn,qx,qy)+1.0
           endif
        endfor
     endfor
     chipx(nn,*,*)=chipx(nn,*,*)/(count(nn,*,*)>1.0)
     chipy(nn,*,*)=chipy(nn,*,*)/(count(nn,*,*)>1.0)
     zeropoint(nn,*,*)=zeropoint(nn,*,*)/(count(nn,*,*)>1.0)
     seeingcommon(nn,*,*)=seeingcommon(nn,*,*)/(count(nn,*,*)>1.0)
     e1common(nn,*,*)=e1common(nn,*,*)/(count(nn,*,*)>1.0)
     e2common(nn,*,*)=e2common(nn,*,*)/(count(nn,*,*)>1.0)
     for ix=0,134 do begin
        for iy=0,134 do begin
           if sigorig(nn,ix,iy) ne 0 then begin
              qx=floor(ix/9.0) & qy=floor(iy/9.0)
              usdx(nn,ix,iy)=sdx(nn,ix,iy)
              usdy(nn,ix,iy)=sdy(nn,ix,iy)
              usmag(nn,ix,iy)=smag(nn,ix,iy)
              use1(nn,ix,iy)=se1(nn,ix,iy)
              use2(nn,ix,iy)=se2(nn,ix,iy)
              ussigma(nn,ix,iy)=ssigma(nn,ix,iy)
              sdx(nn,ix,iy)=sdx(nn,ix,iy)-chipx(nn,qx,qy)
              sdy(nn,ix,iy)=sdy(nn,ix,iy)-chipy(nn,qx,qy)
              smag(nn,ix,iy)=smag(nn,ix,iy)-zeropoint(nn,qx,qy)
              se1(nn,ix,iy)=se1(nn,ix,iy)-e1common(nn,qx,qy)
              se2(nn,ix,iy)=se2(nn,ix,iy)-e2common(nn,qx,qy)
              ssigma(nn,ix,iy)=ssigma(nn,ix,iy)-seeingcommon(nn,qx,qy)
           endif
        endfor
     endfor
  endfor

 ;figure out global means
 for nn=0,11 do begin
    ccx=0L
     for ix=0,134 do begin
        for iy=0,134 do begin
           if sigorig(nn,ix,iy) ne 0 then begin
              musdx(nn)=musdx(nn)+usdx(nn,ix,iy)
              musdy(nn)=musdy(nn)+usdy(nn,ix,iy)
              musmag(nn)=musmag(nn)+usmag(nn,ix,iy)
              muse1(nn)=muse1(nn)+use1(nn,ix,iy)
              muse2(nn)=muse2(nn)+use2(nn,ix,iy)
              mussigma(nn)=mussigma(nn)+ussigma(nn,ix,iy)
              ccx=ccx+1L
           endif
        endfor
     endfor
     musdx(nn)=musdx(nn)/double(ccx>1)
     musdy(nn)=musdy(nn)/double(ccx>1)
     musmag(nn)=musmag(nn)/double(ccx>1)
     muse1(nn)=muse1(nn)/double(ccx>1)
     muse2(nn)=muse2(nn)/double(ccx>1)
     mussigma(nn)=mussigma(nn)/double(ccx>1)
  endfor


  for nn=0,11 do begin
     for ix=0,134 do begin
        for iy=0,134 do begin
           if sigorig(nn,ix,iy) ne 0 then begin
              tusdx(nn,ix,iy)=usdx(nn,ix,iy)-musdx(nn)
              tusdy(nn,ix,iy)=usdy(nn,ix,iy)-musdy(nn)
              tusmag(nn,ix,iy)=usmag(nn,ix,iy)-musmag(nn)
              tuse1(nn,ix,iy)=use1(nn,ix,iy)-muse1(nn)
              tuse2(nn,ix,iy)=use2(nn,ix,iy)-muse2(nn)
              tussigma(nn,ix,iy)=ussigma(nn,ix,iy)-mussigma(nn)
           endif
        endfor
     endfor
  endfor


;figure out ranges
  low=dblarr(6)+1e30
  high=dblarr(6)-1e30
  for nn=0,11 do begin
     for ix=0,134 do begin
        for iy=0,134 do begin
           if tusmag(nn,ix,iy) lt low(0) and sigorig(nn,ix,iy) ne 0 then low(0)=tusmag(nn,ix,iy)
           if tusmag(nn,ix,iy) gt high(0) and sigorig(nn,ix,iy) ne 0 then high(0)=tusmag(nn,ix,iy)
           if usdx(nn,ix,iy) lt low(1) and sigorig(nn,ix,iy) ne 0 then low(1)=usdx(nn,ix,iy)
           if usdx(nn,ix,iy) gt high(1) and sigorig(nn,ix,iy) ne 0 then high(1)=usdx(nn,ix,iy)
           if usdy(nn,ix,iy) lt low(2) and sigorig(nn,ix,iy) ne 0 then low(2)=usdy(nn,ix,iy)
           if usdy(nn,ix,iy) gt high(2) and sigorig(nn,ix,iy) ne 0 then high(2)=usdy(nn,ix,iy)
           if tussigma(nn,ix,iy) lt low(3) and sigorig(nn,ix,iy) ne 0 then low(3)=tussigma(nn,ix,iy)
           if tussigma(nn,ix,iy) gt high(3) and sigorig(nn,ix,iy) ne 0 then high(3)=tussigma(nn,ix,iy)
           if use1(nn,ix,iy) lt low(4) and sigorig(nn,ix,iy) ne 0 then low(4)=use1(nn,ix,iy)
           if use1(nn,ix,iy) gt high(4) and sigorig(nn,ix,iy) ne 0 then high(4)=use1(nn,ix,iy)
           if use2(nn,ix,iy) lt low(5) and sigorig(nn,ix,iy) ne 0 then low(5)=use2(nn,ix,iy)
           if use2(nn,ix,iy) gt high(5) and sigorig(nn,ix,iy) ne 0 then high(5)=use2(nn,ix,iy)
        endfor
     endfor
  endfor
  print,low
  print,high

;contour plots
  loadct,5,/silent
  !p.multi=[0,4,3]
  aaa=findgen(33)/32.*2*!Pi
  usersym,cos(aaa),sin(aaa),/fill
;goto,skipcont

  for ii=0,5 do begin
     ll=low(ii)+(high(ii)-low(ii))*findgen(20)/19.
     if ii eq 0 then for nn=0,11 do contour,tusmag(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 1 then for nn=0,11 do contour,usdx(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 2 then for nn=0,11 do contour,usdy(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 3 then for nn=0,11 do contour,tussigma(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 4 then for nn=0,11 do contour,use1(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 5 then for nn=0,11 do contour,use2(nn,*,*),xxx,yyy,nlevels=20,/fill,/xstyle,/ystyle,levels=ll
     if ii eq 0 then ss='Photometry'
     if ii eq 1 then ss='Astrometry (x Component)'
     if ii eq 2 then ss='Astrometry (y Component)'
     if ii eq 3 then ss='PSF FWHM'
     if ii eq 4 then ss='Ellipticity (First Component)'
     if ii eq 5 then ss='Ellipticity (Second Component)'
     xyouts,0.1,0.98,ss,/normal
     xyouts,0.7,0.98,ss2,/normal
  endfor

skipcont:
  !p.multi=[0,2,2]
  loadct,39,/silent
  ;histograms
  for hist=0,7 do begin
     if hist eq 0 then begin
        nm='Relative Zeropoint (mag)'
        hh=zeropoint-mean(zeropoint(where(count ne 0)))
        hn=count
        nx=15
        ny=15
        mx=max([abs(max(hh(where(hn ne 0)))),abs(min(hh(where(hn ne 0))))])*1.1
        mn=-mx
     endif
     if hist eq 1 then begin
        nm='Relative Photometry (mag)'
        hh=smag
        hn=sigorig
        nx=135
        ny=135
        mx=max([abs(max(hh(where(hn ne 0)))),abs(min(hh(where(hn ne 0))))])*1.1
        mn=-mx
     endif
     if hist eq 2 then begin
        nm='Global Astrometry (arcsecond)'
        hh=sqrt((chipx)^2+(chipy)^2)
        hn=count
        nx=15
        ny=15
        mx=max(hh(where(hn ne 0)))*1.1
        mn=min(hh(where(hn ne 0)))/1.1
     endif
     if hist eq 3 then begin
        nm='Differential Astrometry (arcsecond)'
        hh=sqrt(sdx*sdx+sdy*sdy)
        hn=sigorig
        nx=135
        ny=135
        mx=max(hh(where(hn ne 0)))*1.1
        mn=min(hh(where(hn ne 0)))/1.1
     endif
     if hist eq 4 then begin
        nm='Global PSF FWHM (arcsecond)'
        hh=seeingcommon
        hn=count
        nx=15
        ny=15
        mx=max(hh(where(hn ne 0)))*1.1
        mn=min(hh(where(hn ne 0)))/1.1
     endif
     if hist eq 5 then begin
        nm='Differential PSF FWHM (arcsecond)'
        hh=ssigma
        hn=sigorig
        nx=135
        ny=135
        mx=max([abs(max(hh(where(hn ne 0)))),abs(min(hh(where(hn ne 0))))])*1.1
        mn=-mx
     endif
     if hist eq 6 then begin
        nm='Common Mode Ellipticity'
        hh=sqrt(e1common*e1common+e2common*e2common)
        hn=count
        nx=15
        ny=15
        mx=max(hh(where(hn ne 0)))*1.1
        mn=min(hh(where(hn ne 0)))/1.1
     endif
     if hist eq 7 then begin
        nm='Differential Ellipticity'
        hh=sqrt(se1*se1+se2*se2)
        hn=sigorig
        nx=135
        ny=135
        mx=max(hh(where(hn ne 0)))*1.1
        mn=min(hh(where(hn ne 0)))/1.1
     endif
     bb=(mx-mn)/20.0

     q1=histogram(hh(where(hn ne 0)),min=mn,bin=bb,max=mx)
     xx=findgen(N_elements(q1))*bb+bb/2.0+mn
     mmm=max(q1/total(q1))
     plot,xx,q1/total(q1),xtitle=nm,psym=10,yr=[0,mmm*1.2],/ystyle,ytitle='Fraction/bin',/xstyle
     result=moment(hh(where(hn ne 0)))
     ss=string(result(0),format='(F6.3)')+' +/- '+string(sqrt(result(1)),format='(F6.3)')
     xyouts,mn+(mx-mn)/20.0,mmm*1.1,ss,/data
     xyouts,0.7,0.98,ss2,/normal

     for nn=0,11 do begin
        if nx eq 1 and ny eq 1 then begin
           oplot,[hh(nn,0,0)],[0.05],psym=8,color=ccolor(nn),linestyle=lline(nn)
        endif else begin
        kk=0L
        for ix=0,nx-1 do begin
           for iy=0,ny-1 do begin
              if hn(nn,ix,iy) gt 0.0 then begin
                 if kk eq 0 then sss=[hh(nn,ix,iy)]
                 if kk ne 0 then sss=[sss,hh(nn,ix,iy)]
                 kk=kk+1
              endif
           endfor
        endfor
        if kk ne 0 then begin
           q1=histogram(sss,min=mn,bin=bb,max=mx)
           xx=findgen(N_elements(q1))*bb+bb/2.0+mn
           oplot,xx,q1/total(q1)/12.0,color=ccolor(nn),psym=10,linestyle=lline(nn)
        endif
        endelse
     endfor


  endfor


  !p.multi=[0,2,2]




;    cos2phi=(dx*dx-dy*dy)/r2
;    sin2phi=2.0*dx*dy/r2
;    if (r2 eq 0.0) then begin
;        cos2phi=1.0
;        sin2phi=0.0
;    endif
;    epar1=e1atm(ii,jj)*cos2phi+e2atm(ii,jj)*sin2phi
;    epar2=e1atm(kk,ll)*cos2phi+e2atm(kk,ll)*sin2phi
;    eper1=e2atm(ii,jj)*cos2phi-e1atm(ii,jj)*sin2phi
;    eper2=e2atm(kk,ll)*cos2phi-e1atm(kk,ll)*sin2phi
;    corr=epar1*epar2+eper1*eper2





;STRUCTURE FUNCTIONS
  goto,finished
for aaa=0,1 do begin
  for rrr=0,3 do begin

  strfs=fltarr(12,20)
  strfse=fltarr(12,20)
  strfn=fltarr(12,20)
     print,aaa,rrr
  for nn=0,11 do begin
        for ix=0,134 do begin
              for iy=0,134 do begin
                 if (sigorig(nn,ix,iy) ne 0.0) then begin
                    for jy=0,134 do begin
                       for jx=0,134 do begin
                          if (jy lt iy  and jx lt ix) then begin
                             if (sigorig(nn,jx,jy) ne 0.0) then begin
                                dx=xxx(ix)-xxx(jx)
                                dy=yyy(iy)-yyy(jy)
                                r2=dx*dx+dy*dy
                                if rrr eq 0 and aaa eq 0 then corr=(usmag(nn,ix,iy)-usmag(nn,jx,jy))^2
                                if rrr eq 1 and aaa eq 0 then corr=(usdx(nn,ix,iy)-usdx(nn,jx,jy))^2+(usdy(nn,ix,iy)-usdy(nn,jx,jy))^2
                                if rrr eq 2 and aaa eq 0 then corr=(ussigma(nn,ix,iy)-ussigma(nn,jx,jy))^2
                                if rrr eq 3 and aaa eq 0 then corr=(use1(nn,ix,iy)-use1(nn,jx,jy))^2+(use2(nn,ix,iy)-use2(nn,jx,jy))^2
                                if rrr eq 0 and aaa eq 1 then corr=(smag(nn,ix,iy)-smag(nn,jx,jy))^2
                                if rrr eq 1 and aaa eq 1 then corr=(sdx(nn,ix,iy)-sdx(nn,jx,jy))^2+(sdy(nn,ix,iy)-sdy(nn,jx,jy))^2
                                if rrr eq 2 and aaa eq 1 then corr=(ssigma(nn,ix,iy)-ssigma(nn,jx,jy))^2
                                if rrr eq 3 and aaa eq 1 then corr=(se1(nn,ix,iy)-se1(nn,jx,jy))^2+(se2(nn,ix,iy)-se2(nn,jx,jy))^2
                                xi=round(alog10((sqrt(r2))>0.01)*4.0+8.0)
                                strfs(nn,xi)=strfs(nn,xi)+corr
                                strfse(nn,xi)=strfse(nn,xi)+corr*corr
                                strfn(nn,xi)=strfn(nn,xi)+1.0
                             endif
                          endif
                       endfor
                    endfor
                 endif
              endfor
           endfor
     endfor


     for nn=0,11 do begin
        strfse(nn,*)=sqrt(strfse(nn,*)/(strfn(nn,*)>1)-(strfs(nn,*))^2/(strfn(nn,*)>1)^2)
        strfs(nn,*)=strfs(nn,*)/(strfn(nn,*)>1)
     endfor
     for nn=0,11 do begin
        strfs(nn,*)=sqrt(strfs(nn,*))
        strfse(nn,*)=0.5/strfs(nn,*)*strfse(nn,*)
     endfor

     xx=10.^((-8.0+findgen(20))/4.0)
     strfst=strfs(0,*)+strfs(1,*)+strfs(2,*)+strfs(3,*)+strfs(4,*)+strfs(5,*)+strfs(6,*)+strfs(7,*)+strfs(8,*)+strfs(9,*)+strfs(10,*)+strfs(11,*)
     strfnt=strfn(0,*)+strfn(1,*)+strfn(2,*)+strfn(3,*)+strfn(4,*)+strfn(5,*)+strfn(6,*)+strfn(7,*)+strfn(8,*)+strfn(9,*)+strfn(10,*)+strfn(11,*)
     strfsst=sqrt((strfs(0,*)-strfst/12.0)^2+(strfs(1,*)-strfst/12.0)^2+(strfs(2,*)-strfst/12.0)^2+(strfs(3,*)-strfst/12.0)^2+(strfs(4,*)-strfst/12.0)^2+(strfs(5,*)-strfst/12.0)^2+(strfs(6,*)-strfst/12.0)^2+(strfs(7,*)-strfst/12.0)^2+(strfs(8,*)-strfst/12.0)^2+(strfs(9,*)-strfst/12.0)^2+(strfs(10,*)-strfst/12.0)^2+(strfs(11,*)-strfst/12.0)^2)/sqrt(12.0)

     maxy=0.0

     for nn=0,11 do begin
        for ii=0,19 do begin
           if strfn(nn,ii) ne 0 and strfs(nn,ii) gt maxy then maxy=strfs(nn,ii)
        endfor
     endfor

     gg=where(strfnt gt 0)
     if rrr eq 0 then sss='Photometry Structure Function (mag)'
     if rrr eq 1 then sss='Astrometry Structure Function (arcsec)'
     if rrr eq 2 then sss='PSF Size Structure Function (arcsec)'
     if rrr eq 3 then sss='Ellipticity Structure Function'
     plot,xx(gg),strfst(gg)/12.0,psym=4,xtitle='Angle (degrees)',ytitle=sss,xr=[0.02,6.0],/xstyle,/xlog,yr=[0,maxy*1.1],/ystyle

     
     ;CFHT DATA
     if qqq eq 3 then begin
     if rrr eq 2 then begin
        for ttt=3,35 do begin
           if ttt lt 10 then filename='cfht_10'+string(ttt,format='(I1)')+'_size.txt'
           if ttt ge 10 then filename='cfht_1'+string(ttt,format='(I2)')+'_size.txt'
           readcol,filename,sa,sb,/silent
           oplot,sa,sb,linestyle=1,color=40,thick=1
;           if aaa eq 0 then print,'CFHT ',sa,sb
        endfor
     endif
      if rrr eq 3 then begin
        for ttt=3,35 do begin
           if ttt lt 10 then filename='cfht_10'+string(ttt,format='(I1)')+'_ellip.txt'
           if ttt ge 10 then filename='cfht_1'+string(ttt,format='(I2)')+'_ellip.txt'
           readcol,filename,sa,sb,/silent
           oplot,sa,sb,linestyle=1,color=40,thick=1
        endfor
     endif
   endif


     errplot,xx(gg),strfst(gg)/12.0-strfsst(gg),strfst(gg)/12.0+strfsst(gg),width=0.0


     if rrr eq 2 and aaa eq 0 then begin
        print,xx(gg)
        print,strfst(gg)/12.0
        print,strfsst(gg)
     endif

     ;fit asymptotic function
     xxq=xx(gg)
     yyq=strfst(gg)/12.0
     eeq=strfsst(gg)
     gg=where(eeq ne 0.0)
     xxq=xxq(gg)
     yyq=yyq(gg)
     eeq=eeq(gg)
     minchi=1e30
     for angs=-2.0,1.0,0.05 do begin
        for asym=0.0,max(yyq),max(yyq)/10.0 do begin
           for ind=0.0,2.0,0.1 do begin
              ang=10.^angs
              yfit=asym*(xxq/ang)^ind/(1+(xxq/ang)^ind)
              chi=total(((yyq-yfit)/eeq)^2)
              if chi lt minchi then begin
                 minchi=chi
                 bestang=ang
                 bestasym=asym
                 bestind=ind
              endif
           endfor
           endfor
     endfor
     xxs=10.^(findgen(100)/100.0*3.0-2.0)
     yys=bestasym*(xxs/bestang)^bestind/(1+(xxs/bestang)^bestind)
     oplot,xxs,yys,linestyle=3
     ss1='Angle: '+string(bestang,format='(F7.3)')
     xyouts,1.0,maxy*0.15,ss1,/data
     ss1='Index: '+string(bestind,format='(F7.3)')
     xyouts,1.0,maxy*0.09,ss1,/data
     ss1='Max:  '+string(asym,format='(F7.3)')
     xyouts,1.0,maxy*0.03,ss1,/data


     for nn=0,11 do begin
        gg=where(strfn(nn,*) ne 0)
        if N_elements(gg) gt 1 then oplot,xx(gg),strfs(nn,gg),color=ccolor(nn),linestyle=lline(nn)
     endfor

 
  endfor
endfor



finished:
  endfor

END





  name(nnn,0)='Atm X=1 PSF size'
  name(nnn,1)='  Average Ellipticity'
  name(nnn,2)='  Stdev Ellipticity'
  name(nnn,3)='  Ellipticity Decorrelation Length'
  name(nnn,4)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,5)='  Stdev Diff Astrometry, 3 arcmin'
  name(nnn,6)='Atm X=1.4 PSF size'
  name(nnn,7)='  Average Ellipticity'
  name(nnn,8)='  Stdev Ellipticity'
  name(nnn,9)='  Ellipticity Decorrelation Length'
  name(nnn,10)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,11)='  Stdev Diff Astrometry, 3 arcmin'

  name(nnn,12)='Atm+Opt PSF size'
  name(nnn,13)='  Average Ellipticity'
  name(nnn,14)='  Stdev Ellipticity'
  name(nnn,15)='  Ellipticity Decorrelation Length'
  name(nnn,16)='  Avg Diff Astrometry, 3 arcmin'
  name(nnn,17)='  Stdev Diff Astrometry, 3 arcmin'

  tolerance_low(nnn,0)=0.0
  tolerance_high(nnn,0)=2.0
  tolerance_high(nnn,4)=1e30
  tolerance_high(nnn,5)=6.0
  tolerance_high(nnn,10)=1e30
  tolerance_high(nnn,11)=6.0
  tolerance_high(nnn,16)=1e30
  tolerance_high(nnn,17)=6.0
  tolerance_low(nnn,0)=0.65-0.05
  tolerance_high(nnn,0)=0.65+0.05

  tolerance_low(nnn,6)=0.65-0.05
  tolerance_high(nnn,6)=0.65+0.05

  tolerance_low(nnn,4)=10.0
  tolerance_low(nnn,10)=10.0
  tolerance_low(nnn,16)=10.0

  tolerance_low(nnn,12)=sqrt(0.65^2+0.39^2)-0.05
  tolerance_high(nnn,12)=sqrt(0.65^2+0.39^2)+0.05

  unit(nnn,0)=' arcsec'
  unit(nnn,1)=' '
  unit(nnn,2)=' '
  unit(nnn,3)=' degrees'
  unit(nnn,4)=' milliarcsec'
  unit(nnn,5)=' milliarcsec'

  unit(nnn,6)=' arcsec'
  unit(nnn,7)=' '
  unit(nnn,8)=' '
  unit(nnn,9)=' degrees'
  unit(nnn,10)=' milliarcsec'
  unit(nnn,11)=' milliarcsec'

  unit(nnn,12)=' arcsec'
  unit(nnn,13)=' '
  unit(nnn,14)=' '
  unit(nnn,15)=' degrees'
  unit(nnn,16)=' milliarcsec'
  unit(nnn,17)=' milliarcsec'
  comparison(nnn,0)='Input value'
  comparison(nnn,1)='CFHT (some opt included)'
  comparison(nnn,2)='CFHT (some opt included)'
  comparison(nnn,3)='Estimate (of order degree)'
  comparison(nnn,4)='Monet best case estimate for 15s'
  comparison(nnn,5)='Monet !4r!3 of !4r!3 Estimate'

  comparison(nnn,6)='Input value'
  comparison(nnn,7)='CFHT (some opt included)'
  comparison(nnn,8)='CFHT (some opt included)'
  comparison(nnn,9)='Estimate (of order degree)'
  comparison(nnn,10)='Monet best case estimate for 15s'
  comparison(nnn,11)='Monet !4r!3 of !4r!3 Estimate'

  comparison(nnn,12)='Input value+LSST Design'
  tolerance_high(nnn,13)=0.1
  tolerance_high(nnn,14)=0.1
  comparison(nnn,13)='LSST Design'
  comparison(nnn,14)='LSST Design'
  comparison(nnn,15)='Estimate (should be ~degree)'
  comparison(nnn,16)='Monet best case estimate for 15s'
  comparison(nnn,17)='Monet !4r!3 of !4r!3 Estimate'
  task(nnn,0)='1C Grids of Stars'

END
