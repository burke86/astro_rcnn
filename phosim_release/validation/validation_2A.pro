;;
;; @package phosim
;; @file validation_2A.pro
;; @brief validation task 2A
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

pro validation_2A,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 2A'
  !p.multi=[0,4,3]

  undersample=10.0
  rebinfactor=10
  mrebin=1
  rebin1=4000./float(rebinfactor)
  rebin2=4000./float(rebinfactor)
  measurerebin=4000./float(rebinfactor*mrebin)
  plotscale=1000.0


  avgoverlap=0.0
  for ii=0,4 do begin

  if ii eq 0 then begin & fname1='lsst_e_2000_f2_R42_S11_E000.fits.gz' & fname2='spot_540.txt' & wv=0.54 & endif
  if ii eq 1 then begin & fname1='lsst_e_2001_f2_R42_S11_E000.fits.gz' & fname2='spot_580.txt' & wv=0.58 & endif
  if ii eq 2 then begin & fname1='lsst_e_2002_f2_R42_S11_E000.fits.gz' & fname2='spot_620.txt' & wv=0.62 & endif
  if ii eq 3 then begin & fname1='lsst_e_2003_f2_R42_S11_E000.fits.gz' & fname2='spot_660.txt' & wv=0.66 & endif
  if ii eq 4 then begin & fname1='lsst_e_2004_f2_R42_S11_E000.fits.gz' & fname2='spot_700.txt' & wv=0.70 & endif


     readcol,fname2,a,b,c,d,e,/silent
;     rrr=sqrt(e*e+d*d)
;     phi=atan(d,e)
;     ggg=where(rrr lt 4150 or (rrr ge 4150 and (abs(phi) gt 0.2 or abs(phi) lt 0.5)))
;     a=a(ggg) & b=b(ggg) & c=c(ggg) & d=d(ggg) & e=e(ggg)
     imagea=hist_2d(c*1000.0,b*1000.0,min1=253911.0-2000*0.01,max1=253911.0+1999*0.01,bin1=0.01*undersample,min2=-2000*0.01,max2=1999*0.01,bin2=0.01*undersample)
     sssa=size(imagea)
     xxa=(findgen(sssa(1))-sssa(1)/2.0)*0.01*undersample
     yya=(findgen(sssa(2))-sssa(2)/2.0)*0.01*undersample

     contour,alog((imagea*plotscale)+1.0),xxa,yya,/fill,xtitle='X - 253911.0 (microns)',ytitle='Y (microns)',nlevels=11,/xstyle,/ystyle

     measurepsf,rebin(imagea,measurerebin,measurerebin),rms,e1,e2,medx,medy,flux
     xyouts,-17.0,17.0,'FWHM: '+string(rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
     xyouts,-17.0,14.0,'!9S!3e!D1!N*FWHM: '+string(e1/abs(e1)*sqrt(abs(e1))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
     xyouts,-17.0,11.0,'!9S!3e!D2!N*FWHM: '+string(e2/abs(e2)*sqrt(abs(e2))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
     xyouts,-17.0,8.0,'x: '+string(medx*0.01*undersample*mrebin+253911.0,format='(F10.2)')+' !4l!3m'
     xyouts,-17.0,5.0,'y: '+string(medy*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'

     ss='ZEMAX !4k!3='+string(wv,format='(F4.2)')+' !4l!3m'
     xyouts,-18.0,-18.0,ss,/data

     data=mrdfits(fname1,0,/silent)
     data=data(0:3999,36:4035)
     image=data/total(data)*18841.0*undersample*undersample
     sss=size(image)
     xx=(findgen(sss(1))-sss(1)/2.0)*0.01
     yy=(findgen(sss(2))-sss(2)/2.0)*0.01

     imagec=rebin(data,rebin1,rebin2)
     imagec=imagec/total(imagec)*18841.0*undersample*undersample/float(rebinfactor)/float(rebinfactor)
     sssc=size(imagec)
     xxc=(findgen(sssc(1))-sssc(1)/2.0)*0.01*rebinfactor
     yyc=(findgen(sssc(2))-sssc(2)/2.0)*0.01*rebinfactor


     contour,alog((imagec*plotscale)+1.0),xxc,yyc,/fill,xtitle='X - 253911.0 (microns)',ytitle='Y (microns)',nlevels=11,/xstyle,/ystyle


     measurepsf,rebin(imagec,measurerebin,measurerebin),rms,e1,e2,medx,medy,flux
     xyouts,-17.0,17.0,'FWHM: '+string(rms*2.35*0.01*undersample*mrebin,format='(F4.2)')+' !4l!3m'
     xyouts,-17.0,14.0,'!9S!3e!D1!N*FWHM: '+string(e1/abs(e1)*sqrt(abs(e1))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
     xyouts,-17.0,11.0,'!9S!3e!D2!N*FWHM: '+string(e2/abs(e2)*sqrt(abs(e2))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
     xyouts,-17.0,8.0,'x: '+string(medx*0.01*undersample*mrebin+253911.0,format='(F10.2)')+' !4l!3m'
     xyouts,-17.0,5.0,'y: '+string(medy*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'


     overlap=0.0 & notoverlap=0.0
     for jj=0L,sssc(1)-1 do begin
        for kk=0L,sssc(2)-1 do begin
           if (imagec(jj,kk) gt 10.0) then begin
              jn=round(jj/undersample*float(rebinfactor))
              kn=round(kk/undersample*float(rebinfactor))
              if (imagea(jn,kn) gt 0.0) then overlap=overlap+1.0 else notoverlap=notoverlap+1.0
           endif
        endfor
     endfor
     overlap=overlap/(overlap+notoverlap)
     avgoverlap=avgoverlap+overlap

     ss='phoSim !4k!3='+string(wv,format='(F4.2)')+' !4l!3m'
     xyouts,-18.0,-18.0,ss,/data

  endfor
  avgoverlap=avgoverlap/5.0

  readcol,'spot_540.txt',a,b,c,/silent
  btotal=b & ctotal=c
  readcol,'spot_580.txt',a,b,c,/silent
  btotal=[b,btotal] & ctotal=[c,ctotal]
  readcol,'spot_620.txt',a,b,c,/silent
  btotal=[b,btotal] & ctotal=[c,ctotal]
  readcol,'spot_660.txt',a,b,c,/silent
  btotal=[b,btotal] & ctotal=[c,ctotal]
  readcol,'spot_700.txt',a,b,c,/silent
  btotal=[b,btotal] & ctotal=[c,ctotal]
  imagea=hist_2d(ctotal*1000.0,btotal*1000.0,min1=253911.0-2000*0.01,max1=253911.0+1999*0.01,bin1=0.01*undersample,min2=-2000*0.01,max2=1999*0.01,bin2=0.01*undersample)
  sss=size(imagea)
  xxa=(findgen(sss(1))-sss(1)/2.0)*0.01*undersample
  yya=(findgen(sss(2))-sss(2)/2.0)*0.01*undersample
  contour,alog((imagea*plotscale)+1.0),xxa,yya,/fill,xtitle='X - 253911.0 (microns)',ytitle='Y (microns)',nlevels=11,/xstyle,/ystyle

  measurepsf,rebin(imagea,measurerebin,measurerebin),rms,e1,e2,medx,medy,flux
  xyouts,-17.0,17.0,'FWHM: '+string(rms*2.35*0.01*undersample*mrebin,format='(F4.2)')+' !4l!3m'
  xyouts,-17.0,14.0,'!9S!3e!D1!N*FWHM: '+string(e1/abs(e1)*sqrt(abs(e1))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
  xyouts,-17.0,11.0,'!9S!3e!D2!N*FWHM: '+string(e2/abs(e2)*sqrt(abs(e2))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
  xyouts,-17.0,8.0,'x: '+string(medx*0.01*undersample*mrebin+253911.0,format='(F10.2)')+' !4l!3m'
  xyouts,-17.0,5.0,'y: '+string(medy*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'

  mm1=rms*2.35*0.01*undersample*mrebin
  mm2=e1
  mm3=e2
  mm4=medx*0.01*undersample*mrebin
  mm5=medy*0.01*undersample*mrebin

  xyouts,-18.0,-18.0,'ZEMAX Composite',/data

  data=mrdfits('lsst_e_2000_f2_R42_S11_E000.fits.gz',0,/silent)
  datatotal=18841.0*data/total(data)*undersample*undersample/float(rebinfactor)/float(rebinfactor)
  data=mrdfits('lsst_e_2001_f2_R42_S11_E000.fits.gz',0,/silent)
  datatotal=18841.0*data/total(data)*undersample*undersample/float(rebinfactor)/float(rebinfactor)+datatotal
  data=mrdfits('lsst_e_2002_f2_R42_S11_E000.fits.gz',0,/silent)
  datatotal=18841.0*data/total(data)*undersample*undersample/float(rebinfactor)/float(rebinfactor)+datatotal
  data=mrdfits('lsst_e_2003_f2_R42_S11_E000.fits.gz',0,/silent)
  datatotal=18841.0*data/total(data)*undersample*undersample/float(rebinfactor)/float(rebinfactor)+datatotal
  data=mrdfits('lsst_e_2004_f2_R42_S11_E000.fits.gz',0,/silent)
  datatotal=18841.0*data/total(data)*undersample*undersample/float(rebinfactor)/float(rebinfactor)+datatotal

  datatotal=datatotal(0:3999,36:4035)
  image=datatotal
  sss=size(image)
  xx=(findgen(sss(1))-sss(1)/2.0)*0.01
  yy=(findgen(sss(2))-sss(2)/2.0)*0.01

  imagec=rebin(datatotal,rebin1,rebin2)
  imagec=imagec/total(imagec)*18841.0*undersample*undersample/float(rebinfactor)/float(rebinfactor)
  sssc=size(imagec)
  xxc=(findgen(sssc(1))-sssc(1)/2.0)*0.01*rebinfactor
  yyc=(findgen(sssc(2))-sssc(2)/2.0)*0.01*rebinfactor


  contour,alog((imagec*plotscale)+1.0),xxc,yyc,/fill,xtitle='X - 253911.0 (microns)',ytitle='Y (microns)',nlevels=11,/xstyle,/ystyle

  measurepsf,rebin(imagec,measurerebin,measurerebin),rms,e1,e2,medx,medy,flux
  xyouts,-17.0,17.0,'FWHM: '+string(rms*2.35*0.01*undersample*mrebin,format='(F4.2)')+' !4l!3m'
  xyouts,-17.0,14.0,'!9S!3e!D1!N*FWHM: '+string(e1/abs(e1)*sqrt(abs(e1))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
  xyouts,-17.0,11.0,'!9S!3e!D2!N*FWHM: '+string(e2/abs(e2)*sqrt(abs(e2))*rms*2.35*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'
  xyouts,-17.0,8.0,'x: '+string(medx*0.01*undersample*mrebin+253911.0,format='(F10.2)')+' !4l!3m'
  xyouts,-17.0,5.0,'y: '+string(medy*0.01*undersample*mrebin,format='(F5.2)')+' !4l!3m'

  xyouts,-18.0,-18.0,'PhoSim Composite',/data

  md1=rms*2.35*0.01*undersample*mrebin
  md2=e1
  md3=e2
  md4=medx*0.01*undersample*mrebin
  md5=medy*0.01*undersample*mrebin



  ss='Raytrace PSF & Astrometry'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 2A; '+vers
  xyouts,0.7,0.98,ss,/normal

  value(nnn,0)=avgoverlap*100.0
  value(nnn,1)=abs(md1-mm1)
  value(nnn,2)=sqrt(sqrt((md2-mm2)^2+(md3-mm3)^2))*abs(md1-mm1)
  value(nnn,3)=sqrt((md4-mm4)^2+(md5-mm5)^2)

  name(nnn,0)='Spot Diagram Overlap'
  name(nnn,1)='PSF FWHM error'
  name(nnn,2)='!9S!3ellip*FWHM error'
  name(nnn,3)='Centroid error'

  tolerance_low(nnn,0)=95.0
  tolerance_high(nnn,0)=1e30
  tolerance_high(nnn,1)=0.25
  tolerance_high(nnn,2)=0.25
  tolerance_high(nnn,3)=0.25

  task(nnn,0)='2A Spot Diagrams'

  unit(nnn,0)='% overlap'
  comparison(nnn,0)='ZEMAX Code'

  unit(nnn,1)=' microns'
  comparison(nnn,1)='ZEMAX Code'

  unit(nnn,2)=' microns'
  comparison(nnn,2)='ZEMAX Code'

  unit(nnn,3)=' microns'
  comparison(nnn,3)='ZEMAX Code'

END
