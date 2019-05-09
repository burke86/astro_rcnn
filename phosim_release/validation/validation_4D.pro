;;
;; @package phosim
;; @file validation_4D.pro
;; @brief validation task 4D
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

pro validation_4D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 4D'

  for mode=0,1 do begin
  for chip=0,1 do begin
  !p.multi=[0,2,2]

  if mode eq 0 then begin
     if chip eq 0 then begin
        file1='lsst_e_4300_f2_R10_S00_E000.fits.gz'
        file2='lsst_e_4301_f2_R10_S00_E000.fits.gz'
        file3='lsst_e_4302_f2_R10_S00_E000.fits.gz'
     endif else begin
        file1='lsst_e_4300_f2_R22_S11_E000.fits.gz'
        file2='lsst_e_4301_f2_R22_S11_E000.fits.gz'
        file3='lsst_e_4302_f2_R22_S11_E000.fits.gz'
     endelse
  endif else begin
     if chip eq 0 then begin
        file1='lsst_e_4303_f2_R10_S00_E000.fits.gz'
        file2='lsst_e_4304_f2_R10_S00_E000.fits.gz'
        file3='lsst_e_4305_f2_R10_S00_E000.fits.gz'
     endif else begin
        file1='lsst_e_4303_f2_R22_S11_E000.fits.gz'
        file2='lsst_e_4304_f2_R22_S11_E000.fits.gz'
        file3='lsst_e_4305_f2_R22_S11_E000.fits.gz'
     endelse
  endelse



  
  data1=mrdfits(file1,0,/silent)
  image1=rebin(data1(0:3999,0:3999),400,400)

  data2=mrdfits(file2,0,/silent)
  image2=rebin(data2(0:3999,0:3999),400,400)

  data3=mrdfits(file3,0,/silent)
  image3=rebin(data3(0:3999,0:3999),400,400)

  nl=100
  result=moment([image1,image2,image3])
  sig=2.0
  ll=(2.0*sig*sqrt(result(1)))*findgen(nl)/double(nl)+(result(0)-sig*sqrt(result(1)))

  
  sss=size(image1)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image1,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image1=image1*100.0

  sss=size(image2)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image2,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image2=image2*100.0

  sss=size(image3)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image3,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image3=image3*100.0

  q1=histogram(data1,min=0,bin=5)
  q2=histogram(data2,min=0,bin=5)
  q3=histogram(data3,min=0,bin=5)
  v1=moment(data1) & print,v1(0),sqrt(v1(1)),max(data1),total(data1)
  v2=moment(data2) & print,v2(0),sqrt(v2(1)),max(data2),total(data2)
  v3=moment(data3) & print,v3(0),sqrt(v3(1)),max(data3),total(data3)
  x1=findgen(N_elements(q1))*5.+5./2.
  x2=findgen(N_elements(q2))*5.+5./2.
  x3=findgen(N_elements(q3))*5.+5./2.
  plot,x2,q2,psym=10
  oplot,x1,q1,linestyle=2,psym=10,color=50
  oplot,x3,q3,linestyle=1,psym=10,color=250

  xyouts,0.15,0.9,'Normal Optimization',/norm
  xyouts,0.15,0.45,'Quick Optimization',/norm
  xyouts,0.6,0.9,'Single Photon',/norm

  ss='Background Optimiztion Accuracy'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4D; '+vers
  xyouts,0.7,0.98,ss,/normal


        !p.multi=0
  column1=dblarr(400)
  column2=dblarr(400)
  column3=dblarr(400)
  for i=0L,400-1 do column1(i)=mean(image1(i,*))
  for i=0L,400-1 do column2(i)=mean(image2(i,*))
  for i=0L,400-1 do column3(i)=mean(image3(i,*))
  if chip eq 0 then plot,column1,psym=10,xtitle='Relative Flux',linestyle=2,color=50,yr=[0.5*max(column1),1.02*max(column1)],/ystyle
  if chip eq 1 then plot,column1,psym=10,xtitle='Relative Flux',linestyle=2,color=50,yr=[0.9*max(column1),1.02*max(column1)],/ystyle
  oplot,column2,psym=10
  oplot,column3,linestyle=1,psym=10,color=250
  !p.multi=[0,2,2]

  f1=total(image1) & f2=total(image2) & f3=total(image3)
  good=where(image1 ge 1 or image2 ge 1)
  chi2=1e30
  if N_elements(good) gt 1 then begin
     chi2=total((image1(good)/f1-image2(good)/f2)^2/((image1(good)/f1/f1+image2(good)/f2/f2)))/float(N_elements(good))
  endif

  if chip eq 0 then begin
     name(nnn,0)='Normal Opt Relative Pixel'
     value(nnn,0)=chi2
     tolerance_low(nnn,0)=0.0
     tolerance_high(nnn,0)=2.0
     unit(nnn,0)=' (!4V!3!U2!N/dof)'
     comparison(nnn,0)='Exact Calculation'

     name(nnn,1)='Normal Opt Absolute Pixel'
     value(nnn,1)=abs(f1-f2)/f2*100.0
     tolerance_low(nnn,1)=0.0
     tolerance_high(nnn,1)=10.0
     unit(nnn,1)=' %'
     comparison(nnn,1)='Exact Calculation'

     good=where(image3 ge 1 or image2 ge 1)
     chi2=1e30
     if N_elements(good) gt 1 then begin
        chi2=total((image3(good)/f3-image2(good)/f2)^2/((image3(good)/f3/f3+image2(good)/f2/f2)))/float(N_elements(good))
     endif
     name(nnn,2)='Quick Opt Relative Pixel'
     value(nnn,2)=chi2
     tolerance_low(nnn,2)=0.0
     tolerance_high(nnn,2)=3.0
     unit(nnn,2)=' (!4V!3!U2!N/dof)'
     comparison(nnn,2)='Exact Calculation'

     name(nnn,3)='Quick Opt Absolute Pixel'
     value(nnn,3)=abs(f3-f2)/f2*100.0
     tolerance_low(nnn,3)=0.0
     tolerance_high(nnn,3)=10.0
     unit(nnn,3)=' %'
     comparison(nnn,3)='Exact Calculation'

     task(nnn,0)='4D Corner chip Back Opt on/off'
  endif

  loadct,39

  data1=mrdfits(file1,0,/silent)
  image1=rebin(data1(0:3999,0:3999),400,400)
  image1=image1/total(image1)

  data2=mrdfits(file2,0,/silent)
  image2=rebin(data2(0:3999,0:3999),400,400)
  image2=image2/total(image2)

    data3=mrdfits(file3,0,/silent)
  image3=rebin(data3(0:3999,0:3999),400,400)
  image3=image3/total(image3)

  nl=100
  result=moment([image1,image2,image3])
  sig=2.0
  ll=(2.0*sig*sqrt(result(1)))*findgen(nl)/double(nl)+(result(0)-sig*sqrt(result(1)))

  sss=size(image1)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image1,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image1=image1*100.0

  sss=size(image2)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image2,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image2=image2*100.0

 sss=size(image3)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*100.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*100.0
  contour,image3,xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=nl,levels=ll
  image3=image3*100.0

  
loadct,39
  bb=2e-2
  ss=size(data1)
  data1=data1+randomu(seed,ss(1),ss(2))
  data2=data2+randomu(seed,ss(1),ss(2))
  data3=data3+randomu(seed,ss(1),ss(2))
  q1=histogram(data1/mean(data1),min=0,bin=bb)
  q2=histogram(data2/mean(data2),min=0,bin=bb)
  q3=histogram(data3/mean(data3),min=0,bin=bb)
  v1=moment(data1/mean(data1)) & print,v1(0),sqrt(v1(1)),max(data1/mean(data1)),total(data1/mean(data1))
  v2=moment(data2/mean(data2)) & print,v2(0),sqrt(v2(1)),max(data2/mean(data2)),total(data2/mean(data2))
  v3=moment(data3/mean(data3)) & print,v3(0),sqrt(v3(1)),max(data3/mean(data3)),total(data3/mean(data3))
  x1=findgen(N_elements(q1))*bb+bb/2.
  x2=findgen(N_elements(q2))*bb+bb/2.
  x3=findgen(N_elements(q3))*bb+bb/2.
  plot,x2,q2,psym=10,xtitle='Relative Flux'
  oplot,x1,q1,linestyle=2,psym=10,color=50
  oplot,x3,q3,linestyle=1,psym=10,color=250


  xyouts,0.15,0.9,'Normal Optimization',/norm
  xyouts,0.15,0.45,'Quick Optimization',/norm
  xyouts,0.6,0.9,'Single Photon',/norm

  ss='Background Optimiztion Accuracy'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4D; '+vers
  xyouts,0.7,0.98,ss,/normal

  !p.multi=0
  for i=0L,400-1 do column1(i)=mean(image1(i,*))/mean(image1)
  for i=0L,400-1 do column2(i)=mean(image2(i,*))/mean(image2)
  for i=0L,400-1 do column3(i)=mean(image3(i,*))/mean(image3)
  if chip eq 0 then plot,column1,psym=10,xtitle='Relative Flux',linestyle=0,color=50,yr=[0.6,1.3],/ystyle
  if chip eq 1 then plot,column1,psym=10,xtitle='Relative Flux',linestyle=0,color=50,yr=[0.8,1.1],/ystyle
  oplot,column2,psym=10
  oplot,column3,linestyle=0,psym=10,color=250

endfor
endfor
  
END
