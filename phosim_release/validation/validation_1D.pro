;;
;; @package phosim
;; @file validation_1D.pro
;; @brief validation task 1D
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

pro validation_1D,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 1D'
  !p.multi=[0,1,2]

  xpos=fltarr(10,4072)
  for i=0L,9 do begin
     filename='lsst_e_130'+string(i,format='(I1)')+'_f2_R22_S11_E000.fits.gz'
     data=mrdfits(filename,0,/silent)
     for j=0L,4071 do begin
        center=1999.5
        x=double(0.0)
        n=double(0.0)
        for k=1990,2010 do begin
           x=x+data(k,j)*double(k)
           n=n+data(k,j)
        endfor
        xpos(i,j)=x/(n>1)-center
     endfor
  endfor
  plot,[0,4071*0.2/10.0],[-0.2,0.2],/nodata,xtitle='Time (seconds)',ytitle='1D centroid offset (arcseconds)',/xstyle
  timev=findgen(4072)*0.2/10.
  for i=0L,9 do oplot,timev,xpos(i,*)*0.2,color=255.0*i/10.0

  shifts=fltarr(10,4)
  for i=0L,9 do begin
     good1=where(timev gt 5.0 and timev le 20.0)
     good2=where(timev gt 20.0 and timev le 35.0)
     good3=where(timev gt 35.0 and timev le 50.0)
     good4=where(timev gt 50.0 and timev le 65.0)
     good5=where(timev gt 65.0 and timev le 80.0)
     shifts(i,0)=abs(mean(xpos(i,good1))-mean(xpos(i,good2)))*0.2
     shifts(i,1)=abs(mean(xpos(i,good2))-mean(xpos(i,good3)))*0.2
     shifts(i,2)=abs(mean(xpos(i,good3))-mean(xpos(i,good4)))*0.2
     shifts(i,3)=abs(mean(xpos(i,good4))-mean(xpos(i,good5)))*0.2
  endfor
  q=histogram(shifts,min=0,bin=0.01,max=0.1)
  xx=findgen(N_elements(q))*0.01

  plot,xx,q/total(q),psym=10,xtitle='Shifts between 15 seconds intervals',ytitle='Fraction/bin',linestyle=1,yr=[0,0.5],thick=3

  lbt=[0.007,0.007,0.007,0.007,0.007,0.011,0.011,0.011,0.016,0.016,0.016,0.016,0.021,0.027,0.027,0.027,0.05]
  q2=histogram(lbt,min=0,bin=0.01,max=0.1)
  xx=findgen(N_elements(q))*0.01
  oplot,xx,q2/total(q2),psym=10,color=250,linestyle=2,thick=3

  legend,linestyle=[1,2],['phoSim','LBT Data'],color=[0,250],/top,/right

  imshift=moment(shifts)
  lbshift=moment(lbt)
  ss='phoSim: '+string(imshift(0),format='(F6.3)')+' +/- '+string(sqrt(imshift(1)),format='(F6.3)')+' arcseconds'
  xyouts,0.6,0.3,ss,/normal
  ss='LBT:   '+string(lbshift(0),format='(F6.3)')+' +/- '+string(sqrt(lbshift(1)),format='(F6.3)')+' arcseconds'
  xyouts,0.6,0.27,ss,/normal

  ss='Atmospheric Astrometry'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 1D; '+vers
  xyouts,0.7,0.98,ss,/normal
  ss='Trails from 10 different atmospheres with 0.67 arcsecond seeing'
  xyouts,0.48,0.59,ss,/normal
  kstwo,shifts,lbt,D,prob
  ss='KS p-value from same distribution:   '+string(prob,format='(F6.3)')
  xyouts,0.6,0.24,ss,/normal

;value(nnn,0)=(prob)*100.0

  value(nnn,0)=imshift(0)
  tolerance_low(nnn,0)=lbshift(0)*0.0
  tolerance_high(nnn,0)=lbshift(0)*2.0

  value(nnn,1)=sqrt(imshift(1))


  lratio=1e30 & hratio=0.0
  for ratio=0.0,10.0,0.01 do begin
     kstwo,shifts*ratio,lbt,D,prob
     if prob gt 0.01 and prob lt 1.0 then begin
        if ratio lt lratio then lratio=ratio
        if ratio gt hratio then hratio=ratio
     endif
  endfor
  ss='Ratio:   '+string(imshift(0)/lbshift(0),format='(F5.2)')+' +/- '+string(0.5*(hratio-lratio),format='(F5.2)')
  xyouts,0.6,0.21,ss,/normal

  name(nnn,0)='Avg 1-D astrometric shift betw 15s exp'
  name(nnn,1)='Stdev 1-D astrometric shift betw 15s exp'

  tolerance_low(nnn,1)=sqrt(lbshift(1))*0.0
  tolerance_high(nnn,1)=sqrt(lbshift(1))*2.0

  task(nnn,0)='1D Star Trails'

  unit(nnn,0)=' arcsec'
  comparison(nnn,0)='LBT Trailed data (factor of 2)'

  unit(nnn,1)=' arcsec'
  comparison(nnn,1)='LBT Trailed data (factor of 2)'

END
