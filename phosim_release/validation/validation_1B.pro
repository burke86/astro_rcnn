;;
;; @package phosim
;; @file validation_1B.pro
;; @brief validation task 1B
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

pro validation_1B,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 1B'
  !p.multi=0

  sarr=fltarr(9)
  e1arr=fltarr(9)
  e2arr=fltarr(9)
  xarr=fltarr(9)
  yarr=fltarr(9)
  for nn=0,8 do begin
     filename='lsst_e_110'+string(nn,format='(I1)')+'_f2_R22_S11_E000.fits.gz'
     data=mrdfits(filename,0,/silent)
     measurepsf,data,rms,e1,e2,medx,medy,flux
     sarr(nn)=rms/50.0*2.35
     e1arr(nn)=e1
     e2arr(nn)=e2
     xarr(nn)=medx/50.0
     yarr(nn)=medy/50.0
  endfor

;grid in m of screen
  gridv=float([65536,16384,4096,1024,256,64,16,4,1])/100.0

  plot,gridv,sarr,psym=-4,ytitle='Value',/xlog,yr=[-0.5,1.2],xr=[1e3,5e-3],/xstyle,/ystyle,xtitle='Screen Pixel Size (m)'
  oplot,gridv,xarr,psym=-4,linestyle=2,color=250
  oplot,gridv,yarr,psym=-4,linestyle=2,color=70
  oplot,gridv,e1arr,psym=-4,linestyle=1,color=40
  oplot,gridv,e2arr,psym=-4,linestyle=1,color=130

  legend,linestyle=[0,2,2,1,1],['PSF FWHM (arcsec)','X Centroid (arcsec)','Y Centroid (arcsec)','Ellipticity (e1)','Ellipticity (e2)'],/top,/left,color=[0,250,70,40,130]

  pass1=sqrt((sarr(8)-sarr(7))^2)
  pass2=sqrt((xarr(8)-xarr(7))^2+(yarr(8)-yarr(7))^2)
  pass3=sqrt(sqrt((e1arr(8)-e1arr(7))^2+(e2arr(8)-e2arr(7))^2))*pass1

  xyouts,0.6,0.5,'FWHM Accuracy',/normal
  xyouts,0.6,0.47,'Centroid Accuracy',/normal
  xyouts,0.6,0.44,'Sqrt(Ellip)*FWHM Accuracy',/normal


  xyouts,0.8,0.5,string(pass1,format='(F6.4)')+' arcsec',/normal
  xyouts,0.8,0.47,string(pass2,format='(F6.4)')+' arcsec',/normal
  xyouts,0.8,0.44,string(pass3,format='(F6.4)')+' arcsec',/normal

  ss='Atmospheric Screen Convergence'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 1B; '+vers
  xyouts,0.7,0.98,ss,/normal


  value(nnn,0)=pass1
  value(nnn,1)=pass2
  value(nnn,2)=pass3
  tolerance_high(nnn,0)=0.01
  tolerance_high(nnn,1)=0.01
  tolerance_high(nnn,2)=0.01

  name(nnn,0)='PSF FWHM Error'
  name(nnn,1)='PSF Centroid Error'
  name(nnn,2)='PSF Sqrt(Ellip)*FWHM Error'
  unit(nnn,0)=' arcsec'
  unit(nnn,1)=' arcsec'
  unit(nnn,2)=' arcsec'
  comparison(nnn,0)='PhoSim double pixel size screens'
  comparison(nnn,1)='PhoSim double pixel size screens'
  comparison(nnn,2)='PhoSim double pixel size screens'
  task(nnn,0)='1B Screen Convergence'

END
