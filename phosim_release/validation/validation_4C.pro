;;
;; @package phosim
;; @file validation_4C.pro
;; @brief validation task 4C
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

pro validation_4C,nnn,vers,value,tolerance_low,tolerance_high,task,name,unit,comparison

  print,'Task 4C'
  !p.multi=[0,2,1]

  data=mrdfits('lsst_e_4200_f0_R22_S11_E000.fits.gz',0,/silent)
  image=data(2050:2250,1896:2176)
  sss=size(image)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*10.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*10.0
  contour,alog(image>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]
  data1=data
  imagea=image

  data=mrdfits('lsst_e_4201_f0_R22_S11_E000.fits.gz',0,/silent)
  image=data(2050:2250,1896:2176)
  sss=size(image)
  xx=(findgen(sss(1))-(sss(1)-2)/2.0)*10.0
  yy=(findgen(sss(2))-(sss(2)-2)/2.0)*10.0
  contour,alog(image>1),xx,yy,/fill,xtitle='X Position (microns)',ytitle='Y Position (microns)',/xstyle,/ystyle,nlevels=11,levels=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5]
  data2=data

  ss='Bright Star Optimization Accuracy'
  xyouts,0.1,0.98,ss,/normal
  ss='Validation Task 4C; '+vers
  xyouts,0.7,0.98,ss,/normal

  good=where(data1 ge 100 and data2 ge 100)
  avgv=mean((data1(good)-data2(good))/data2(good))
  value(nnn,0)=avgv*100.0
  name(nnn,0)='PSF Comparison'
  tolerance_low(nnn,0)=-50.0
  tolerance_high(nnn,0)=50.0
  unit(nnn,0)=' %'
  comparison(nnn,0)='Exact Calculation'
  task(nnn,0)='4C m=12 Star; Opt on/off'

END
